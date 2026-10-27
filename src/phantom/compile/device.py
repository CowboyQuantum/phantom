"""`compile_device`: layout + material + seed → a versioned, serialisable `Device` (ADR 0001, 0023).

The compile tier runs every draw and every derivation once per device instance:

1. the scenario draw (interface charges, gate misalignment, wiring faults, the noise draw) for the
   seed and severity (`disorder.draw_scenario`);
2. the electrostatics basis on the clean layout's grid with the draw applied (closed form, or the
   FEM hook of `compile.backends`);
3. the programmed operating point mapped through the wiring faults to the voltages the gates carry;
4. the two dots, the sensor dot and the constant-interaction model from the basis
   (`charge.derive_double_dot`), with the sensor gate moved onto the flank of its Coulomb peak;
5. the barriers to the leads and between the dots at the tunnelling level, their WKB transmissions
   from the action integral along the cut (`electrostatics.barrier_action`), the rates f_attempt·T,
   the interdot tunnel coupling estimate t ≈ ħω₀ √T / π, and the barrier gate's exponential scale
   from two nearby evaluations;
6. the valley splitting of each dot from its harmonic radius (`valley`), and the g-factors (`spin`).

The `Device` (`phantom.device`) holds parameters only: the layout (as data), the scenario draw, the dots, the CI
capacitances, the sensor coupling, the rates, the valley and spin draws, and the model versions of
every layer that produced them. The basis grid is not stored: it is recomputed on demand from the
layout and the draw, which is deterministic, so the file stays small and diffable (plan §8 #3). The
runtime objects every layer consumes are built from the `Device` by its `runtime_*` methods and never
call a solver.
"""

from __future__ import annotations

import dataclasses
from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np

import phantom
from phantom.charge import ci as charge_ci
from phantom.charge import derived as charge_derived
from phantom.charge.derived import derive_double_dot
from phantom.compile.backends import Backend, basis_for
from phantom.constants import HBAR2_OVER_2ME_EV_NM2
from phantom.device import (
    DEVICE_FORMAT_VERSION,
    CompiledBarrier,
    CompiledDot,
    Device,
    compiled_sensor,
)
from phantom.device import MODEL_VERSION as device_model_version
from phantom.disorder import scenario as disorder_scenario
from phantom.disorder.scenario import ScenarioParams, draw_scenario
from phantom.electrostatics import closed_form as electrostatics_closed_form
from phantom.electrostatics.closed_form import GateBasis
from phantom.electrostatics.dots import barrier_action
from phantom.layout.geometry import Layout, MaterialClass
from phantom.readout import rates as readout_rates
from phantom.readout.rates import attempt_rate_hz
from phantom.spin import dynamics as spin_dynamics
from phantom.spin.dynamics import SpinParams, draw_spin
from phantom.valley import splitting as valley_splitting
from phantom.valley.splitting import ValleyParams, draw_valley_splitting_uev

MODEL_VERSION = 2
"""2: the sensor set point uses the sensor with the draw's coupling asymmetry (#43)."""


@dataclass(frozen=True)
class CompileParams:
    """What `compile_device` needs besides the layout and the seed."""

    operating_point_v: dict[str, float]
    """Programmed voltage per gate name."""
    temperature_k: float = 0.15
    grid_pitch_nm: float = 2.0
    severity: float = 1.0
    """Scales the scenario draw; 0 is the clean device."""
    wiring_faults: bool = True
    backend: Backend = "closed_form"
    sensor_gate: str = "SG"
    sensor_flank_widths: float = 0.9
    """Where on the sensor peak to sit, in units of the peak width from its centre."""
    barrier_gate: str = "B_M"
    reservoir_gates: tuple[str, str] = ("R_L", "R_R")
    lead_fermi_ev: float = 2e-3
    valley: ValleyParams | None = None
    """Defaults to `ValleyParams.for_material(layout.material)`."""
    spin: SpinParams = field(default_factory=SpinParams)


# --- compile --------------------------------------------------------------------------------------


def _programmed_vector(layout: Layout, op: dict[str, float]) -> np.ndarray:
    missing = [g for g in layout.gate_names if g not in op]
    if missing:
        raise ValueError(f"operating point lacks gates {missing}")
    return np.array([float(op[g]) for g in layout.gate_names])


def _sensor_set_point(layout: Layout, derived, v: np.ndarray, p: CompileParams, basis: GateBasis, beta_factors: tuple[float, float] = (1.0, 1.0)):
    """Move the sensor gate so the sensor sits `sensor_flank_widths` from a peak centre, then re-derive.
    The sensor is the one the runtime reads (`device.compiled_sensor`, with the draw's asymmetry)."""
    sensor = compiled_sensor(derived.sensor_alpha_v_per_v, beta_factors)
    if p.sensor_gate not in layout.gate_names:
        return derived, v
    j = layout.gate_names.index(p.sensor_gate)
    if sensor.alpha_v_per_v[j] <= 0:
        return derived, v
    n0 = derived.ci.mean_occupation(v)
    phi0 = float(sensor.potential(v, n0))
    k = np.rint((phi0 - sensor.phase_v) / sensor.period_v)
    target = sensor.phase_v + k * sensor.period_v + p.sensor_flank_widths * sensor.width_v
    v = v.copy()
    v[j] += (target - phi0) / sensor.alpha_v_per_v[j]
    return derive_double_dot(layout, v, temperature_k=p.temperature_k, basis=basis), v


def _barrier(basis: GateBasis, v: np.ndarray, d, point_nm: tuple[float, float], mass: float, lead_fermi_ev: float) -> tuple[float, float, float, float]:
    """(height, width, transmission, attempt rate) of the barrier from dot `d` to `point_nm`, the
    transmission exp(−2S) from the WKB action along the cut."""
    h, w, action = barrier_action(basis, v, d, point_nm, mass, lead_fermi_ev=lead_fermi_ev)
    return float(h), float(w), float(np.exp(-2.0 * action)), attempt_rate_hz(d.curvature_ev_per_nm2, mass)


def interdot_tunnel_coupling_ev(basis: GateBasis, v: np.ndarray, d1, d2, mass: float, lead_fermi_ev: float = 2e-3) -> tuple[float, float, float, float]:
    """(t, height, width, rate) for the interdot barrier: t ≈ ħω₀ √T / π with T the WKB transmission."""
    h, w, trans, f_attempt = _barrier(basis, v, d1, (d2.x_nm, d2.y_nm), mass, lead_fermi_ev)
    hbar_omega = np.sqrt(2 * HBAR2_OVER_2ME_EV_NM2 * max(d1.curvature_ev_per_nm2, 1e-12) / mass)
    return float(hbar_omega * np.sqrt(trans) / np.pi), h, w, f_attempt * trans


def compile_device(layout: Layout, params: CompileParams, seed: int) -> Device:
    """Run every draw and derivation once; return the `Device` (ADR 0023)."""
    material = layout.material
    mat = MaterialClass.by_name(material)
    draw = draw_scenario(layout, ScenarioParams.for_material(material), seed, params.severity, params.wiring_faults)
    perturbed, basis = basis_for(params.backend, layout, draw, params.grid_pitch_nm)
    v_prog = _programmed_vector(layout, params.operating_point_v)
    v = draw.wiring.apply(v_prog) if draw.wiring is not None else v_prog.copy()
    derived = derive_double_dot(perturbed, v, temperature_k=params.temperature_k, basis=basis)
    beta_factors = tuple(draw.wiring.sensor_beta_factors) if draw.wiring is not None else (1.0, 1.0)
    derived, v = _sensor_set_point(perturbed, derived, v, params, basis, beta_factors)
    d1, d2 = derived.dots
    res = perturbed.reservoir_points_nm()
    r_l, r_r = params.reservoir_gates
    h_l, w_l, t_l, f_l = _barrier(basis, v, d1, res[r_l], mat.effective_mass_me, params.lead_fermi_ev)
    h_r, w_r, t_r, f_r = _barrier(basis, v, d2, res[r_r], mat.effective_mass_me, params.lead_fermi_ev)
    t, h_m, w_m, g_m = interdot_tunnel_coupling_ev(basis, v, d1, d2, mat.effective_mass_me, params.lead_fermi_ev)
    g_l, g_r = f_l * t_l, f_r * t_r
    # the barrier gate's exponential scale: t at ±10 mV on the barrier gate, same dots
    jb = layout.gate_names.index(params.barrier_gate)
    dv = 0.01
    vp, vm = v.copy(), v.copy()
    vp[jb] += dv
    vm[jb] -= dv
    t_p = interdot_tunnel_coupling_ev(basis, vp, d1, d2, mat.effective_mass_me, params.lead_fermi_ev)[0]
    t_m = interdot_tunnel_coupling_ev(basis, vm, d1, d2, mat.effective_mass_me, params.lead_fermi_ev)[0]
    v_scale = float(2 * dv / np.log(max(t_p, 1e-300) / max(t_m, 1e-300))) if t_p > t_m > 0 else 0.05
    valley_params = params.valley or ValleyParams.for_material(material)
    radii = np.array([d1.radius_nm, d2.radius_nm] + ([derived.sensor_dot.radius_nm] if derived.sensor_dot else []))
    valleys = draw_valley_splitting_uev(valley_params, radii, np.random.default_rng(seed + 1))
    spin = draw_spin(params.spin, np.random.default_rng(seed + 2), seed)

    def cdot(d, ev: float) -> CompiledDot:
        return CompiledDot(d.name, float(d.x_nm), float(d.y_nm), float(d.energy_ev), float(d.curvature_ev_per_nm2), float(d.radius_nm), tuple(float(a) for a in d.lever_arms), float(ev))

    params_dict = {**asdict(params), "valley": asdict(valley_params), "spin": asdict(params.spin)}
    provenance = {
        "phantom_version": phantom.__version__,
        "format_version": DEVICE_FORMAT_VERSION,
        "compile_model_version": MODEL_VERSION,
        "model_versions": {
            "electrostatics.closed_form": electrostatics_closed_form.MODEL_VERSION,
            "disorder.scenario": disorder_scenario.MODEL_VERSION,
            "charge.ci": charge_ci.MODEL_VERSION,
            "charge.derived": charge_derived.MODEL_VERSION,
            "readout.rates": readout_rates.MODEL_VERSION,
            "valley.splitting": valley_splitting.MODEL_VERSION,
            "spin.dynamics": spin_dynamics.MODEL_VERSION,
            "device": device_model_version,
        },
    }
    return Device(
        seed=seed,
        material=material,
        layout=asdict(layout),
        scenario=draw.to_dict(),
        programmed_point_v=tuple(float(x) for x in v_prog),
        operating_point_v=tuple(float(x) for x in v),
        gate_names=layout.gate_names,
        dots=(cdot(d1, valleys[0]), cdot(d2, valleys[1])),
        sensor_dot=cdot(derived.sensor_dot, valleys[2]) if derived.sensor_dot else None,
        c_dot_af=tuple(tuple(float(x) for x in row) for row in derived.ci.c_dot_af),
        c_gate_af=tuple(tuple(float(x) for x in row) for row in derived.ci.c_gate_af),
        temperature_k=params.temperature_k,
        sensor_alpha_v_per_v=tuple(float(x) for x in derived.sensor_alpha_v_per_v),
        barriers=(
            CompiledBarrier(("dot1", "lead_left"), float(h_l), float(w_l), g_l),
            CompiledBarrier(("dot2", "lead_right"), float(h_r), float(w_r), g_r),
            CompiledBarrier(("dot1", "dot2"), h_m, w_m, g_m),
        ),
        tunnel_coupling_ev=t,
        barrier_gate=params.barrier_gate,
        barrier_v_scale_v=v_scale,
        g_factors=spin.g_factors,
        params=params_dict,
        provenance=provenance,
    )


def replace_params(params: CompileParams, **changes: Any) -> CompileParams:
    return dataclasses.replace(params, **changes)
