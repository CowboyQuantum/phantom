"""`compile_device`: layout + material + seed → a versioned, serialisable `Device` (ADR 0001, 0023).

The compile tier runs every draw and every derivation once per device instance:

1. the scenario draw (interface charges, gate misalignment, wiring faults, the noise draw) for the
   seed and severity (`disorder.draw_scenario`);
2. the electrostatics basis on the clean layout's grid with the draw applied (closed form, or the
   FEM hook of `compile.backends`);
3. the programmed operating point mapped through the wiring faults to the voltages the gates carry;
4. the two dots, the sensor dot and the constant-interaction model from the basis
   (`charge.derive_double_dot`, with the sensor's charging energy and coupling to the dots), with the
   sensor gate moved onto the flank of its Coulomb peak;
4a. the sensor's leads (ADR 0026): the rate through each sensor barrier to its bath from the WKB action,
   its slope against every gate at the compiled tunnelling level, and each bath's lever arms;
5. the barriers to the leads and between the dots at the tunnelling level, their WKB transmissions
   from the action integral along the cut (`electrostatics.barrier_action`), the rates f_attempt·T
   and their slopes against every gate at that level held fixed (#53, as the sensor's in 4a),
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
from phantom.electrostatics.dots import barrier_action, line_cut, tunnelling_level_ev
from phantom.layout.geometry import Layout, MaterialClass
from phantom.readout import rates as readout_rates
from phantom.readout.rates import attempt_rate_hz
from phantom.readout.sensor import WING_WIDTHS, SensorLeads, SensorSETParams
from phantom.spin import dynamics as spin_dynamics
from phantom.spin.dynamics import SpinParams, draw_spin
from phantom.valley import splitting as valley_splitting
from phantom.valley.splitting import ValleyParams, draw_valley_splitting_uev

MODEL_VERSION = 4
"""2: the sensor set point uses the sensor with the draw's coupling asymmetry (#43).
3: the sensor dot's leads, charging energy and coupling; the set point on its broadened flank (ADR 0026).
4: the qubit dots' barriers carry their attempt rate, tunnelling level and gate slopes (#53, ADR 0027)."""

SENSOR_SLOPE_STEP_V = 0.01
"""Half-step of the central differences that give the rates' slopes against each gate (the sensor's
leads and the qubit dots' barriers), V."""
SET_POINT_WIDTH_PASSES = 3
"""Passes of the last set-point step, each taking the peak width at the previous pass's target."""


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
    sensor_barrier_gates: tuple[str, str] = ("SB_L", "SB_R")
    """The sensor's barrier gates, source side first (ADR 0026)."""
    sensor_reservoir_gates: tuple[str, str] = ("SR_L", "SR_R")
    """The sensor's bath gates, source side first; a layout without them gives the typed sensor."""
    bath_threshold_ev: float = -0.05
    """Electron energy below which a bath is accumulated, eV (the turn-on protocols' placeholder)."""
    bath_width_ev: float = 0.004
    """Width of a bath's accumulation step, eV."""
    valley: ValleyParams | None = None
    """Defaults to `ValleyParams.for_material(layout.material)`."""
    spin: SpinParams = field(default_factory=SpinParams)


# --- compile --------------------------------------------------------------------------------------


def _programmed_vector(layout: Layout, op: dict[str, float]) -> np.ndarray:
    missing = [g for g in layout.gate_names if g not in op]
    if missing:
        raise ValueError(f"operating point lacks gates {missing}")
    return np.array([float(op[g]) for g in layout.gate_names])


def sensor_tunnelling_levels_ev(basis: GateBasis, layout: Layout, v: np.ndarray, s, p: CompileParams, mass: float) -> tuple[float, float] | None:
    """The level each sensor lead's electrons tunnel at, from `tunnelling_level_ev` on the cut from
    the sensor dot to the centroid of that lead's reservoir gate at carried voltages `v`."""
    names = tuple(p.sensor_reservoir_gates)
    if s is None or any(n not in layout.gate_names for n in names):
        return None
    levels = []
    for n in names:
        point = layout.gate(n).centroid_nm
        _, u = line_cut(basis, v, s.x_nm, s.y_nm, point[0], point[1], 401)
        levels.append(float(tunnelling_level_ev(s, float(u[-1]), mass, p.lead_fermi_ev)))
    return levels[0], levels[1]


def _log_rate_slopes(basis: GateBasis, v: np.ndarray, d, point_nm: tuple[float, float], mass: float, level_ev: float) -> tuple[float, ...]:
    """d ln Γ/dV = −2 dS/dV of the barrier from dot `d` to `point_nm` against every gate, at the fixed
    tunnelling level `level_ev`: central differences of ±`SENSOR_SLOPE_STEP_V` with the dot held where
    it is, 1/V."""
    row = np.zeros(len(v))
    for g in range(len(v)):
        vp, vm = v.copy(), v.copy()
        vp[g] += SENSOR_SLOPE_STEP_V
        vm[g] -= SENSOR_SLOPE_STEP_V
        s_p = barrier_action(basis, vp, d, point_nm, mass, level_ev=level_ev)[2]
        s_m = barrier_action(basis, vm, d, point_nm, mass, level_ev=level_ev)[2]
        row[g] = -(s_p - s_m) / SENSOR_SLOPE_STEP_V
    return tuple(float(x) for x in row)


def _sensor_leads(
    basis: GateBasis, layout: Layout, v: np.ndarray, s, p: CompileParams, mass: float, levels: tuple[float, float] | None, slopes: bool = True
) -> SensorLeads | None:
    """The sensor dot's barriers and baths at carried voltages `v` (ADR 0026); None without a sensor
    dot or without the sensor's reservoir gates in the layout.

    Each lead's rate is f_a·exp(−2S), the WKB action S through the barrier from the sensor dot to the
    centroid of its reservoir gate, at the fixed tunnelling level `levels` (compile takes them at the
    wired programmed point, before the sensor set point). On a Coulomb peak the electrons tunnel at
    the leads' electrochemical potential, which no gate moves; the closed-form potential has no
    absolute Fermi level, so the level is placed once by `tunnelling_level_ev` and held there, and a
    gate changes a rate through the barrier's profile alone, at the compile tier and at runtime alike.
    The slopes d ln Γ/dV = −2 dS/dV are central differences of ±`SENSOR_SLOPE_STEP_V` on each gate.
    Each bath's lever arms and offset are the basis and the draw's energy at that centroid."""
    names = tuple(p.sensor_reservoir_gates)
    if s is None or levels is None or any(n not in layout.gate_names for n in names):
        return None
    points = [layout.gate(n).centroid_nm for n in names]
    f_a = attempt_rate_hz(s.curvature_ev_per_nm2, mass)
    heights, widths, rates, rows = [], [], [], []
    for point, level in zip(points, levels, strict=True):
        h, w, action = barrier_action(basis, v, s, point, mass, level_ev=level)
        heights.append(float(h))
        widths.append(float(w))
        rates.append(float(f_a * np.exp(-2.0 * action)))
        rows.append(_log_rate_slopes(basis, v, s, point, mass, level) if slopes else (0.0,) * len(v))
    return SensorLeads(
        barrier_gates=tuple(p.sensor_barrier_gates),
        reservoir_gates=names,
        v_ref_v=tuple(float(x) for x in v),
        gamma_ref_hz=(rates[0], rates[1]),
        attempt_rate_hz=float(f_a),
        log_rate_slopes_per_v=(rows[0], rows[1]),
        bath_lever_arms=tuple(tuple(float(x) for x in basis.at(*pt)) for pt in points),
        bath_offset_ev=tuple(basis.offset_at(*pt) for pt in points),
        bath_threshold_ev=float(p.bath_threshold_ev),
        bath_width_ev=float(p.bath_width_ev),
        tunnelling_level_ev=(float(levels[0]), float(levels[1])),
        barrier_height_ev=(heights[0], heights[1]),
        barrier_width_nm=(widths[0], widths[1]),
    )


def _derived_sensor(derived, beta_factors: tuple[float, float], p: CompileParams, leads: SensorLeads | None) -> SensorSETParams:
    beta = None if derived.sensor_beta_v_per_e is None else tuple(float(x) for x in derived.sensor_beta_v_per_e)
    return compiled_sensor(derived.sensor_alpha_v_per_v, beta_factors, derived.sensor_charging_energy_ev, beta, p.temperature_k, leads)


def _sensor_set_point(
    layout: Layout,
    derived,
    v: np.ndarray,
    p: CompileParams,
    basis: GateBasis,
    beta_factors: tuple[float, float] = (1.0, 1.0),
    mass: float = 0.19,
    levels: tuple[float, float] | None = None,
):
    """Move the sensor gate so the sensor sits `sensor_flank_widths` from a peak centre, then re-derive.
    The sensor is the one the runtime reads (`device.compiled_sensor`, with the draw's asymmetry).
    Returns the derivation, the voltages and the sensor's leads (None for the typed sensor).

    A sensor with leads then takes one more step with the re-derived sensor and its leads compiled
    there, without deriving again: its peak (a fraction of a millivolt of φ_s) is narrow against the
    jump of its potential α·V when a re-derivation moves the dot minimum by a grid step, so only a step
    with the lever arms the runtime keeps lands on the target. Its peak width is that of the sech² with
    the same wing as its broadened peak (`width_v` when ħΓ ≪ kT), at the rates the leads give at the
    step's own target (`SET_POINT_WIDTH_PASSES` fixed passes, the same on every platform). The
    derivation returned is the sensor's, a fraction of a period of the sensor gate from the final
    point, and the leads are anchored there; `compile_device` derives the double dot again at the end."""
    sensor = _derived_sensor(derived, beta_factors, p, _sensor_leads(basis, layout, v, derived.sensor_dot, p, mass, levels, slopes=False))
    if p.sensor_gate not in layout.gate_names:
        return derived, v, None
    j = layout.gate_names.index(p.sensor_gate)
    if sensor.alpha_v_per_v[j] <= 0:
        return derived, v, None

    def step(sensor: SensorSETParams, derived, v: np.ndarray, passes: int) -> np.ndarray:
        phi0 = float(sensor.potential(v, derived.ci.mean_occupation(v)))
        k = np.rint((phi0 - sensor.phase_v) / sensor.period_v)
        out = v
        for _ in range(passes):
            width = sensor.width_v if sensor.leads is None else sensor.wing_offset_v(out) / WING_WIDTHS
            target = sensor.phase_v + k * sensor.period_v + p.sensor_flank_widths * width
            out = v.copy()
            out[j] += (target - phi0) / sensor.alpha_v_per_v[j]
        return out

    v = step(sensor, derived, v, 1)
    derived = derive_double_dot(layout, v, temperature_k=p.temperature_k, basis=basis)
    if sensor.leads is None:
        return derived, v, None
    leads = _sensor_leads(basis, layout, v, derived.sensor_dot, p, mass, levels)
    v = step(_derived_sensor(derived, beta_factors, p, leads), derived, v, SET_POINT_WIDTH_PASSES)
    return derived, v, leads


def _barrier(basis: GateBasis, v: np.ndarray, d, point_nm: tuple[float, float], mass: float, lead_fermi_ev: float) -> tuple[float, float, float, float]:
    """(height, width, transmission, attempt rate) of the barrier from dot `d` to `point_nm`, the
    transmission exp(−2S) from the WKB action along the cut."""
    h, w, action = barrier_action(basis, v, d, point_nm, mass, lead_fermi_ev=lead_fermi_ev)
    return float(h), float(w), float(np.exp(-2.0 * action)), attempt_rate_hz(d.curvature_ev_per_nm2, mass)


def _gated_barrier(basis: GateBasis, v: np.ndarray, d, point_nm: tuple[float, float], mass: float, lead_fermi_ev: float, between: tuple[str, str]) -> CompiledBarrier:
    """The barrier from dot `d` to `point_nm` with what the runtime needs to follow the gates (#53): its
    rate f_a·exp(−2S) at `v`, the attempt rate, the tunnelling level (`tunnelling_level_ev`, placed at
    `v` and then held fixed, as for the sensor's leads, ADR 0026) and the slopes of ln Γ against every
    gate at that level."""
    _, u = line_cut(basis, v, d.x_nm, d.y_nm, point_nm[0], point_nm[1], 401)
    level = float(tunnelling_level_ev(d, float(u[-1]), mass, lead_fermi_ev))
    h, w, trans, f_a = _barrier(basis, v, d, point_nm, mass, lead_fermi_ev)
    return CompiledBarrier(between, float(h), float(w), f_a * trans, float(f_a), level, _log_rate_slopes(basis, v, d, point_nm, mass, level))


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
    levels = sensor_tunnelling_levels_ev(basis, perturbed, v, derived.sensor_dot, params, mat.effective_mass_me)
    derived, v, sensor_leads = _sensor_set_point(perturbed, derived, v, params, basis, beta_factors, mat.effective_mass_me, levels)
    sensor = derived  # the sensor's linearisation: the one its set point was placed with
    if sensor_leads is not None:
        derived = derive_double_dot(perturbed, v, temperature_k=params.temperature_k, basis=basis)
    d1, d2 = derived.dots
    res = perturbed.reservoir_points_nm()
    r_l, r_r = params.reservoir_gates
    mass, fermi = mat.effective_mass_me, params.lead_fermi_ev
    barriers = (
        _gated_barrier(basis, v, d1, res[r_l], mass, fermi, ("dot1", "lead_left")),
        _gated_barrier(basis, v, d2, res[r_r], mass, fermi, ("dot2", "lead_right")),
        _gated_barrier(basis, v, d1, (d2.x_nm, d2.y_nm), mass, fermi, ("dot1", "dot2")),
    )
    t = interdot_tunnel_coupling_ev(basis, v, d1, d2, mass, fermi)[0]
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
    radii = np.array([d1.radius_nm, d2.radius_nm] + ([sensor.sensor_dot.radius_nm] if sensor.sensor_dot else []))
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
        sensor_dot=cdot(sensor.sensor_dot, valleys[2]) if sensor.sensor_dot else None,
        c_dot_af=tuple(tuple(float(x) for x in row) for row in derived.ci.c_dot_af),
        c_gate_af=tuple(tuple(float(x) for x in row) for row in derived.ci.c_gate_af),
        temperature_k=params.temperature_k,
        sensor_alpha_v_per_v=tuple(float(x) for x in sensor.sensor_alpha_v_per_v),
        barriers=barriers,
        tunnel_coupling_ev=t,
        barrier_gate=params.barrier_gate,
        barrier_v_scale_v=v_scale,
        g_factors=spin.g_factors,
        params=params_dict,
        provenance=provenance,
        sensor_charging_energy_ev=sensor.sensor_charging_energy_ev,
        sensor_beta_v_per_e=None if sensor.sensor_beta_v_per_e is None else (float(sensor.sensor_beta_v_per_e[0]), float(sensor.sensor_beta_v_per_e[1])),
        sensor_leads=sensor_leads,
    )


def replace_params(params: CompileParams, **changes: Any) -> CompileParams:
    return dataclasses.replace(params, **changes)
