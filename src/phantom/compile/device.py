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
   and their slopes against every gate at that level held fixed (#53, as the sensor's in 4a);
5a. the double dot re-derived on a grid along each barrier gate (`shape_gates`): its lever arms,
   offsets, capacitances, positions, radii, barrier actions and the interdot barrier's top, as changes
   from the operating point (#59, ADR 0030);
5b. the interdot barrier's top ΔU above its tunnelling level, its slopes against every gate and its
   curvature ħω_b (#66), and the screening factor per gate (#70, ADR 0034): Kemble's transmission of
   the top is the barrier's one transmission, so the interdot rate is f_a·T, the tunnel coupling
   t = (ħω₀/π)·√T, the interdot barrier's slopes are −2πσ·dΔU/dV/ħω_b and the barrier gate's e-fold
   scale of t is 2 over its slope; without a top between the dots the WKB estimates stand;
5c. the puddles (#72, ADR 0035): the bare potential along the cut between the dots above the interdot
   tunnelling level, each dot's turning point on it (its puddle radius at the operating point) and
   each dot's gate-induced charge there, from which the runtime scales the dots with their occupancy;
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
from phantom.charge import barrier as charge_barrier
from phantom.charge import ci as charge_ci
from phantom.charge import derived as charge_derived
from phantom.charge import puddle as charge_puddle
from phantom.charge import shape as charge_shape
from phantom.charge.barrier import kemble_transmission
from phantom.charge.derived import capacitance_matrix_af, derive_double_dot
from phantom.charge.puddle import PuddleGrowth
from phantom.charge.shape import ShapeGrid, quantity_layout
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
from phantom.electrostatics.dots import barrier_action, line_cut, refined_dot, tunnelling_level_ev
from phantom.layout.geometry import Layout, MaterialClass
from phantom.readout import rates as readout_rates
from phantom.readout.rates import attempt_rate_hz
from phantom.readout.sensor import WING_WIDTHS, SensorLeads, SensorSETParams
from phantom.spin import dynamics as spin_dynamics
from phantom.spin.dynamics import SpinParams, draw_spin
from phantom.valley import splitting as valley_splitting
from phantom.valley.splitting import ValleyParams, draw_valley_splitting_uev

MODEL_VERSION = 8
"""2: the sensor set point uses the sensor with the draw's coupling asymmetry (#43).
3: the sensor dot's leads, charging energy and coupling; the set point on its broadened flank (ADR 0026).
4: the qubit dots' barriers carry their attempt rate, tunnelling level and gate slopes (#53, ADR 0027).
5: the double dot's shape along the barrier gates, on a grid (`Device.shape`, #59, ADR 0030).
6: the interdot barrier's top, its slopes and curvature, and its top on the shape grid (#66, ADR 0032).
7: the interdot barrier's screening per gate; t, the interdot rate, its slopes and the barrier gate's
scale from Kemble's transmission of the top (#70, ADR 0034).
8: the puddles: the compiled cut between the dots, each dot's turning point on it and its gate-induced
charge at the operating point (`Device.puddle`, #72, ADR 0035)."""

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
    shape_gates: tuple[str, ...] = ("B_L", "B_M", "B_R")
    """The barrier gates along which the dots' shape is compiled on a grid (#59, ADR 0030); gates the
    layout lacks are skipped."""
    shape_span_v: float = 0.2
    """Half-span of each shape grid around the operating point, V."""
    shape_points: int = 17
    """Nodes per shape grid (odd, so the operating point is one)."""
    interdot_barrier_screening: float = 0.4
    """The fraction of the barrier gate's action on the interdot barrier's top, in the bare closed-form
    potential, that reaches it (#70, ADR 0034). A placeholder for the screening by the electrons that
    the solver will add (`calibration: none`); 1 is the bare potential. At 0.4 the reference device's
    honeycomb goes from hexagons to a superdot over about 50 mV of B_M, and t changes by a factor e
    every 23 mV of it."""
    interdot_cross_screening: float = 0.3
    """The action of every other gate (the plungers above all) on the interdot barrier's top, relative
    to the barrier gate's, after `interdot_barrier_screening` (the two multiply). A placeholder
    (`calibration: none`); 1 is the bare potential, in which a plunger moves the interdot barrier
    nearly as much as the barrier gate does (#48, #70). At 0.3 a plunger acts about four times more
    weakly than B_M on the reference device."""
    puddles: bool = True
    """Whether the dots grow with their occupancy (#72, ADR 0035): each dot's capacitances scale as
    √(n/N_ref) with its gate-induced charge (the classical N ∝ R³ of a puddle in a parabolic well) and
    the interdot barrier is seen from the puddles' edges. False keeps the dots at their compiled size."""


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


TOP_SAMPLES = 401
"""Points along the straight cut between the dots for the interdot barrier's top."""
TOP_FIT_NM = 4.0
"""Half-width of the parabola fitted at the barrier's top for its curvature, nm."""


def _bilinear(u: np.ndarray, basis: GateBasis, xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    """The map `u` (ny, nx) at points (xs, ys) by bilinear interpolation."""
    fx = np.clip((xs - basis.x_nm[0]) / (basis.x_nm[1] - basis.x_nm[0]), 0, len(basis.x_nm) - 1.000001)
    fy = np.clip((ys - basis.y_nm[0]) / (basis.y_nm[1] - basis.y_nm[0]), 0, len(basis.y_nm) - 1.000001)
    ix, iy = fx.astype(int), fy.astype(int)
    tx, ty = fx - ix, fy - iy
    return u[iy, ix] * (1 - tx) * (1 - ty) + u[iy, ix + 1] * tx * (1 - ty) + u[iy + 1, ix] * (1 - tx) * ty + u[iy + 1, ix + 1] * tx * ty


def interdot_cut(basis: GateBasis, v: np.ndarray, d1, d2, level_ev: float) -> tuple[np.ndarray, np.ndarray]:
    """(distance from dot 1 (nm), U − `level_ev` (eV)) along the straight cut between the dots at carried
    `v`, sampled bilinearly at `TOP_SAMPLES` points."""
    u = basis.energy_ev(v)
    t = np.linspace(0.0, 1.0, TOP_SAMPLES)
    xs, ys = d1.x_nm + t * (d2.x_nm - d1.x_nm), d1.y_nm + t * (d2.y_nm - d1.y_nm)
    return t * float(np.hypot(d2.x_nm - d1.x_nm, d2.y_nm - d1.y_nm)), _bilinear(u, basis, xs, ys) - level_ev


def puddle_turning_points_nm(s: np.ndarray, cut: np.ndarray) -> tuple[float, float]:
    """Each dot's turning point at the level on the cut (`interdot_cut`), measured from its own centre
    toward the other dot: where U − level first rises through zero from dot 1's side and last falls
    through it on dot 2's, by linear interpolation; the saddle's position when the cut never rises
    above the level (#72)."""
    above = np.flatnonzero(cut > 0)
    k_top = int(np.argmax(cut))
    if not len(above):
        return float(s[k_top]), float(s[-1] - s[k_top])
    k = int(above[0])
    first = s[k] if k == 0 else float(s[k - 1] + (0.0 - cut[k - 1]) / (cut[k] - cut[k - 1]) * (s[k] - s[k - 1]))
    k = int(above[-1])
    last = s[k] if k == len(s) - 1 else float(s[k] + (0.0 - cut[k]) / (cut[k + 1] - cut[k]) * (s[k + 1] - s[k]))
    return float(first), float(s[-1] - last)


def interdot_barrier_top(basis: GateBasis, v: np.ndarray, d1, d2, level_ev: float, mass: float) -> tuple[float, float]:
    """(ΔU, ħω_b) of the barrier between the dots at carried `v` (#66): the highest U on the straight cut
    between them (bilinear, `TOP_SAMPLES` points) above the fixed tunnelling level `level_ev`, negative
    once the barrier falls below it; ħω_b = ħ√(|U''|/m*) from a parabola fitted ±`TOP_FIT_NM` around the
    top, NaN when the top sits at an end of the cut (no barrier between the dots)."""
    s, cut = interdot_cut(basis, v, d1, d2, 0.0)
    k = int(np.argmax(cut))
    near = np.abs(s - s[k]) <= TOP_FIT_NM
    if k in (0, len(cut) - 1) or near.sum() < 5:
        return float(cut[k] - level_ev), float("nan")
    a = np.polyfit(s[near] - s[k], cut[near], 2)[0]
    return float(cut[k] - level_ev), float(np.sqrt(2 * HBAR2_OVER_2ME_EV_NM2 * abs(2 * a) / mass))


SHAPE_JUMP_NM = 12.0
"""A dot moving farther than this between neighbouring nodes is taken as another minimum: the grid ends."""


def _shape_row(layout: Layout, basis: GateBasis, v: np.ndarray, params: CompileParams, mass: float, ends: tuple, levels: tuple[float, float, float]) -> np.ndarray | None:
    """One node of the shape grid at carried voltages `v` (`charge.shape.quantity_layout`), or None
    when the potential does not form two dots there. Positions are refined between grid points
    (`refined_dot`), and the capacitances and actions are taken with the refined dots."""
    try:
        derived = derive_double_dot(layout, v, temperature_k=params.temperature_k, basis=basis)
    except ValueError:
        return None
    u = basis.energy_ev(v)
    d1, d2 = (refined_dot(basis, u, d) for d in derived.dots)
    c = capacitance_matrix_af(d1, d2, layout.stack.epsilon_well)
    starts = (d1, d2, d1)
    points = (ends[0], ends[1], (d2.x_nm, d2.y_nm))
    log_rate = [-2.0 * float(barrier_action(basis, v, d, pt, mass, level_ev=lv)[2]) for d, pt, lv in zip(starts, points, levels, strict=True)]
    return np.concatenate(
        [
            np.stack([d1.lever_arms, d2.lever_arms]).ravel(),
            [basis.offset_at(d1.x_nm, d1.y_nm), basis.offset_at(d2.x_nm, d2.y_nm)],
            [c[0, 0], c[1, 1], -c[0, 1]],
            [d1.x_nm, d2.x_nm, d1.y_nm, d2.y_nm],
            [d1.radius_nm, d2.radius_nm],
            log_rate,
            [interdot_barrier_top(basis, v, d1, d2, levels[2], mass)[0]],
        ]
    )


def _shape_grid(layout: Layout, basis: GateBasis, v: np.ndarray, params: CompileParams, mass: float, ends: tuple, levels: tuple[float, float, float]) -> ShapeGrid | None:
    """The double dot re-derived on `shape_points` nodes along each shape gate, ±`shape_span_v` around
    `v` with the other gates there (#59, ADR 0030). Each gate's grid runs outward from the operating
    point until the potential stops forming two dots, or a dot jumps by more than `SHAPE_JUMP_NM` (it
    merged and another minimum took its place); the changes are taken from the operating point's node."""
    gates = [g for g in params.shape_gates if g in layout.gate_names]
    if not gates:
        return None
    n = int(params.shape_points)
    centre = n // 2
    offsets = np.linspace(-params.shape_span_v, params.shape_span_v, n)
    xy = quantity_layout(len(v))
    idx, nodes, changes, below, above = [], [], [], [], []
    for g in gates:
        j = layout.gate_names.index(g)
        rows: dict[int, np.ndarray] = {}
        ended = {}
        for step in (1, -1):
            k, last = centre, None
            while 0 <= k < n:
                w = v.copy()
                w[j] += offsets[k]
                row = _shape_row(layout, basis, w, params, mass, ends, levels)
                if row is not None and last is not None:
                    jump = np.hypot(row[xy["x_nm"]] - last[xy["x_nm"]], row[xy["y_nm"]] - last[xy["y_nm"]])
                    row = row if float(np.max(jump)) <= SHAPE_JUMP_NM else None
                if row is None:
                    break
                rows[k], last = row, row
                k += step
            ended[step] = 0 <= k < n
        keys = sorted(rows)
        if len(keys) < 2:
            continue
        ref = rows[centre]
        idx.append(j)
        nodes.append(tuple(float(v[j] + offsets[k]) for k in keys))
        changes.append(tuple(tuple(float(x) for x in rows[k] - ref) for k in keys))
        below.append(bool(ended[-1]))
        above.append(bool(ended[1]))
    if not idx:
        return None
    return ShapeGrid(tuple(idx), len(v), tuple(nodes), tuple(changes), tuple(below), tuple(above))


def compile_provenance() -> dict[str, Any]:
    """What a device compiled now records about the code that compiled it: the phantom version, the
    device format and the model versions of every module compile draws on."""
    return {
        "phantom_version": phantom.__version__,
        "format_version": DEVICE_FORMAT_VERSION,
        "compile_model_version": MODEL_VERSION,
        "model_versions": {
            "electrostatics.closed_form": electrostatics_closed_form.MODEL_VERSION,
            "disorder.scenario": disorder_scenario.MODEL_VERSION,
            "charge.ci": charge_ci.MODEL_VERSION,
            "charge.derived": charge_derived.MODEL_VERSION,
            "charge.shape": charge_shape.MODEL_VERSION,
            "charge.barrier": charge_barrier.MODEL_VERSION,
            "charge.puddle": charge_puddle.MODEL_VERSION,
            "readout.rates": readout_rates.MODEL_VERSION,
            "valley.splitting": valley_splitting.MODEL_VERSION,
            "spin.dynamics": spin_dynamics.MODEL_VERSION,
            "device": device_model_version,
        },
    }


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
    shape = _shape_grid(perturbed, basis, v, params, mass, (res[r_l], res[r_r]), tuple(b.tunnelling_level_ev for b in barriers))
    # the interdot barrier's top and its slopes against every gate, the dots held where they are (#66)
    level_m = barriers[2].tunnelling_level_ev
    top, hbar_omega_b = interdot_barrier_top(basis, v, d1, d2, level_m, mass)
    top_slopes = []
    for g in range(len(v)):
        vp, vm = v.copy(), v.copy()
        vp[g] += SENSOR_SLOPE_STEP_V
        vm[g] -= SENSOR_SLOPE_STEP_V
        top_slopes.append((interdot_barrier_top(basis, vp, d1, d2, level_m, mass)[0] - interdot_barrier_top(basis, vm, d1, d2, level_m, mass)[0]) / (2 * SENSOR_SLOPE_STEP_V))
    # the barrier's one transmission (#70, ADR 0034): Kemble's of the screened top. The interdot rate,
    # t, the interdot barrier's slopes and the barrier gate's scale all come from it; the WKB estimates
    # above stand only when the potential has no top between the dots
    screening = tuple(float(params.interdot_barrier_screening * (1.0 if g == params.barrier_gate else params.interdot_cross_screening)) for g in layout.gate_names)
    if np.isfinite(hbar_omega_b):
        transmission = float(kemble_transmission(top, hbar_omega_b))
        hbar_omega_0 = np.sqrt(2 * HBAR2_OVER_2ME_EV_NM2 * max(d1.curvature_ev_per_nm2, 1e-12) / mass)
        t = float(hbar_omega_0 / np.pi * np.sqrt(transmission))
        kemble_slopes = tuple(float(-2 * np.pi * s * ds / hbar_omega_b) for s, ds in zip(screening, top_slopes, strict=True))
        inter = barriers[2]
        barriers = (barriers[0], barriers[1], dataclasses.replace(inter, rate_hz=float(inter.attempt_rate_hz * transmission), log_rate_slopes_per_v=kemble_slopes))
        v_scale = float(2.0 / kemble_slopes[jb]) if kemble_slopes[jb] > 0 else v_scale
    valley_params = params.valley or ValleyParams.for_material(material)
    radii = np.array([d1.radius_nm, d2.radius_nm] + ([sensor.sensor_dot.radius_nm] if sensor.sensor_dot else []))
    valleys = draw_valley_splitting_uev(valley_params, radii, np.random.default_rng(seed + 1))
    spin = draw_spin(params.spin, np.random.default_rng(seed + 2), seed)

    def cdot(d, ev: float) -> CompiledDot:
        return CompiledDot(d.name, float(d.x_nm), float(d.y_nm), float(d.energy_ev), float(d.curvature_ev_per_nm2), float(d.radius_nm), tuple(float(a) for a in d.lever_arms), float(ev))

    params_dict = {**asdict(params), "valley": asdict(valley_params), "spin": asdict(params.spin)}
    provenance = compile_provenance()
    device = Device(
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
        shape=shape,
        interdot_barrier_top_ev=top,
        interdot_barrier_top_slopes_ev_per_v=tuple(float(x) for x in top_slopes),
        interdot_barrier_hbar_omega_ev=None if not np.isfinite(hbar_omega_b) else hbar_omega_b,
        interdot_barrier_screening=screening,
    )
    # the puddles (#72, ADR 0035): the compiled cut, each dot's turning point on it, and each dot's
    # gate-induced charge at the operating point, taken through the runtime CI model so that the runtime
    # finds exactly these numbers there
    if params.puddles and np.isfinite(hbar_omega_b):
        s_cut, cut = interdot_cut(basis, v, d1, d2, level_m)
        n_ref = device.runtime_ci().induced_charge(v)
        if np.all(n_ref > 0):
            device = dataclasses.replace(
                device,
                puddle=PuddleGrowth(
                    n_ref=(float(n_ref[0]), float(n_ref[1])),
                    turning_nm=puddle_turning_points_nm(s_cut, cut),
                    cut_nm=tuple(float(x) for x in s_cut),
                    cut_ev=tuple(float(x) for x in cut),
                ),
            )
    return device


def replace_params(params: CompileParams, **changes: Any) -> CompileParams:
    return dataclasses.replace(params, **changes)
