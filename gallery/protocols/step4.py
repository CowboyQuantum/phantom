"""Step 4 protocols: dynamics from the layout's barriers. WKB lead rates and a master equation give the
charge actually present along a sweep; latching and direction dependence emerge from rate × dwell."""

from __future__ import annotations

import numpy as np

from gallery.protocols.step2 import _axis, _layout, _operating_point
from gallery.runner import Context, Protocol
from phantom.charge import derive_double_dot
from phantom.electrostatics import GateBasis, barrier_action
from phantom.layout import MaterialClass
from phantom.readout import RateEquation, SensorSETParams, wkb_rate_from_action_hz


def derive_dynamics(lay, v, basis, temperature_k: float = 0.15):
    """The derived double dot at `v` plus the WKB rates to the left and right reservoirs and across the
    middle barrier, all read off the landscape at the operating point, from the action integral along
    each barrier (the same estimate as the compiled device, ADR 0023)."""
    d = derive_double_dot(lay, v, temperature_k=temperature_k, basis=basis)
    m = MaterialClass.by_name(lay.material).effective_mass_me
    res = lay.reservoir_points_nm()
    d1, d2 = d.dots
    h_l, w_l, s_l = barrier_action(basis, v, d1, res["R_L"], m)
    h_r, w_r, s_r = barrier_action(basis, v, d2, res["R_R"], m)
    h_m, w_m, s_m = barrier_action(basis, v, d1, (d2.x_nm, d2.y_nm), m)
    g_l = float(wkb_rate_from_action_hz(s_l, d1.curvature_ev_per_nm2, m))
    g_r = float(wkb_rate_from_action_hz(s_r, d2.curvature_ev_per_nm2, m))
    g_m = float(wkb_rate_from_action_hz(s_m, d1.curvature_ev_per_nm2, m))
    return d, RateEquation(d.ci, g_l, g_r, g_m), {"left": (h_l, w_l, g_l), "right": (h_r, w_r, g_r), "interdot": (h_m, w_m, g_m)}


def _sensor_from(d):
    return SensorSETParams(alpha_v_per_v=d.sensor_alpha_v_per_v)


def _map_grid(v0, i1, i2, half, n):
    p1 = v0[i1] + np.linspace(-half, half, n)
    p2 = v0[i2] + np.linspace(-half, half, n)
    vv = np.broadcast_to(v0, (n, n, len(v0))).copy()
    vv[..., i1] = p1[:, None]
    vv[..., i2] = p2[None, :]
    return p1, p2, vv


def _raster_mean_occupation(rate_eq: RateEquation, vv: np.ndarray, dwell_s: float, reverse_fast: bool = False) -> np.ndarray:
    """⟨N⟩ (n1, n2, 2) scanning the fast axis (second index) for each slow value in order; the
    distribution is carried across the flyback so latched charge persists between rows."""
    n1, n2, _ = vv.shape
    out = np.empty((n1, n2, 2))
    p = states = None
    for i in range(n1):
        row = vv[i, ::-1] if reverse_fast else vv[i]
        mean, p, states = rate_eq.sweep(row, dwell_s, p0=p, states0=states)
        out[i] = mean[::-1] if reverse_fast else mean
    return out


def stability_vs_sweep_rate(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """The coarse map at three dwell times per point with the lead rates derived from the layout at the
    operating point: when Γ·dwell ≫ 1 the map is the equilibrium honeycomb; as the dwell shrinks below
    1/Γ the charge lags, lines fade and shift along the fast axis, and finally the charge is frozen."""
    with ctx.compile():
        lay = _layout(p)
        v0 = _operating_point(lay, p)
        basis = GateBasis.from_layout(lay)
        d, rate_eq, bars = derive_dynamics(lay, v0, basis, float(p.instance.get("temperature_k", 0.15)))
        sensor = _sensor_from(d)
        dwells = [float(x) for x in p.instance["dwell_s"]]
    ax1, ax2 = p.sweep
    i1, i2 = lay.gate_names.index(ax1["gate"]), lay.gate_names.index(ax2["gate"])
    n = int(ax1["points"])
    p1, p2, vv = _map_grid(v0, i1, i2, float(ax1["half_span_v"]), n)
    out: dict[str, np.ndarray] = {
        "v_p1_v": p1,
        "v_p2_v": p2,
        "gamma_left_hz": np.array(rate_eq.gamma_left_hz),
        "gamma_right_hz": np.array(rate_eq.gamma_right_hz),
        "gamma_interdot_hz": np.array(rate_eq.gamma_interdot_hz),
        "lead_barrier_mev_nm": np.array([[bars["left"][0] * 1e3, bars["left"][1]], [bars["right"][0] * 1e3, bars["right"][1]]]),
    }
    with ctx.run(n * n * len(dwells)):
        for k, dwell in enumerate(dwells):
            mean = _raster_mean_occupation(rate_eq, vv, dwell)
            phi = sensor.potential(vv, mean)
            out[f"current_dwell{k}_a"] = sensor.current_a(phi)
            out[f"n1_dwell{k}"] = mean[..., 0]
        ctx.count(model_evals=n * n * len(dwells))
    return out


def stability_direction(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """The same map at one intermediate dwell, scanned with the fast axis forward and backward: the
    lagging transitions sit on different sides of their equilibrium position, the hysteresis of a latched
    stability diagram."""
    with ctx.compile():
        lay = _layout(p)
        v0 = _operating_point(lay, p)
        basis = GateBasis.from_layout(lay)
        d, rate_eq, _ = derive_dynamics(lay, v0, basis, float(p.instance.get("temperature_k", 0.15)))
        sensor = _sensor_from(d)
        dwell = float(p.instance["dwell_s"])
    ax1, ax2 = p.sweep
    i1, i2 = lay.gate_names.index(ax1["gate"]), lay.gate_names.index(ax2["gate"])
    n = int(ax1["points"])
    p1, p2, vv = _map_grid(v0, i1, i2, float(ax1["half_span_v"]), n)
    out: dict[str, np.ndarray] = {"v_p1_v": p1, "v_p2_v": p2, "gamma_left_hz": np.array(rate_eq.gamma_left_hz)}
    with ctx.run(2 * n * n):
        fwd = _raster_mean_occupation(rate_eq, vv, dwell)
        bwd = _raster_mean_occupation(rate_eq, vv, dwell, reverse_fast=True)
        out["current_forward_a"] = sensor.current_a(sensor.potential(vv, fwd))
        out["current_backward_a"] = sensor.current_a(sensor.potential(vv, bwd))
        out["n1_difference"] = fwd[..., 0] - bwd[..., 0]
        ctx.count(model_evals=2 * n * n)
    return out


def lead_rate_time_trace(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """P1 stepped across a dot-1 transition at t = 0 for three left-barrier settings; the mean charge on
    dot 1 relaxes exponentially toward its new equilibrium at the lead rate, which the layout predicts
    through the WKB barrier. The relaxation rate read off the trace is the lead rate."""
    with ctx.compile():
        lay = _layout(p)
        v0 = _operating_point(lay, p)
        basis = GateBasis.from_layout(lay)
        settings = [float(x) for x in p.instance["barrier_v"]]
        j_bl = lay.gate_names.index("B_L")
        j_p1 = lay.gate_names.index("P1")
        n_t = int(p.sweep[0]["points"])
        cases = []
        for vb in settings:
            v = v0.copy()
            v[j_bl] = vb
            d, rate_eq, _ = derive_dynamics(lay, v, basis, float(p.instance.get("temperature_k", 0.15)))
            # step P1 by a quarter of a transition period: start just below a transition, end above it
            period = d.ci.transition_period_v(0, "P1")
            cases.append((vb, v, rate_eq, period))
    out: dict[str, np.ndarray] = {}
    with ctx.run(n_t * len(settings)):
        for k, (vb, v, rate_eq, period) in enumerate(cases):
            # place the step across the next dot-1 transition: scan P1 upward for the ground state to gain
            # an electron, then step from 0.2 period below that point to 0.2 period above it
            scan = np.repeat(v[None, :], 400, axis=0)
            scan[:, j_p1] += np.linspace(0.0, 1.5 * period, 400)
            n1 = rate_eq.charge.ground_state_occupation(scan)[:, 0]
            cross = int(np.argmax(n1 > n1[0]))
            v_before = scan[cross].copy()
            v_before[j_p1] -= 0.2 * period
            v_after = scan[cross].copy()
            v_after[j_p1] += 0.2 * period
            gamma = rate_eq.gamma_left_hz
            t = np.linspace(0.0, 6.0 / max(gamma, 1e-3), n_t)
            dt = t[1] - t[0]
            states0, p0 = rate_eq.stationary(v_before)
            seq = np.repeat(v_after[None, :], n_t, axis=0)
            mean, _, _ = rate_eq.sweep(seq, dt, p0=p0, states0=states0)
            n_eq = rate_eq.charge.mean_occupation(v_after)[0]
            out[f"t_case{k}_s"] = t
            out[f"n1_case{k}"] = mean[:, 0]
            out[f"n1_deficit_case{k}"] = np.maximum(n_eq - mean[:, 0], 1e-12)
            out[f"gamma_case{k}_hz"] = np.array(gamma)
            out[f"barrier_case{k}_v"] = np.array(vb)
        ctx.count(model_evals=n_t * len(settings))
    return out


# --- RF reflectometry and single-shot readout ---------------------------------------------------

from phantom.readout import (
    TankCircuit,
    noisy_device_for_sampling,
    sensor_resistance_ohm,
    single_shots,
)
from phantom.readout.measure import SensedDoubleDot as _Sensed


def _tank(p: Protocol) -> TankCircuit:
    return TankCircuit(**p.instance.get("tank", {}))


def rf_resonance(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Reflected amplitude and phase against carrier frequency for several sensor resistances: the
    resonance of the tank, the dip to zero at the matched resistance, and the phase flipping by π as the
    sensor crosses it."""
    with ctx.compile():
        tank = _tank(p)
        resistances = np.array([float(r) for r in p.instance["sensor_resistances_ohm"]])
    f = _axis(p.sweep[0])
    out: dict[str, np.ndarray] = {"f_mhz": f / 1e6, "f0_mhz": np.array(tank.resonance_hz / 1e6), "r_match_ohm": np.array(tank.matched_resistance_ohm)}
    with ctx.run(len(f) * len(resistances)):
        for k, r in enumerate(resistances):
            g = tank.reflection(f, r)
            out[f"mag_r{k}"] = np.abs(g)
            out[f"phase_r{k}_deg"] = np.degrees(np.angle(g))
        ctx.count(model_evals=len(f) * len(resistances))
    return out


def rf_sensor_flank(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """The sensor plunger swept across one Coulomb peak: the DC current, the sensor resistance it implies,
    and the reflected amplitude at the carrier frequency. The RF flank mirrors the DC peak, steepest
    where the resistance crosses the matching value."""
    with ctx.compile():
        d = _Sensed()
        tank = _tank(p)
        f_carrier = float(p.instance.get("carrier_hz", tank.resonance_hz))
        centre = d.sensor.phase_v
    half = float(p.sweep[0].get("half_span_v", 0.15))
    sg = np.linspace(centre - half, centre + half, int(p.sweep[0]["points"]))
    v = np.zeros((len(sg), len(d.charge.gate_names)))
    v[:, 2] = sg
    with ctx.run(len(sg)):
        i = d.current_a(v)
        r = sensor_resistance_ohm(i, d.sensor.bias_v)
        g = tank.reflection(f_carrier, r)
        ctx.count(model_evals=len(sg))
    return {"v_sg_v": sg, "current_a": i, "resistance_ohm": r, "rf_magnitude": np.abs(g), "rf_phase_deg": np.degrees(np.angle(g))}


def _readout_point(d: _Sensed, flank_widths: float):
    v = np.array([0.0, 0.0, d.sensor.phase_v + flank_widths * d.sensor.width_v])
    n0 = d.charge.mean_occupation(v)
    states = np.array([n0, n0 + np.array([1.0, 0.0])])
    return v, states


def single_shot_histogram(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Many consecutive shots at a readout point on the sensor flank, the dot in one of two charge states
    at random: the histogram of integrated signals is bimodal, its separation set by one electron's shift
    of the sensor and its widths by the readout noise over the integration time, with 1/f drift and
    telegraph switching moving both peaks between shots."""
    inst = p.instance
    with ctx.compile():
        base = _Sensed.default_noisy(seed=int(inst.get("seed", 0)))
        v, states = _readout_point(base, float(inst.get("sensor_flank_widths", 0.9)))
        tau = float(inst["tau_s"])
        dt = float(inst.get("dt_s", tau / 20))
        dev = noisy_device_for_sampling(base, dt, seed=int(inst.get("seed", 0)) + 100)
        n_shots = int(inst.get("n_shots", 2000))
        rng = np.random.default_rng(int(inst.get("seed", 0)) + 200)
    with ctx.run(n_shots):
        res = single_shots(dev, v, states, tau, dt, n_shots, rng)
        ctx.count(model_evals=n_shots)
    bins = np.linspace(res.signals_a.min(), res.signals_a.max(), 61)
    centres = 0.5 * (bins[1:] + bins[:-1])
    h0, _ = np.histogram(res.signals_a[res.true_state == 0], bins)
    h1, _ = np.histogram(res.signals_a[res.true_state == 1], bins)
    return {
        "signal_pa": centres * 1e12,
        "counts_state0": h0.astype(float),
        "counts_state1": h1.astype(float),
        "shot_index": np.arange(n_shots, dtype=float),
        "shot_signal_pa": res.signals_a * 1e12,
        "snr": np.array(res.snr),
        "fidelity": np.array(res.fidelity),
    }


def snr_vs_integration_time(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Signal-to-noise ratio and readout fidelity of the two-state discrimination against integration
    time: square-root growth while white noise dominates, then saturation and decline once the shots
    are long enough for 1/f drift and telegraph switching to move the distributions."""
    inst = p.instance
    with ctx.compile():
        base = _Sensed.default_noisy(seed=int(inst.get("seed", 0)))
        v, states = _readout_point(base, float(inst.get("sensor_flank_widths", 0.9)))
        taus = np.array([float(x) for x in inst["tau_s"]])
        n_shots = int(inst.get("n_shots", 400))
    snr = np.empty(len(taus))
    fid = np.empty(len(taus))
    with ctx.run(len(taus) * n_shots):
        for k, tau in enumerate(taus):
            dt = min(tau / 20, 1e-4)
            dev = noisy_device_for_sampling(base, dt, seed=int(inst.get("seed", 0)) + 100 + k)
            res = single_shots(dev, v, states, tau, dt, n_shots, np.random.default_rng(int(inst.get("seed", 0)) + 200 + k))
            snr[k] = res.snr
            fid[k] = res.fidelity
        ctx.count(model_evals=len(taus) * n_shots)
    return {"tau_s": taus, "snr": snr, "fidelity": fid, "sqrt_law": snr[0] * np.sqrt(taus / taus[0])}
