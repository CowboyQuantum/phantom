"""Step 8 protocols: the consumer session on the compiled device — the charge a sample sees, the charge after a hold, and a closing barrier that slows the dynamics until they latch and telegraph (#53)."""

from __future__ import annotations

import numpy as np

from gallery.protocols.step7 import _compiled
from gallery.runner import Context, Protocol
from phantom.device import Device
from phantom.optimum import sensor_optimum


def _tuned_point(dev: Device) -> np.ndarray:
    """The programmed point with the sensor gate on the optimum's flank, as a consumer tunes the sensor
    before measuring: the sensor dot's peaks are a few hundred µV wide, so the programmed sensor gate
    alone usually sits between them (ADR 0026). The gallery device's sensor gate is not faulted, so its
    programmed and carried voltages agree."""
    v = np.array(dev.programmed_point_v)
    v[dev.gate_names.index("SG")] = sensor_optimum(dev).v_sensor_gate_v
    return v


def _dot1_degeneracy(dev: Device, half_v: float = 0.03, n: int = 4001, v0: np.ndarray | None = None) -> np.ndarray:
    """The programmed point along P1 nearest `v0` (the tuned point when None) where dot 1 is half-way
    through a transition to its lead."""
    v0 = _tuned_point(dev) if v0 is None else v0
    j = dev.gate_names.index("P1")
    grid = v0 + np.linspace(-half_v, half_v, n)[:, None] * np.eye(len(v0))[j]
    occ = dev.at(grid).n_mean
    frac = occ[:, 0] - np.floor(occ[:, 0])
    flat2 = np.abs(occ[:, 1] - np.rint(occ[:, 1])) < 1e-6
    return grid[int(np.argmin(np.where(flat2, np.abs(frac - 0.5), np.inf)))]


def session_charge_trace(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """At a dot-1 degeneracy, `Session.sense` in `sampled` mode for dwells of 0.02, 0.5 and 50 times
    1/Γ_L (noise off, so the current follows the charge), and with the device's noise at a 1 ms dwell.
    Short dwells give a telegraph trace correlated from sample to sample; long dwells give independent
    Boltzmann draws."""
    samples = int(p.sweep[0]["points"])
    with ctx.compile():
        dev = _compiled(p)
        v = _dot1_degeneracy(dev)
        gamma = float(dev.at(v).gamma_left_hz)
    out: dict[str, np.ndarray] = {"sample": np.arange(samples, dtype=float)}
    with ctx.run(4 * samples):
        for label, gd in (("fast", 0.02), ("mid", 0.5), ("slow", 50.0)):
            m = dev.session(noise=False).sense(v, gd / gamma, samples=samples)
            out[f"current_{label}_a"] = m.current_a
            out[f"n1_{label}"] = m.n_present[:, 0]
        m = dev.session().sense(v, 1e-3, samples=samples)
        out["current_noisy_1ms_a"] = m.current_a
        out["n1_noisy_1ms"] = m.n_present[:, 0]
        ctx.count(model_evals=4 * samples)
    base = np.floor(out["n1_slow"].min())
    for k, label in enumerate(("fast", "mid", "slow")):
        out[f"n1_{label}_offset"] = out[f"n1_{label}"] - base + 1.5 * k
    out["gamma_left_hz"] = np.array(gamma)
    out["v_degeneracy_v"] = v
    return out


def session_hold_relaxation(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """`hold` at a point A long enough to settle, then `sense` at a point B one electron away along P1:
    the mean charge of dot 1 (noise off) relaxes as n_B + (n_A − n_B)·exp(−Γ_out·t) (drawn beside it),
    Γ_out the total rate out of A's charge state at B: the left lead and the interdot path, which empties
    through dot 2's fast lead (the rates follow the gates, #53). One sampled trajectory steps once at a
    random time."""
    samples = int(p.sweep[0]["points"])
    with ctx.compile():
        dev = _compiled(p)
        v0 = _tuned_point(dev)
        j = dev.gate_names.index("P1")
        per_e = dev.truth().charge_lattice_e_per_v[0][0]
        a = _dot1_degeneracy(dev) + 0.5 / per_e * np.eye(len(v0))[j]
        b = a - 1.0 / per_e * np.eye(len(v0))[j]
        gamma_a = float(dev.at(a).gamma_left_hz)
        vb = dev.effective_voltages(b)
        n_a_state = np.rint(dev.at(a).n_ground).astype(int)
        states = np.unique(np.concatenate([dev.runtime_rates().window_states(vb), n_a_state[None]]), axis=0)
        _, w = dev.runtime_rates().rate_matrix(vb, states)
        i = int(np.flatnonzero((states == n_a_state).all(axis=1))[0])
        gamma = float(-w[i, i])  # every rate out of A's charge state at B
        dwell = 6.0 / gamma / samples
    with ctx.run(2 * samples):
        run = dev.session(noise=False, charge="mean")
        run.hold(a, 200.0 / gamma_a)
        t0 = run.time_s
        m = run.sense(b, dwell, samples=samples)
        traj = dev.session(noise=False)
        traj.hold(a, 200.0 / gamma_a)
        s = traj.sense(b, dwell, samples=samples)
        ctx.count(model_evals=2 * samples)
    t = m.times_s - t0
    n_a, n_b = dev.at(a).n_mean[0], dev.at(b).n_mean[0]
    return {
        "t_s": t,
        "n1_mean": m.n_present[:, 0],
        "n1_analytic": n_b + (n_a - n_b) * np.exp(-gamma * t),
        "n1_trajectory": s.n_present[:, 0],
        "gamma_out_hz": np.array(gamma),
        "gamma_left_b_hz": np.asarray(dev.at(b).gamma_left_hz),
    }


def session_closing_barrier(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """The dots' rates follow the gates (#53). Along B_L below the tuned point, `Device.at` gives Γ_L,
    Γ_R and Γ_m. At dot 1's transition with B_L open and closed, P1 is swept up and then down at
    `dwell_s` per point in `mean` mode (noise off): open, the two directions agree; closed, the charge
    lags the gate both ways and the transition opens into a loop. Parked at the closed transition, a
    `sampled` session telegraphs at the same dwell's scale, noise-free and with the device's noise."""
    rates_sweep, p1_sweep, trace = p.sweep
    dwell = float(p.instance["dwell_s"])
    closed_dv = float(p.instance["closed_b_l_v"])
    trace_dwell = float(p.instance["trace_dwell_s"])
    with ctx.compile():
        dev = _compiled(p)
        v0 = _tuned_point(dev)
        jb, jp = dev.gate_names.index("B_L"), dev.gate_names.index("P1")
        db = np.linspace(float(rates_sweep["from"]), float(rates_sweep["to"]), int(rates_sweep["points"]))
        centres = {}
        for label, dv in (("open", 0.0), ("closed", closed_dv)):
            v = v0.copy()
            v[jb] += dv
            centres[label] = _dot1_degeneracy(dev, v0=v)
    n_p1, n_trace = int(p1_sweep["points"]), int(trace["points"])
    out: dict[str, np.ndarray] = {"v_b_l_v": v0[jb] + db, "dp1_v": np.linspace(-float(p1_sweep["half_span_v"]), float(p1_sweep["half_span_v"]), n_p1), "sample": np.arange(n_trace, dtype=float)}
    with ctx.run(len(db) + 4 * n_p1 + 2 * n_trace):
        line = np.broadcast_to(v0, db.shape + v0.shape).copy()
        line[:, jb] += db
        s = dev.at(line)
        out["gamma_left_hz"], out["gamma_right_hz"], out["gamma_interdot_hz"] = s.gamma_left_hz, s.gamma_right_hz, s.gamma_interdot_hz
        for label, c in centres.items():
            up = np.broadcast_to(c, (n_p1,) + c.shape).copy()
            up[:, jp] += out["dp1_v"]
            run = dev.session(noise=False, charge="mean")
            run.hold(up[0], 1.0)
            m_up = run.sense(up, dwell)
            m_down = run.sense(up[::-1], dwell)
            out[f"n1_up_{label}"], out[f"n1_down_{label}"] = m_up.n_present[:, 0], m_down.n_present[::-1, 0]
            out[f"current_up_{label}_a"], out[f"current_down_{label}_a"] = m_up.current_a, m_down.current_a[::-1]
        c = centres["closed"]
        out["current_trace_a"] = dev.session(noise=False).sense(c, trace_dwell, samples=n_trace).current_a
        out["current_trace_noisy_a"] = dev.session().sense(c, trace_dwell, samples=n_trace).current_a
        ctx.count(model_evals=len(db) + 4 * n_p1 + 2 * n_trace)
    out["gamma_left_closed_hz"] = np.asarray(dev.at(centres["closed"]).gamma_left_hz)
    out["gamma_left_open_hz"] = np.asarray(dev.at(centres["open"]).gamma_left_hz)
    return out
