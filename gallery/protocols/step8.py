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


def _dot1_degeneracy(dev: Device, half_v: float = 0.03, n: int = 4001, v0: np.ndarray | None = None, n1_from: int | None = None) -> np.ndarray:
    """The programmed point along P1 nearest `v0` (the tuned point when None) where dot 1 is half-way
    through a transition to its lead; with `n1_from`, the transition from that many electrons to one
    more, so a shifted point keeps following the same line."""
    v0 = _tuned_point(dev) if v0 is None else v0
    j = dev.gate_names.index("P1")
    grid = v0 + np.linspace(-half_v, half_v, n)[:, None] * np.eye(len(v0))[j]
    occ = dev.at(grid).n_mean
    frac = occ[:, 0] - np.floor(occ[:, 0])
    flat2 = np.abs(occ[:, 1] - np.rint(occ[:, 1])) < 1e-2  # dot 2 settled, up to its interdot admixture (#58)
    if n1_from is not None:
        flat2 &= np.floor(occ[:, 0]) == n1_from
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
    Γ_R and Γ_m. At dot 1's transition, P1 is swept up and then down at `dwell_s` per point in `mean`
    mode (noise off) three ways: barriers open, the two directions agree; B_L closed alone, they agree
    too, because dot 2's lead loads dot 1 through the interdot coupling at Γ_R(t/Δ)² (#58); B_L and B_M
    closed, dot 1 is isolated, the charge lags the gate both ways and the transition opens into a loop.
    Parked at the isolated transition, a `sampled` session telegraphs at the same dwell's scale,
    noise-free and with the device's noise."""
    rates_sweep, p1_sweep, trace = p.sweep
    dwell = float(p.instance["dwell_s"])
    alone = {"B_L": float(p.instance["closed_b_l_v"])}
    isolated = {g: float(x) for g, x in p.instance["isolated_v"].items()}
    trace_dwell = float(p.instance["trace_dwell_s"])
    with ctx.compile():
        dev = _compiled(p)
        v0 = _tuned_point(dev)
        jb, jp = dev.gate_names.index("B_L"), dev.gate_names.index("P1")
        db = np.linspace(float(rates_sweep["from"]), float(rates_sweep["to"]), int(rates_sweep["points"]))
        centres: dict[str, np.ndarray] = {}
        for label, shift in (("open", {}), ("closed", alone), ("isolated", isolated)):
            v = v0.copy()
            for gate, dv in shift.items():
                v[dev.gate_names.index(gate)] += dv
            # the same transition of dot 1 each time: the barriers shift it along P1 (cross-capacitance)
            n1 = None if not centres else int(np.floor(dev.at(centres["open"]).n_mean[0]))
            centres[label] = _dot1_degeneracy(dev, half_v=0.06, n=12001, v0=v, n1_from=n1)
    n_p1, n_trace = int(p1_sweep["points"]), int(trace["points"])
    out: dict[str, np.ndarray] = {"v_b_l_v": v0[jb] + db, "dp1_v": np.linspace(-float(p1_sweep["half_span_v"]), float(p1_sweep["half_span_v"]), n_p1), "sample": np.arange(n_trace, dtype=float)}
    with ctx.run(len(db) + 6 * n_p1 + 2 * n_trace):
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
        c = centres["isolated"]
        out["current_trace_a"] = dev.session(noise=False).sense(c, trace_dwell, samples=n_trace).current_a
        out["current_trace_noisy_a"] = dev.session().sense(c, trace_dwell, samples=n_trace).current_a
        ctx.count(model_evals=len(db) + 6 * n_p1 + 2 * n_trace)
    for label, c in centres.items():
        s = dev.at(c)
        out[f"gamma_left_{label}_hz"] = np.asarray(s.gamma_left_hz)
        out[f"tunnel_coupling_{label}_ev"] = np.asarray(s.tunnel_coupling_ev)
    return out


def _interdot_centre(dev: Device, v0: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(the middle of the interdot segment of the honeycomb cell around `v0`, the detuning direction as
    a unit vector in programmed volts): the CI energies depend on N − n_g alone, so the segment's middle
    is floor(n_g) + ½(tp₁ + tp₂) in induced charge, reached with P1 and P2; along the direction n_g moves
    by (½, −½), one electron from dot 2 to dot 1 at a fixed total."""
    ci = dev.runtime_ci()
    j = [dev.gate_names.index(g) for g in ("P1", "P2")]

    def ng(v: np.ndarray) -> np.ndarray:
        return ci.induced_charge(dev.effective_voltages(v))

    jac = np.stack([ng(v0 + np.eye(len(v0))[k]) - ng(v0) for k in j], axis=1)  # electrons per volt, linear
    tp1, tp2 = ci.triple_points_ng()
    v = v0.copy()
    v[j] += np.linalg.solve(jac, np.floor(ng(v0)) + 0.5 * (tp1 + tp2) - ng(v0))
    d = np.zeros(len(v0))
    d[j] = np.linalg.solve(jac, np.array([0.5, -0.5]))
    return v, d / np.linalg.norm(d)


def interdot_line_vs_barrier(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """The interdot coupling follows the gates and hybridizes the charge states (#58, ADR 0029). Along
    B_M around the tuned point, `Device.at` gives t and Γ_m. Detuning cuts through the middle of the
    interdot segment of the honeycomb cell at three B_M give dot 1's charge beside the pair's two-level
    polarization ½(1 − ε/Ω·tanh(Ω/2kT)), Ω = √(ε² + 4t²), and the sensor current. A map gives dot 1's
    charge over detuning and B_M, each row centred on its cell's segment."""
    bm_sweep, cut, rows = p.sweep
    with ctx.compile():
        dev = _compiled(p)
        v0 = _tuned_point(dev)
        jb = dev.gate_names.index("B_M")
        ci = dev.runtime_ci()
        kt = 8.617333262e-5 * dev.temperature_k
        dbm = np.linspace(float(bm_sweep["from"]), float(bm_sweep["to"]), int(bm_sweep["points"]))
        x = np.linspace(-float(cut["half_span_v"]), float(cut["half_span_v"]), int(cut["points"]))
        cuts = [float(b) for b in cut["at"]]
        rows_v = np.linspace(float(rows["from"]), float(rows["to"]), int(rows["points"]))
        xm = np.linspace(-float(rows["half_span_v"]), float(rows["half_span_v"]), int(rows["columns"]))

        def centred(db: float) -> tuple[np.ndarray, np.ndarray]:
            v = v0.copy()
            v[jb] += db
            return _interdot_centre(dev, v)

        cut_centres = [centred(b) for b in cuts]
        row_centres = [centred(b) for b in rows_v]
    n_points = len(dbm) + len(cuts) * len(x) + len(rows_v) * len(xm)
    out: dict[str, np.ndarray] = {"v_b_m_v": v0[jb] + dbm, "x_v": x, "row_b_m_v": v0[jb] + rows_v, "xm_v": xm}
    with ctx.run(n_points):
        line = np.broadcast_to(v0, dbm.shape + v0.shape).copy()
        line[:, jb] += dbm
        s = dev.at(line)
        out["tunnel_coupling_ev"], out["gamma_interdot_hz"] = s.tunnel_coupling_ev, s.gamma_interdot_hz
        for k, (c, d) in enumerate(cut_centres):
            s = dev.at(c + x[:, None] * d)
            base = np.floor(s.n_mean[len(x) // 2]).astype(int)  # the pair shares one electron above it
            a, b = np.array([base[0] + 1, base[1]]), np.array([base[0], base[1] + 1])
            eps = np.array([float(np.diff(ci.state_energies_ev(v, np.stack([b, a])))[0]) for v in s.v_carried_v])
            omega = np.sqrt(eps**2 + 4 * s.tunnel_coupling_ev**2)
            out[f"n1_cut{k}"] = s.n_mean[:, 0] - base[0]
            out[f"two_level_cut{k}"] = 0.5 * (1 - eps / omega * np.tanh(omega / (2 * kt)))
            out[f"current_cut{k}_a"] = s.sensor_current_a
            out[f"tunnel_coupling_cut{k}_ev"] = np.asarray(s.tunnel_coupling_ev[len(x) // 2])
        grid = np.stack([c + xm[:, None] * d for c, d in row_centres])
        n = dev.at(grid).n_mean
        out["n1_map"] = (n[..., 0] - np.floor(n[:, len(xm) // 2, 0])[:, None]).T  # [detuning, B_M]
        ctx.count(model_evals=n_points)
    return out


def _sensor_level_at(dev: Device, v: np.ndarray, offset_kt: float) -> np.ndarray:
    """`v` with the sensor gate moved so the sensor's nearest level sits `offset_kt` kT above the
    leads, with the dots in their ground state (the gallery device's sensor gate is not faulted)."""
    sensor, j = dev.runtime_sensor(), dev.gate_names.index("SG")
    out = v.copy()
    for _ in range(2):
        s = dev.at(out)
        x = (float(sensor.potential(s.v_carried_v, s.n_ground)) - sensor.phase_v) / sensor.period_v
        out[j] += (np.rint(x) - offset_kt * sensor.kt_ev / sensor.charging_energy_ev - x) * sensor.period_v / sensor.alpha_v_per_v[j]
    return out


def session_slow_sensor(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """The sensor dot's charge in time (#54, ADR 0028). Its two barriers close together along the first
    sweep (`State.sensor_gamma_hz`). Parked with the sensor's level `offset_kt` kT above the leads, a
    `sampled` session (noise off) reads a smooth current with shot noise when the barriers are open and
    single-electron bursts with the sensor's own telegraph when they are closed. Then P1 is swept across
    dot 1's transition at `sweep_dwell_s` per point, read by the open and by the closed sensor."""
    closing, trace, p1_sweep = p.sweep
    inst = p.instance
    closed = {g: float(dv) for g, dv in inst["closed_v"].items()}
    offset, dwell, sweep_dwell = float(inst["offset_kt"]), float(inst["trace_dwell_s"]), float(inst["sweep_dwell_s"])
    with ctx.compile():
        dev = _compiled(p)
        v0 = _tuned_point(dev)
        frac = np.linspace(float(closing["from"]), float(closing["to"]), int(closing["points"]))
        names = dev.gate_names
        line = np.broadcast_to(v0, frac.shape + v0.shape).copy()
        for g, dv in closed.items():
            line[:, names.index(g)] += frac * dv
        shut = v0.copy()
        for g, dv in closed.items():
            shut[names.index(g)] += dv
        parked = {"open": _sensor_level_at(dev, v0, offset), "closed": _sensor_level_at(dev, shut, offset)}
        jp = names.index("P1")
        centre = {k: _dot1_degeneracy(dev, v0=x) for k, x in parked.items()}
    n_trace, n_p1 = int(trace["points"]), int(p1_sweep["points"])
    dp1 = np.linspace(-float(p1_sweep["half_span_v"]), float(p1_sweep["half_span_v"]), n_p1)
    out: dict[str, np.ndarray] = {"closing_fraction": frac, "sample": np.arange(n_trace, dtype=float), "dp1_v": dp1}
    with ctx.run(len(frac) + 2 * n_trace + 2 * n_p1):
        rates = dev.at(line).sensor_gamma_hz
        out["gamma_source_hz"], out["gamma_drain_hz"] = rates[:, 0], rates[:, 1]
        for k, x in parked.items():
            m = dev.session(noise=False).sense(x, dwell, samples=n_trace)
            out[f"current_{k}_a"], out[f"sensor_charge_{k}"] = m.current_a, m.sensor_charge - np.min(m.sensor_charge)
            up = np.broadcast_to(centre[k], (n_p1,) + v0.shape).copy()
            up[:, jp] += dp1
            s = dev.session(noise=False).sense(up, sweep_dwell)
            out[f"sweep_current_{k}_a"], out[f"sweep_n1_{k}"] = s.current_a, s.n_present[:, 0]
        ctx.count(model_evals=len(frac) + 2 * n_trace + 2 * n_p1)
    out["sensor_rates_closed_hz"] = dev.at(parked["closed"]).sensor_gamma_hz
    out["sensor_rates_open_hz"] = dev.at(parked["open"]).sensor_gamma_hz
    return out
