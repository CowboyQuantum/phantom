"""Step 8 protocols: the consumer session on the compiled device — the charge a sample sees, and the charge after a hold."""

from __future__ import annotations

import numpy as np

from gallery.protocols.step7 import _compiled
from gallery.runner import Context, Protocol
from phantom.device import Device


def _dot1_degeneracy(dev: Device, half_v: float = 0.03, n: int = 4001) -> np.ndarray:
    """The programmed point nearest the programmed point where dot 1 is half-way through a transition to its lead."""
    v0 = np.array(dev.programmed_point_v)
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
        gamma = dev.at(v).gamma_left_hz
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
    the mean charge of dot 1 (noise off) relaxes as n_B + (n_A − n_B)·exp(−Γ_L·t) (drawn beside it), and
    one sampled trajectory steps once at a random time."""
    samples = int(p.sweep[0]["points"])
    with ctx.compile():
        dev = _compiled(p)
        v0 = np.array(dev.programmed_point_v)
        j = dev.gate_names.index("P1")
        per_e = dev.truth().charge_lattice_e_per_v[0][0]
        a = _dot1_degeneracy(dev) + 0.5 / per_e * np.eye(len(v0))[j]
        b = a - 1.0 / per_e * np.eye(len(v0))[j]
        gamma = dev.at(a).gamma_left_hz
        dwell = 6.0 / gamma / samples
    with ctx.run(2 * samples):
        run = dev.session(noise=False, charge="mean")
        run.hold(a, 200.0 / gamma)
        t0 = run.time_s
        m = run.sense(b, dwell, samples=samples)
        traj = dev.session(noise=False)
        traj.hold(a, 200.0 / gamma)
        s = traj.sense(b, dwell, samples=samples)
        ctx.count(model_evals=2 * samples)
    t = m.times_s - t0
    n_a, n_b = dev.at(a).n_mean[0], dev.at(b).n_mean[0]
    return {
        "t_s": t,
        "n1_mean": m.n_present[:, 0],
        "n1_analytic": n_b + (n_a - n_b) * np.exp(-gamma * t),
        "n1_trajectory": s.n_present[:, 0],
        "gamma_left_hz": np.array(gamma),
    }
