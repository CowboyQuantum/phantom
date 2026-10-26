"""Single-shot spin readout by spin-dependent tunnelling (Elzerman readout) (ADR 0022).

The dot level is pulsed so that the Fermi level of the lead lies between the two spin states: an
electron in the upper spin state (|↑⟩ here) can tunnel out at rate Γ_out, after which a lower-state
electron tunnels in at rate Γ_in; the charge sensor sees a blip, empty for the interval between the two
events. A lower-state electron stays put and the trace is flat. Spin relaxation at rate 1/T₁ before the
tunnel-out event turns a blip into no blip. The ensemble-average probability that the dot is empty at
time t after the pulse, for an upper-state electron, is

    P_empty(t) = Γ_out (e^{−Γ_in t} − e^{−λ t}) / (λ − Γ_in),    λ = Γ_out + 1/T₁,

which is the curve a consumer thresholds to set the read window. The rates are the lead rates of the
layout (`readout.rates`) at the read point; here they are parameters.
"""

from __future__ import annotations

import numpy as np

MODEL_VERSION = 1


def empty_probability(t_s: np.ndarray, gamma_out_hz: float, gamma_in_hz: float, t1_s: float) -> np.ndarray:
    """P(dot empty at t | upper spin state at t = 0)."""
    t = np.asarray(t_s, dtype=float)
    lam = gamma_out_hz + 1.0 / t1_s
    if abs(lam - gamma_in_hz) < 1e-9 * lam:
        return gamma_out_hz * t * np.exp(-lam * t)
    return gamma_out_hz * (np.exp(-gamma_in_hz * t) - np.exp(-lam * t)) / (lam - gamma_in_hz)


def elzerman_shots(
    rng: np.random.Generator,
    n_shots: int,
    p_upper: float,
    gamma_out_hz: float,
    gamma_in_hz: float,
    t1_s: float,
    dt_s: float,
    n_samples: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Charge traces (n_shots, n_samples) of 1 (occupied) and 0 (empty) sampled every `dt_s`, and the
    prepared spin of each shot (1 upper, 0 lower). Event times are exponential draws in a fixed order
    (spin, relaxation, tunnel-out, tunnel-in), so the shots are a deterministic function of the seed."""
    t = np.arange(n_samples) * dt_s
    traces = np.ones((n_shots, n_samples), dtype=np.int8)
    spins = (rng.uniform(size=n_shots) < p_upper).astype(np.int8)
    for k in range(n_shots):
        t_relax = rng.exponential(t1_s)
        t_out = rng.exponential(1.0 / gamma_out_hz)
        t_in = rng.exponential(1.0 / gamma_in_hz)
        if spins[k] == 1 and t_out < t_relax:
            traces[k, (t >= t_out) & (t < t_out + t_in)] = 0
    return traces, spins
