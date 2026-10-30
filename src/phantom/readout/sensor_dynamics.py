"""The sensor dot's own charge as a process in time (#54, ADR 0028).

A sensor with leads (ADR 0026) holds M electrons, counted on the level numbering of φ_s: adding the
(M+1)-th costs μ_M = E_c (M − x), x = (φ_s − phase)/period, so M − 1 < x < M in its ground state and a
Coulomb peak is the degeneracy x = M. From M the dot gains an electron from the source at
Γ_S n_S(μ_M) and from the drain at Γ_D n_D(μ_M), and loses one to the source at Γ_S (1 − n_S(μ_{M−1}))
and to the drain at Γ_D (1 − n_D(μ_{M−1})), with n_S(μ) = n_γ(μ − eV_sd/2), n_D(μ) = n_γ(μ + eV_sd/2)
the leads' occupations convolved with the level's Lorentzian of width γ = ħ(Γ_S + Γ_D)
(`readout.broadening`). Through one level the stationary current of this chain is then exactly the
sensor's current of ADR 0026, e Γ_SΓ_D/(Γ_S + Γ_D) [n_S − n_D]; the levels the chain does not hold add
their stationary current as a constant (they carry the Lorentzian tails, elastic cotunnelling). This
is sequential tunnelling with the lifetime-broadened occupations put into the rates, a phenomenological
rate equation, exact in its stationary mean and in its thermal limit (ħΓ ≪ kT).

The charge a sample reads is the net number Q of electrons that leave to the drain during its dwell T,
I = eQ/T plus the other levels' current. When the chain makes at most `EXACT_JUMPS` jumps in a dwell
(`T · max out-rate`), Q and the end state come from an exact trajectory (one exponential waiting time
and one channel per jump); otherwise Q is Gaussian with the exact mean given the start state (the
integral of e^{Ws} by an augmented matrix exponential) and the long-time variance of full counting
statistics, σ² = 1ᵀJ₂p − 2·1ᵀJ₁ W^D J₁ p (Flindt et al. 2005: J₁, J₂ the counted jump rates weighted by
±1 and by 1, W^D the Drazin inverse), which for one level is the Fano factor (Γ_S² + Γ_D²)/(Γ_S + Γ_D)²
at large bias; the end state is drawn from e^{WT}. Calibration: none.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm

from phantom.constants import E_CHARGE_C, HBAR_EV_S
from phantom.readout.broadening import occupation
from phantom.readout.sensor import SensorSETParams

MODEL_VERSION = 1

EXACT_JUMPS = 200.0
"""Above this many expected jumps in a dwell the counted charge is Gaussian; at or below it, exact."""
RATE_FLOOR_HZ = 1e-300


@dataclass(frozen=True)
class SensorChain:
    """The sensor's charge states M = `m_lo` … `m_lo` + k − 1 at one gate point and dot charge."""

    m_lo: int
    w: np.ndarray
    """(k, k) generator: w[j, i] the rate from state i to j, Hz."""
    j1: np.ndarray
    """(k, k) counted rates: +rate for a jump that sends an electron into the drain, −rate for one
    that takes one from it, Hz."""
    j2: np.ndarray
    """(k, k) the counted rates' magnitudes, Hz."""
    background_a: float
    """Stationary current of the sensor levels the chain does not hold, A."""

    @property
    def states(self) -> np.ndarray:
        return self.m_lo + np.arange(len(self.w))


def _level_terms(sensor: SensorSETParams, phi: float, v: np.ndarray, v_sd: float):
    """(x, E_c, rates (Γ_S, Γ_D), γ in eV, kT, bias/2) at one point."""
    gs, gd = (float(r) for r in sensor.leads.rates_hz(v))
    return (phi - sensor.phase_v) / sensor.period_v, sensor.charging_energy_ev, gs, gd, HBAR_EV_S * (gs + gd), sensor.kt_ev, v_sd / 2


def _occupations(mu: np.ndarray, gamma_ev: float, kt: float, half: float) -> tuple[np.ndarray, np.ndarray]:
    return occupation(mu - half, gamma_ev, kt), occupation(mu + half, gamma_ev, kt)


def jump_rates_hz(sensor: SensorSETParams, phi: float, v: np.ndarray, v_sd: float, m: np.ndarray) -> np.ndarray:
    """(..., 4) rates out of states `m`: (in from the source, in from the drain, out to the source,
    out to the drain), Hz."""
    x, e_c, gs, gd, gamma_ev, kt, half = _level_terms(sensor, phi, v, v_sd)
    m = np.asarray(m, dtype=float)
    up_s, up_d = _occupations(e_c * (m - x), gamma_ev, kt, half)
    dn_s, dn_d = _occupations(e_c * (m - 1 - x), gamma_ev, kt, half)
    return np.stack([gs * up_s, gd * up_d, gs * (1 - dn_s), gd * (1 - dn_d)], axis=-1)


def window_states(sensor: SensorSETParams, phi: float, v_sd: float, margin_ev: float | None = None) -> tuple[int, int]:
    """(lowest, highest) charge the chain needs at potential `phi` and bias `v_sd`: both sides of every
    level within the bias window widened by `margin_ev` (20 kT by default), and at least both sides of
    the level nearest its centre."""
    x = (phi - sensor.phase_v) / sensor.period_v
    reach = (abs(v_sd) / 2 + (20 * sensor.kt_ev if margin_ev is None else margin_ev)) / sensor.charging_energy_ev
    k = int(np.rint(x))
    return min(k, int(np.ceil(x - reach))), max(k, int(np.floor(x + reach))) + 1


def chain(sensor: SensorSETParams, phi: float, v: np.ndarray, v_sd: float, m_lo: int, m_hi: int) -> SensorChain:
    """The generator over M = `m_lo` … `m_hi` at potential `phi`, carried voltages `v` (gates,) and bias
    `v_sd` (V), with the counted rates and the other levels' current."""
    x, e_c, gs, gd, gamma_ev, kt, half = _level_terms(sensor, phi, v, v_sd)
    k = m_hi - m_lo + 1
    w, j1 = np.zeros((k, k)), np.zeros((k, k))
    held = np.arange(m_lo, m_hi)  # the levels whose filling the chain follows: M → M + 1 at μ_M
    if len(held):
        n_s, n_d = _occupations(e_c * (held - x), gamma_ev, kt, half)
        i = np.arange(len(held))
        w[i + 1, i] = gs * n_s + gd * n_d
        w[i, i + 1] = gs * (1 - n_s) + gd * (1 - n_d)
        j1[i + 1, i] = -gd * n_d
        j1[i, i + 1] = gd * (1 - n_d)
    w[np.diag_indices(k)] = -w.sum(axis=0)
    # the other levels among the sensor's five, as `SensorSETParams.current_a` sums them
    mu = sensor.levels_ev(np.asarray(phi))
    levels = np.rint(x) + np.arange(-2, 3)
    others = mu[(levels < m_lo) | (levels >= m_hi)]
    n_s, n_d = _occupations(others, gamma_ev, kt, half)
    window = n_s - n_d
    window = np.maximum(window, 0.0) if half >= 0 else np.minimum(window, 0.0)
    background = E_CHARGE_C * gs * gd / max(gs + gd, RATE_FLOOR_HZ) * float(window.sum())
    return SensorChain(int(m_lo), w, j1, np.abs(j1), background)


def stationary(c: SensorChain) -> np.ndarray:
    """The chain's stationary distribution (k,): the detailed-balance product along the chain."""
    k = len(c.w)
    log_p = np.zeros(k)
    for i in range(k - 1):
        up, down = c.w[i + 1, i], c.w[i, i + 1]
        log_p[i + 1] = log_p[i] + np.log(max(up, RATE_FLOOR_HZ)) - np.log(max(down, RATE_FLOOR_HZ))
    p = np.exp(log_p - log_p.max())
    return p / p.sum()


def counting_cumulants(c: SensorChain, p: np.ndarray | None = None) -> tuple[float, float]:
    """(mean, variance) of the counted charge per second in the stationary state, electrons/s and
    electrons²/s: 1ᵀJ₁p and 1ᵀJ₂p − 2·1ᵀJ₁x with W x = J₁p − (1ᵀJ₁p) p, 1ᵀx = 0."""
    p = stationary(c) if p is None else p
    mean = float(c.j1.sum(axis=0) @ p)
    rhs = c.j1 @ p - mean * p
    k = len(p)
    a = np.vstack([c.w, np.ones((1, k))])
    x = np.linalg.lstsq(a, np.concatenate([rhs, [0.0]]), rcond=None)[0]
    return mean, float(c.j2.sum(axis=0) @ p - 2.0 * c.j1.sum(axis=0) @ x)


def propagate(c: SensorChain, q: np.ndarray, dwell_s: float) -> tuple[np.ndarray, float]:
    """(distribution after the dwell, expected counted charge in electrons) from distribution `q` (k,):
    e^{WT}q and 1ᵀJ₁∫₀ᵀe^{Ws}q ds from one matrix exponential of the generator augmented by q."""
    k = len(q)
    a = np.zeros((k + 1, k + 1))
    a[:k, :k] = c.w * dwell_s
    a[:k, k] = q * dwell_s
    e = expm(a)
    after = np.maximum(e[:k, :k] @ q, 0.0)
    return after / after.sum(), float(c.j1.sum(axis=0) @ e[:k, k])


def expected_jumps(c: SensorChain, dwell_s: float) -> float:
    return float(dwell_s * np.max(-np.diag(c.w)))


def trajectory(sensor: SensorSETParams, phi: float, v: np.ndarray, v_sd: float, m0: int, dwell_s: float, rng: np.random.Generator) -> tuple[int, int]:
    """(counted electrons, charge at the end) of one exact trajectory from charge `m0` over the dwell:
    each jump takes an exponential waiting time at the total out-rate and a channel by its weight."""
    m, q, t = int(m0), 0, 0.0
    while True:
        r = jump_rates_hz(sensor, phi, v, v_sd, np.array([m]))[0]
        total = float(r.sum())
        if not total > 0:
            return q, m
        t += rng.exponential(1.0 / total)
        if t > dwell_s:
            return q, m
        channel = int(np.searchsorted(np.cumsum(r), rng.random() * total, side="right"))
        channel = min(channel, 3)
        m += 1 if channel < 2 else -1
        if channel == 1:  # in from the drain
            q -= 1
        elif channel == 3:  # out to the drain
            q += 1


__all__ = [
    "EXACT_JUMPS",
    "MODEL_VERSION",
    "SensorChain",
    "chain",
    "counting_cumulants",
    "expected_jumps",
    "jump_rates_hz",
    "propagate",
    "stationary",
    "trajectory",
    "window_states",
]
