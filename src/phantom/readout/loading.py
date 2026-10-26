"""Rate-limited loading of the dots from their reservoirs, and the lead-rate landscape over gate space.

A minimal Step 1 model so the stability protocols can show what the dot-to-lead tunnel rate does to a
scan: the rate falls exponentially as a barrier closes, and once it drops below the inverse of the time
the measurement spends at a point (or at one half-cycle of a chopped gate) the charge state can no
longer follow the gate voltages. The measured occupation then stays at its previous value until the dot
gets a chance to equilibrate, including across the flyback of a raster row, which is what a latched or
vanishing transition line in a real scan is. At the other extreme a very fast rate broadens the
transition by the lifetime width ħΓ.

Two forms:
- `LeadBarrier`: one barrier gate, Γ(V_B), used by `plunger-vs-barrier-map` (dot 2 has no bath of its
  own there and freezes with dot 1).
- `LeadRateLandscape`: a rate per dot over the whole gate vector, Γ_i(V) = rate0_i · exp(w_i · V),
  used by the chopped-gate stability maps, where the plungers themselves open and close the lead
  barriers through cross-coupling, so low-voltage corners latch and high-voltage corners broaden.

Lifetime broadening is applied as an energy width sqrt((kT)² + (ħΓ)²) in the Boltzmann occupation,
a thermal stand-in for the Lorentzian line shape. Step 4 replaces all of this with lead tunnel rates
and a master equation on the device clock; this module is `calibration: none` and carries no claim
beyond the exponential barrier dependence.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from phantom.constants import HBAR_EV_S, K_B_EV_PER_K

MODEL_VERSION = 2
"""Bumped whenever any output of this module changes for the same parameters and seed (ADR 0003)."""



@dataclass(frozen=True)
class LeadBarrier:
    """A barrier between a dot and its reservoir: tunnel rate Γ(V) = rate0_hz · exp((V − v_ref_v) / v_scale_v).

    `rate0_hz` is the rate at `v_ref_v`; `v_scale_v` is the gate voltage that changes it by a factor e.
    Range: rates 1 Hz–1 GHz over a few hundred mV of barrier gate; scales 10–100 mV.
    """

    name: str
    rate0_hz: float = 1e3
    v_ref_v: float = 0.0
    v_scale_v: float = 0.03

    def __post_init__(self) -> None:
        if self.rate0_hz <= 0 or self.v_scale_v <= 0:
            raise ValueError("rate0_hz and v_scale_v must be positive")

    def rate_hz(self, v_b: float | np.ndarray) -> np.ndarray:
        return self.rate0_hz * np.exp((np.asarray(v_b, dtype=float) - self.v_ref_v) / self.v_scale_v)


@dataclass(frozen=True)
class LeadRateLandscape:
    """Dot-to-lead tunnel rate of each dot as a function of the full gate vector:
    Γ_i(V) = rate0_hz[i] · exp(weights_per_v[i] · V). A weight of 1/0.01 per volt means the rate grows
    by a factor e every 10 mV of that gate. Range: rate0 1 Hz–1 MHz; weights 20–200 per volt."""

    gate_names: tuple[str, ...]
    rate0_hz: np.ndarray = field(default_factory=lambda: np.array([1e3, 1e3]))
    weights_per_v: np.ndarray = field(default_factory=lambda: np.array([[125.0, 0.0, 0.0], [0.0, 125.0, 0.0]]))

    def __post_init__(self) -> None:
        if self.weights_per_v.shape != (2, len(self.gate_names)) or self.rate0_hz.shape != (2,):
            raise ValueError("weights_per_v must be (2, gates) and rate0_hz (2,)")
        if np.any(self.rate0_hz <= 0):
            raise ValueError("rate0_hz must be positive")

    def rate_hz(self, v: np.ndarray) -> np.ndarray:
        """Γ (..., 2) in Hz for gate voltages `v` (..., gates)."""
        return self.rate0_hz * np.exp(np.asarray(v, dtype=float) @ self.weights_per_v.T)

    def broadened_kt_ev(self, v: np.ndarray, temperature_k: float) -> np.ndarray:
        """Energy width (...,) combining the thermal width with the lifetime width of the faster dot."""
        kt = K_B_EV_PER_K * max(temperature_k, 1e-6)
        gamma = self.rate_hz(v).max(axis=-1)
        return np.sqrt(kt**2 + (HBAR_EV_S * gamma) ** 2)


def latched_occupation(
    n_equilibrium: np.ndarray,
    rate_hz: np.ndarray,
    dwell_s: float,
    rng: np.random.Generator,
    per_dot: bool = False,
) -> np.ndarray:
    """Occupation actually present at each step of a scan, in scan order.

    `n_equilibrium` (steps, dots) is what the dots would hold if they could equilibrate at each step;
    `rate_hz` is the dot-to-bath rate there, (steps,) for one shared bath (`per_dot=False`: all dots
    settle together or not at all) or (steps, dots) with `per_dot=True` (each dot settles on its own).
    A dot equilibrates with probability 1 − exp(−Γ · dwell_s) and otherwise keeps the previous step's
    occupation. The first step always equilibrates (the scan starts from a settled device).
    Deterministic for a given `rng`.
    """
    n_eq = np.asarray(n_equilibrium, dtype=float)
    steps, dots = n_eq.shape
    gamma = np.asarray(rate_hz, dtype=float)
    if not per_dot:
        gamma = np.repeat(gamma[:, None], dots, axis=1)
    p_settle = 1.0 - np.exp(-gamma * dwell_s)
    u = rng.random(p_settle.shape) if per_dot else np.repeat(rng.random(steps)[:, None], dots, axis=1)
    settle = u < p_settle
    settle[0, :] = True
    out = np.empty_like(n_eq)
    current = n_eq[0].copy()
    for i in range(steps):
        current = np.where(settle[i], n_eq[i], current)
        out[i] = current
    return out
