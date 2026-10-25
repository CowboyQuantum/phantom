"""Rate-limited loading of a dot from its reservoir through a lead barrier, in a raster scan.

A minimal Step 1 model so the `plunger-vs-barrier-map` protocol can show what closing a lead barrier
does to a charge-stability scan: the tunnel rate between the dot and its bath falls exponentially as
the barrier closes, and once the rate drops below the inverse dwell time per point the charge state can
no longer follow the gate voltages. The measured occupation then stays at its previous value until the
dot gets a chance to equilibrate, including across the flyback at the end of a raster row, which is
what a latched or vanishing transition line in a real scan is.

Dot 2 is assumed to have no bath of its own: its electrons come and go only through dot 1, so when dot 1
is cut off from its reservoir the whole double dot's charge is frozen. Interdot rearrangement at a fixed
total charge is not modelled here.

Step 4 replaces this with lead tunnel rates and a master equation on the device clock; this module is
`calibration: none` and carries no claim beyond the exponential barrier dependence.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MODEL_VERSION = 1
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


def latched_occupation(
    n_equilibrium: np.ndarray,
    rate_hz: np.ndarray,
    dwell_s: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Occupation actually present at each point of a raster scan, in scan order.

    `n_equilibrium` (points, dots) is what the dots would hold if they could equilibrate at each point;
    `rate_hz` (points,) is the dot-to-bath rate there. At each point the dot equilibrates with
    probability 1 − exp(−Γ · dwell_s) and otherwise keeps the previous point's occupation. The first
    point always equilibrates (the scan starts from a settled device). Deterministic for a given `rng`.
    """
    n_eq = np.asarray(n_equilibrium, dtype=float)
    p_settle = 1.0 - np.exp(-np.asarray(rate_hz, dtype=float) * dwell_s)
    u = rng.random(len(p_settle))
    out = np.empty_like(n_eq)
    current = n_eq[0]
    for i in range(len(p_settle)):
        if i == 0 or u[i] < p_settle[i]:
            current = n_eq[i]
        out[i] = current
    return out
