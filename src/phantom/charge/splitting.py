"""The interdot conductance sets the triple-point splitting, from a grid to a superdot (#66, ADR 0032).

In double-dot transport the interdot segment of the honeycomb, the separation of the triple points
along anti-detuning, grows as the barrier between the dots opens. Two decoupled dots give a grid;
hexagons form as the dots couple, and two merged dots give the parallel diagonal lines of one superdot
(van der Wiel et al., RMP 2002). In the measurements the splitting follows the interdot conductance
(Waugh et al. 1995). For two equivalent dots, Golden & Halperin relate the fractional peak splitting f to
the dimensionless conductance g of each of N_ch interdot channels:

    weak coupling:          f = (2 ln 2/π²)·N_ch·g,
    strong coupling, N_ch = 2:  f = 1 + (16e^γ/π³)·(1 − g)·ln(1 − g),

and they warn that the two limits do not interpolate reliably. `golden_halperin_splitting` joins them
with a monotone cubic (PCHIP) through points of each limit inside its range of validity
(`calibration: none`).

The conductance per channel is Kemble's transmission through a parabolic barrier,
g = 1/(1 + exp(2πΔU/ħω_b)), with ΔU the barrier top above the tunnelling level held fixed (ADR 0027).
ΔU turns negative once the barrier falls below the level, so g goes continuously from 0 (closed) through
½ (the level at the top) towards 1 (the dots merging). `InterdotSplitting` reads ΔU and g from an
`InterdotBarrier` (`charge.barrier`, #70): the compiled top moved by every gate through screened slopes
and, along the barrier gates, the shape grid continued past the merge (ADR 0030, 0034). It is the same
transmission that sets t and the interdot rate.

The splitting enters the CI model as an effective interdot capacitance. For CI capacitances
[[C_t, −C_m], [−C_m, C_t]] the classical fractional splitting is f_C = 2C_m/(C_t + C_m), taken here with
C_t the mean of the two totals. The target f = f_C + (F_MAX − f_C)·f_GH(g) is reached by adding
ΔC_m·[[1, −1], [−1, 1]]: the mutual capacitance and both totals grow together, which is the superdot
limit C_m → ∞, and

    ΔC_m = (F_MAX − f_C)·f_GH(g)·(C_t + C_m) / (2(1 − f)),

which is exactly zero when g = 0. F_MAX < 1, so the superdot is approached, not reached. The
gate-induced charge keeps the geometric capacitances: the tunnelling relaxes the charge quantization, it
does not move the gates' charge.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import PchipInterpolator

from phantom.charge.barrier import InterdotBarrier, kemble_transmission

MODEL_VERSION = 3
"""2: ΔU and g come from the `InterdotBarrier`, screened and continued past the merge (#70, ADR 0034).
3: the puddles' occupancy scale reaches the barrier (#72, ADR 0035)."""

EULER_GAMMA = 0.5772156649015329
F_MAX = 0.98
"""The largest splitting the CI model is taken to: the superdot is approached, not reached (C_m finite)."""
WEAK_UP_TO = 0.15
"""Golden–Halperin's weak-coupling form is used for g up to here."""
STRONG_FROM = 0.85
"""And their strong-coupling form from here."""


def weak_splitting(g: np.ndarray, n_channels: int = 2) -> np.ndarray:
    """f = (2 ln 2/π²)·N_ch·g: Golden & Halperin's weak-coupling limit."""
    return 2 * np.log(2) / np.pi**2 * n_channels * np.asarray(g, dtype=float)


def strong_splitting(g: np.ndarray) -> np.ndarray:
    """f = 1 + (16e^γ/π³)(1 − g)·ln(1 − g): Golden & Halperin's strong-coupling limit for N_ch = 2."""
    one_minus = 1 - np.asarray(g, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(one_minus > 0, one_minus * np.log(np.where(one_minus > 0, one_minus, 1.0)), 0.0)
    return 1 + 16 * np.exp(EULER_GAMMA) / np.pi**3 * term


def _knots() -> tuple[np.ndarray, np.ndarray]:
    weak = np.linspace(0.0, WEAK_UP_TO, 7)
    strong = np.concatenate([np.linspace(STRONG_FROM, 0.99, 8), [0.995, 0.999, 1.0]])
    return np.concatenate([weak, strong]), np.concatenate([weak_splitting(weak), strong_splitting(strong)])


_G, _F = _knots()
_CURVE = PchipInterpolator(_G, _F, extrapolate=False)


def golden_halperin_splitting(g: np.ndarray) -> np.ndarray:
    """The fractional splitting f(g) for N_ch = 2 (spin), g per channel in [0, 1]: Golden & Halperin's
    weak limit up to g = 0.15 and strong limit from g = 0.85, joined monotonically (PCHIP)."""
    return _CURVE(np.clip(np.asarray(g, dtype=float), 0.0, 1.0))


@dataclass(frozen=True)
class InterdotSplitting:
    """The interdot barrier's conductance at the gates and the splitting it sets."""

    barrier: InterdotBarrier
    """The barrier between the dots: its top at the gates and its transmission (`charge.barrier`)."""

    def barrier_top_ev(self, v: np.ndarray, scale: np.ndarray | None = None) -> np.ndarray:
        """ΔU (...) at carried voltages `v` (..., gates), eV; with the puddles' `scale` (..., 2) (#72)."""
        return self.barrier.barrier_top_ev(v, scale)

    def conductance(self, v: np.ndarray, scale: np.ndarray | None = None) -> np.ndarray:
        """g (...) per channel at carried voltages `v`: the barrier's transmission."""
        return self.barrier.transmission(v, scale)

    def splitting(self, v: np.ndarray, c: np.ndarray, scale: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
        """(f (...), ΔC_m (...) in aF): the total fractional splitting and the interdot capacitance that
        sets it, given the geometric capacitances `c` (..., 2, 2) at carried voltages `v` and the
        puddles' `scale` (..., 2) when the dots grow with their occupancy (#72)."""
        c = np.asarray(c, dtype=float)
        c_t = 0.5 * (c[..., 0, 0] + c[..., 1, 1])
        c_m = -c[..., 0, 1]
        f_c = 2 * c_m / (c_t + c_m)
        f_gh = golden_halperin_splitting(self.conductance(v, scale))
        f = f_c + (F_MAX - f_c) * f_gh
        return f, (F_MAX - f_c) * f_gh * (c_t + c_m) / (2 * (1 - f))


__all__ = [
    "F_MAX",
    "MODEL_VERSION",
    "InterdotBarrier",
    "InterdotSplitting",
    "golden_halperin_splitting",
    "kemble_transmission",
    "strong_splitting",
    "weak_splitting",
]
