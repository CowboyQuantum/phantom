"""The interdot barrier at the gates: its top above the tunnelling level, screened, continued past the merge, and the one transmission that t, Γ_m and the conductance share (#70, ADR 0034).

Compile records the top ΔU of the barrier between the dots above the interdot tunnelling level held
fixed (ADR 0027), its slopes against every gate with the dots held, and its curvature ħω_b (#66). The
shape grid (ADR 0030) records the top along each barrier gate with the dots moved. `InterdotBarrier`
gives ΔU at any carried gate point:

    ΔU(V) = ΔU_ref + Σ_g σ_g · δ_g(V_g),

where δ_g is the gate's change of the top, from the shape grid along a barrier gate and from the
compiled slope for every other gate, and σ_g is a *screening factor* per gate. The closed-form potential
has no screening by the electrons (ADR 0013), so in it the plungers move the interdot barrier nearly as
much as the barrier gate does, and every gate moves it fast (#48, #70). Until the solver exists, compile
sets σ_g from two placeholders (`calibration: none`): the barrier gate's own factor and, for every
other gate, that factor times a cross-coupling factor. σ = 1 for every gate is the bare potential.

Beyond an end of the shape grid the dots' shape is held (ADR 0030), but the barrier's top is the
potential on the cut between the held dots, which the gates keep moving: along a barrier gate ΔU is
continued linearly past the grid's end with the interpolant's end slope. So a barrier gate opened far
past the merge keeps lowering the top, and no plunger brings the two dots back (#70). The merged flag
still marks the held shape.

The transmission is Kemble's for a parabolic top, T = 1/(1 + exp(2πΔU/ħω_b)): e^(−2πΔU/ħω_b) deep
below the top (the WKB action of a parabola is πΔU/ħω_b), ½ when the level reaches the top, and
continuously towards 1 as the barrier sinks. It is the interdot barrier's one transmission: the
conductance per channel g = T (`charge.splitting`), the coupling t = (ħω₀/π)·√T (`charge.hybrid`) and
the interdot rate Γ_m = f_a·T (`readout.rates`). On the reference device the WKB action's compiled slopes
agree with 2π·dΔU/dV/ħω_b within 10–30 % on every gate.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phantom.charge.shape import ShapeGrid

MODEL_VERSION = 1


def kemble_transmission(top_ev: np.ndarray, hbar_omega_b_ev: float) -> np.ndarray:
    """1/(1 + exp(2πΔU/ħω_b)): transmission through a parabolic barrier whose top is ΔU above the
    electron's energy (negative below it), ħω_b its curvature frequency (Kemble 1935)."""
    x = np.clip(2 * np.pi * np.asarray(top_ev, dtype=float) / hbar_omega_b_ev, -700, 700)
    return 1.0 / (1.0 + np.exp(x))


@dataclass(frozen=True)
class InterdotBarrier:
    """The barrier between the dots at the carried gate voltages: its top above the fixed interdot
    tunnelling level and Kemble's transmission of it."""

    v_ref_v: tuple[float, ...]
    """Carried gate voltages the barrier is compiled at, V."""
    top_ev: float
    """ΔU at `v_ref_v`: the top above the interdot tunnelling level, eV (negative below it)."""
    top_slopes_ev_per_v: tuple[float, ...]
    """dΔU/dV against every gate in the bare potential, the dots held, eV/V."""
    hbar_omega_b_ev: float
    """ħω_b of the top's curvature, eV."""
    screening: tuple[float, ...] | None = None
    """σ_g per gate, the fraction of each gate's bare action on the top that reaches it; None is 1 for
    every gate (the bare potential)."""
    shape: ShapeGrid | None = None
    """When given, ΔU along its barrier gates follows its `barrier_top_ev` column, continued linearly
    past each end, and the slopes count for the other gates only."""

    def __post_init__(self) -> None:
        n = len(self.v_ref_v)
        if len(self.top_slopes_ev_per_v) != n or (self.screening is not None and len(self.screening) != n):
            raise ValueError("slopes and screening need one entry per gate")
        if not self.hbar_omega_b_ev > 0:
            raise ValueError("hbar_omega_b_ev must be positive")

    def _sigma(self) -> np.ndarray:
        return np.ones(len(self.v_ref_v)) if self.screening is None else np.asarray(self.screening, dtype=float)

    def barrier_top_ev(self, v: np.ndarray) -> np.ndarray:
        """ΔU (...) at carried voltages `v` (..., gates), eV."""
        v = np.asarray(v, dtype=float)
        sigma = self._sigma()
        slopes = np.asarray(self.top_slopes_ev_per_v, dtype=float) * sigma
        along = 0.0
        if self.shape is not None:
            slopes[list(self.shape.gate_indices)] = 0.0
            along = self.shape.quantity_along(v, "barrier_top_ev", weights=sigma, beyond="linear")[..., 0]
        return self.top_ev + (v - np.asarray(self.v_ref_v)) @ slopes + along

    def transmission(self, v: np.ndarray) -> np.ndarray:
        """T (...) at carried voltages `v`: Kemble's transmission of the top there."""
        return kemble_transmission(self.barrier_top_ev(v), self.hbar_omega_b_ev)

    def log_odds_change(self, v: np.ndarray) -> np.ndarray:
        """ln(T/(1 − T)) − ln(T_ref/(1 − T_ref)) = −2π(ΔU(V) − ΔU_ref)/ħω_b (...): the exponent by which a
        rate f_a·T moves from its reference in `readout.rates.saturating_rates_hz`."""
        return -2 * np.pi * (self.barrier_top_ev(v) - self.top_ev) / self.hbar_omega_b_ev


__all__ = ["MODEL_VERSION", "InterdotBarrier", "kemble_transmission"]
