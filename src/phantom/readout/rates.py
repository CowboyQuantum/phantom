"""Tunnel rates through barriers of the landscape, in the WKB approximation.

A barrier of height `height_ev` above the dot level and width `width_nm` at that level transmits with
probability exp(−2 κ w), κ = √(2 m* ΔU) / ħ, and the electron attempts it at the frequency set by its
confinement, so

    Γ = f_attempt · exp(−2 κ w),   f_attempt ≈ ħω / h = ω / 2π   (ω from the dot's curvature).

`wkb_tunnel_rate_hz` is the rectangular-barrier estimate of ADR 0017 (height × width); a smooth barrier
of the same height and width transmits more, by one to two decades at the layouts used here.
`wkb_rate_from_action_hz` takes the WKB action S = ∫ κ ds along the actual barrier
(`electrostatics.barrier_action`) and gives Γ = f_attempt · exp(−2S); the compiled device and the
gallery's dynamics protocols use it (ADR 0023), and the rectangular form remains as its bound and for
the P0 limits. `GateRates` carries compiled rates to other gate voltages at runtime: the exponent
moves with the action's compiled slopes and `saturating_rates_hz` caps the rate at the attempt rate
(ADR 0026 for the sensor's leads, #53 for the qubit dots'). The interdot rate instead follows the
interdot barrier's one transmission, Γ_m = f_a·T (`charge.barrier`, #70): Kemble's transmission of the
barrier's screened top, the same T that sets t and the conductance. The lead barriers' screening waits
for the FEM backend (#48).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import expit

from phantom.charge.barrier import InterdotBarrier
from phantom.charge.shape import ShapeGrid
from phantom.constants import HBAR2_OVER_2ME_EV_NM2, HBAR_EV_S

MODEL_VERSION = 3
"""2: with an `InterdotBarrier`, the interdot rate is f_a·T of the barrier's one transmission (#70).
3: the interdot barrier is seen from the puddles' edges when their scale is given (#72)."""

ATTEMPT_RATE_FLOOR_HZ = 1e9
ATTEMPT_RATE_CEILING_HZ = 1e13


def wkb_transmission(height_ev: np.ndarray, width_nm: np.ndarray, effective_mass_me: float) -> np.ndarray:
    """exp(−2 κ w) for a rectangular barrier; 1 where the barrier is absent (height ≤ 0)."""
    h = np.maximum(np.asarray(height_ev, dtype=float), 0.0)
    w = np.maximum(np.asarray(width_nm, dtype=float), 0.0)
    kappa = np.sqrt(h * effective_mass_me / HBAR2_OVER_2ME_EV_NM2)  # 1/nm: κ² = 2m*ΔU/ħ²
    return np.exp(-2.0 * kappa * w)


def attempt_rate_hz(curvature_ev_per_nm2: float, effective_mass_me: float) -> float:
    """ω/2π for the harmonic dot of curvature k: ω = √(k / m*), in Hz, clipped to a plausible band."""
    hbar_omega_ev = np.sqrt(2 * HBAR2_OVER_2ME_EV_NM2 * max(curvature_ev_per_nm2, 1e-12) / effective_mass_me)
    f = hbar_omega_ev / (2 * np.pi * HBAR_EV_S)
    return float(np.clip(f, ATTEMPT_RATE_FLOOR_HZ, ATTEMPT_RATE_CEILING_HZ))


def wkb_rate_from_action_hz(action: np.ndarray, curvature_ev_per_nm2: float, effective_mass_me: float) -> np.ndarray:
    """Γ = f_attempt · exp(−2S) for the WKB action S (dimensionless) of the barrier seen from a dot of
    the given curvature."""
    return attempt_rate_hz(curvature_ev_per_nm2, effective_mass_me) * np.exp(-2.0 * np.asarray(action, dtype=float))


def wkb_tunnel_rate_hz(
    height_ev: np.ndarray, width_nm: np.ndarray, curvature_ev_per_nm2: float, effective_mass_me: float
) -> np.ndarray:
    """Γ = f_attempt · T for the barrier (height, width) seen from a dot of the given curvature."""
    return attempt_rate_hz(curvature_ev_per_nm2, effective_mass_me) * wkb_transmission(height_ev, width_nm, effective_mass_me)


def saturating_rates_hz(gamma_ref_hz: np.ndarray, attempt_rate_hz: np.ndarray | float, exponent: np.ndarray) -> np.ndarray:
    """f_a · σ(ln(Γ_ref / (f_a − Γ_ref)) + x), σ(x) = 1/(1 + e⁻ˣ): a tunnel rate whose exponent moves by
    `exponent` (..., k) from the reference rate `gamma_ref_hz` (k,), exponential while the barrier is
    closed and saturating at the attempt rate f_a (k,) as it opens; exactly Γ_ref at x = 0 (ADR 0026)."""
    g_ref = np.asarray(gamma_ref_hz, dtype=float)
    f_a = np.asarray(attempt_rate_hz, dtype=float)
    anchor = np.log(g_ref) - np.log(np.maximum(f_a - g_ref, 1e-12 * f_a))
    return f_a * expit(anchor + exponent)


@dataclass(frozen=True)
class GateRates:
    """Tunnel rates that follow the gate voltages (#53): k barriers compiled at a reference point, each
    with the WKB rate there, the attempt rate of the dot it leaves and the slopes d ln Γ/dV of its
    action against every gate at a tunnelling level held fixed (the leads' electrochemical potential,
    as for the sensor's leads, ADR 0026). At other voltages the exponent moves linearly with the slopes
    and `saturating_rates_hz` caps each rate at its attempt rate. Only the prefactors of the rate
    equation move, so detailed balance and the Boltzmann stationary state are untouched."""

    v_ref_v: tuple[float, ...]
    """Carried gate voltages the rates are compiled at, one per gate, V."""
    gamma_ref_hz: tuple[float, ...]
    """Rate through each barrier at `v_ref_v`, Hz."""
    attempt_rate_hz: tuple[float, ...]
    """Attempt frequency of the dot each barrier is seen from, Hz: the ceiling of its rate."""
    log_rate_slopes_per_v: tuple[tuple[float, ...], ...]
    """d ln Γ/dV of each barrier with respect to every gate, 1/V."""
    shape: ShapeGrid | None = None
    """When given (three barriers: left lead, right lead, interdot), the exponents along its barrier
    gates are the changes of the actions on its grid (#59, the dots moved and reshaped), and the slopes
    count for the other gates only."""
    interdot: InterdotBarrier | None = None
    """When given (three barriers), the third barrier's exponent is the barrier's `log_odds_change`, so
    its rate is f_a·T of the one interdot transmission (#70); `gamma_ref_hz[2]` is then f_a·T at the
    reference point."""

    def __post_init__(self) -> None:
        k, n = len(self.gamma_ref_hz), len(self.v_ref_v)
        if len(self.attempt_rate_hz) != k or len(self.log_rate_slopes_per_v) != k:
            raise ValueError("one attempt rate and one row of slopes per barrier")
        if any(len(row) != n for row in self.log_rate_slopes_per_v):
            raise ValueError("slopes need one entry per gate")
        if not all(f > 0 for f in self.attempt_rate_hz) or not all(0 < g for g in self.gamma_ref_hz):
            raise ValueError("rates must be positive")
        if self.interdot is not None and k != 3:
            raise ValueError("an interdot barrier needs the three barriers: left lead, right lead, interdot")

    def rates_hz(self, v: np.ndarray, interdot_scale: np.ndarray | None = None) -> np.ndarray:
        """Γ of each barrier (..., k), Hz, at carried voltages `v` (..., gates); `interdot_scale` (..., 2)
        is the puddles' scale the interdot barrier is seen from (#72)."""
        v = np.asarray(v, dtype=float)
        dv = v - np.asarray(self.v_ref_v)
        slopes = np.array(self.log_rate_slopes_per_v)
        if self.shape is None:
            exponent = dv @ slopes.T
        else:
            slopes[:, list(self.shape.gate_indices)] = 0.0
            exponent = dv @ slopes.T + self.shape.quantity(self.shape.change(v), "log_rate")
        if self.interdot is not None:
            exponent = np.concatenate([exponent[..., :2], self.interdot.log_odds_change(v, interdot_scale)[..., None]], axis=-1)
        return saturating_rates_hz(np.asarray(self.gamma_ref_hz), np.asarray(self.attempt_rate_hz), exponent)
