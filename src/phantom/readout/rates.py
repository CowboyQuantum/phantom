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
the P0 limits. Screening waits for the FEM backend.
"""

from __future__ import annotations

import numpy as np

from phantom.constants import HBAR2_OVER_2ME_EV_NM2, HBAR_EV_S

MODEL_VERSION = 1

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
