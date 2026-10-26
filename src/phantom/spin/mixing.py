"""Singlet–triplet mixing of a separated pair: the spin funnel (ADR 0022).

A singlet separated into (1,1) sits J below the triplets' centre; in a field B the triplet T₋ sits at
−E_Z with E_Z = ḡ μ_B B. Where J(ε) = E_Z the singlet and T₋ cross, and any term that does not conserve
total spin (the hyperfine gradient between the dots, spin–orbit) opens an anticrossing of size Δ and
mixes them. Sweeping ε and B maps out the curve J(ε) = E_Z: the spin funnel. At the same time the
Zeeman difference ΔE_Z = (g₁ − g₂) μ_B B mixes S with T₀, which costs singlet probability wherever
ΔE_Z is comparable to J. Both are two-level mixings; with a wait τ the transferred probability is
Δ²/(Δ² + δ²) sin²(√(Δ² + δ²) τ / 2ħ), and averaged over a jittered wait it is ½ Δ²/(Δ² + δ²).
"""

from __future__ import annotations

import numpy as np

from phantom.constants import HBAR_EV_S

MODEL_VERSION = 1


def two_level_transfer(delta_ev: np.ndarray, coupling_ev: float, tau_s: float | None) -> np.ndarray:
    """Probability transferred across an anticrossing of coupling Δ at energy distance δ after a wait τ;
    the time-averaged ½ Δ²/(Δ² + δ²) when τ is None."""
    delta = np.asarray(delta_ev, dtype=float)
    w2 = coupling_ev**2 + delta**2
    lorentz = coupling_ev**2 / np.maximum(w2, 1e-300)
    if tau_s is None:
        return 0.5 * lorentz
    return lorentz * np.sin(np.sqrt(w2) * tau_s / (2 * HBAR_EV_S)) ** 2


def funnel_singlet_probability(
    j_ev: np.ndarray,
    e_zeeman_ev: np.ndarray,
    delta_mix_ev: float,
    delta_zeeman_ev: np.ndarray | float = 0.0,
    tau_s: float | None = None,
) -> np.ndarray:
    """Singlet probability after separation: 1 − P(S→T₋) − P(S→T₀), with the S–T₋ crossing at J = E_Z
    (coupling Δ_mix) and the S–T₀ mixing at distance J (coupling ΔE_Z, the Zeeman difference)."""
    p_tm = two_level_transfer(np.asarray(j_ev) - np.asarray(e_zeeman_ev), delta_mix_ev, tau_s)
    dz = np.asarray(delta_zeeman_ev, dtype=float)
    # S–T₀ mixing: coupling ΔE_Z/2 in the (S, T₀) basis at energy distance J
    w2 = (dz / 2) ** 2 + np.asarray(j_ev) ** 2
    lorentz = (dz / 2) ** 2 / np.maximum(w2, 1e-300)
    p_t0 = 0.5 * lorentz if tau_s is None else lorentz * np.sin(np.sqrt(w2) * tau_s / (2 * HBAR_EV_S)) ** 2
    return np.clip(1.0 - p_tm - p_t0, 0.0, 1.0)
