"""Lifetime-broadened tunnelling through one level: occupation, linear-response lineshape and its wing.

A level coupled to two leads with tunnel rates Γ_S and Γ_D has a Lorentzian spectral function of full
width γ = ħ(Γ_S + Γ_D) (Breit–Wigner). The current through it at bias V_sd (Landauer, one spin channel)
is

    I = e Γ_S Γ_D / (Γ_S + Γ_D) · [ n_γ(ε − eV_sd/2) − n_γ(ε + eV_sd/2) ],

with ε the level relative to the leads' mean Fermi level and n_γ the Fermi function convolved with the
normalised Lorentzian, in closed form through the digamma function ψ:

    n_γ(ε) = ½ − Im ψ(½ + (γ/2 + iε) / (2π kT)) / π.

For γ ≪ kT this is the Fermi function and the current is the sequential-tunnelling result (reference
entry `beenakker1991-coulomb-blockade`); for γ ≫ kT it is the zero-temperature Lorentzian step and the
peak conductance reaches (e²/h)·4Γ_SΓ_D/(Γ_S + Γ_D)², e²/h for symmetric barriers. The convolution of
the thermal and lifetime widths is the lineshape used to fit Coulomb peaks of small dots (entry
`foxman1993-quantum-levels-coulomb-island`).

In linear response the current is e Γ_SΓ_D/(Γ_S + Γ_D) · eV_sd · L(ε) with the normalised lineshape

    L(ε) = −dn_γ/dε = Re ψ′(z) / (2π² kT),   z = ½ + (γ/2 + iε) / (2π kT),

whose slope involves ψ″ and whose inflection, the *wing* where the slope is steepest, is the zero of
Re ψ‴(z). scipy evaluates ψ for complex arguments but not its derivatives, so `polygamma` evaluates
ψ⁽ᵐ⁾(z) for m ≥ 1 by the recurrence ψ⁽ᵐ⁾(z) = ψ⁽ᵐ⁾(z + 1) + (−1)^{m+1} m!/z^{m+1}, shifting the argument
by `SHIFT`, and the asymptotic series with Bernoulli numbers there. The wing is found by bisection on
the sign of Re ψ‴, which is platform-independent to the last bits of the bracket (ADR 0014). The two
limits are closed forms: ε_wing = 2kT·atanh(1/√3) (the inflection of sech²) as γ → 0, and γ/(2√3) (the
inflection of the Lorentzian) as kT → 0.

Spin degeneracy is not counted (a factor of order one on the peak height and a shift of order kT·ln 2
of its position), as in the rest of the sensor model.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import psi

from phantom.constants import HBAR_EV_S

MODEL_VERSION = 1
"""Bumped whenever any output of this module changes for the same inputs (ADR 0003)."""

SHIFT = 20
"""Recurrence steps before the asymptotic series: |z + SHIFT| ≥ 20.5 for Re z ≥ ½."""

_BERNOULLI = (1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66, -691 / 2730, 7 / 6)
"""B₂, B₄, …, B₁₄."""


def polygamma(m: int, z: np.ndarray) -> np.ndarray:
    """ψ⁽ᵐ⁾(z) for m ≥ 1 and complex z with Re z > 0, to about 1e-15 relative."""
    if m < 1:
        raise ValueError("polygamma takes m >= 1; use scipy.special.psi for m = 0")
    z = np.asarray(z, dtype=complex)
    sign = -1.0 if m % 2 == 0 else 1.0  # (−1)^{m+1}
    fact_m = math.factorial(m)
    head = np.zeros(z.shape, dtype=complex)
    for k in range(SHIFT):
        head = head + 1.0 / (z + k) ** (m + 1)
    w = z + SHIFT
    tail = math.factorial(m - 1) / w**m + fact_m / (2 * w ** (m + 1))
    for k, b in enumerate(_BERNOULLI, start=1):
        tail = tail + b * math.factorial(2 * k + m - 1) / math.factorial(2 * k) / w ** (2 * k + m)
    return sign * (tail + fact_m * head)


def _z(eps_ev: np.ndarray, gamma_ev: np.ndarray, kt_ev: float) -> np.ndarray:
    return 0.5 + (0.5 * np.asarray(gamma_ev, dtype=float) + 1j * np.asarray(eps_ev, dtype=float)) / (2 * np.pi * kt_ev)


def occupation(eps_ev: np.ndarray, gamma_ev: np.ndarray, kt_ev: float) -> np.ndarray:
    """n_γ(ε): the occupation of a level of full width γ at energy ε above the Fermi level (eV)."""
    return 0.5 - psi(_z(eps_ev, gamma_ev, kt_ev)).imag / np.pi


def lineshape_per_ev(eps_ev: np.ndarray, gamma_ev: np.ndarray, kt_ev: float) -> np.ndarray:
    """L(ε) = −dn_γ/dε (1/eV): the normalised linear-response peak; 1/(4kT) at ε = 0 when γ = 0."""
    return polygamma(1, _z(eps_ev, gamma_ev, kt_ev)).real / (2 * np.pi**2 * kt_ev)


def lineshape_slope_per_ev2(eps_ev: np.ndarray, gamma_ev: np.ndarray, kt_ev: float) -> np.ndarray:
    """dL/dε (1/eV²): negative above the peak."""
    return -polygamma(2, _z(eps_ev, gamma_ev, kt_ev)).imag / (4 * np.pi**3 * kt_ev**2)


def wing_ev(gamma_ev: float, kt_ev: float, iterations: int = 200) -> float:
    """ε of the steepest point of L above the peak: the zero of Re ψ‴(z(ε)), by bisection."""
    lo, hi = 0.0, 4.0 * (2.0 * kt_ev + gamma_ev)

    def curvature_sign(eps: float) -> bool:  # True where L is concave (inside the wing)
        return bool(polygamma(3, _z(np.array(eps), gamma_ev, kt_ev)).real > 0)

    while curvature_sign(hi):
        hi *= 2.0
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:
            break
        if curvature_sign(mid):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def peak_conductance_quanta(gamma_s_hz: np.ndarray, gamma_d_hz: np.ndarray, kt_ev: float) -> np.ndarray:
    """Zero-bias peak conductance in units of e²/h: G/(e²/h) = 2πħ · Γ_SΓ_D/(Γ_S + Γ_D) · L(0), which
    is Γ_SΓ_D/(Γ_S + Γ_D) · h/(4kT) for γ ≪ kT and tends to 4Γ_SΓ_D/(Γ_S + Γ_D)² as γ/kT → ∞."""
    gs, gd = np.asarray(gamma_s_hz, dtype=float), np.asarray(gamma_d_hz, dtype=float)
    total = np.maximum(gs + gd, 1e-300)
    gamma_ev = HBAR_EV_S * (gs + gd)
    return 2 * np.pi * HBAR_EV_S * gs * gd / total * lineshape_per_ev(np.zeros_like(gamma_ev), gamma_ev, kt_ev)


__all__ = [
    "MODEL_VERSION",
    "lineshape_per_ev",
    "lineshape_slope_per_ev2",
    "occupation",
    "peak_conductance_quanta",
    "polygamma",
    "wing_ev",
]
