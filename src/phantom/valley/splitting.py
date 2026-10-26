"""Valley splitting per dot as a statistical draw (ADR 0020).

In a silicon quantum well the two low-lying conduction valleys are coupled by the sharp interfaces; the
coupling Δ is a complex number and the valley splitting is E_V = 2|Δ| (here written as |Δ| with the
factor absorbed). With random alloy disorder at a Si/SiGe interface the coupling is the sum of a
deterministic part Δ_det, set by the ideal interface and the vertical field, and a random part whose two
components are Gaussian with standard deviation σ. The splitting then follows a Rice distribution
(Rayleigh when Δ_det = 0), and the theory states how σ scales:

- σ falls as 1/ℓ with the dot radius ℓ, because the random coupling is a sum over the atoms the
  envelope covers (variance ∝ 1/area);
- σ grows as √(x̄(1 − x̄)) with the germanium concentration x̄ the electron samples, the binomial
  variance of the alloy, times the envelope's overlap with the alloy, which grows with the interface
  width;
- a wide interface or a tail into the barrier raises σ; a sharp interface raises Δ_det.

Silicon MOS has no alloy: its coupling is dominated by the deterministic interface term and grows with
the vertical field, with a smaller random part from interface roughness. Both are placeholders
(`calibration: none`) until the reference entries (`paqueletwuetz2022-valley-alloy`, `losert2023-valley-
strategies`, `yang2013-mos-valley`) are confirmed; the P1 row of Step 5 then sets σ and Δ_det.

An interface step inside the dot shifts the valley phase; the hook `step_factor` lets a later model
scale σ for a given step density without changing the draw's shape.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import i0e

MODEL_VERSION = 1


@dataclass(frozen=True)
class ValleyParams:
    """Parameters of the Rice distribution for a material class. Energies in µeV, lengths in nm."""

    material: str = "SiGe"
    delta_det_uev: float = 20.0
    """Deterministic coupling |Δ_det| (µeV): small for a graded SiGe interface, hundreds for MOS."""
    sigma_ref_uev: float = 45.0
    """Random component's standard deviation at the reference radius, concentration and interface width."""
    radius_ref_nm: float = 10.0
    ge_ref: float = 0.3
    """Reference germanium concentration x̄ in the barrier."""
    interface_width_ref_nm: float = 1.0
    field_v_per_nm: float = 0.0
    """MOS: vertical field; Δ_det grows linearly with it from `delta_det_uev` at zero field."""
    delta_det_per_field_uev_nm_per_v: float = 0.0

    def __post_init__(self) -> None:
        if self.sigma_ref_uev < 0 or self.delta_det_uev < 0 or self.radius_ref_nm <= 0 or not 0 < self.ge_ref < 1:
            raise ValueError("sigma and delta_det >= 0, radius_ref > 0, 0 < ge_ref < 1")

    @classmethod
    def for_material(cls, material: str) -> ValleyParams:
        if material == "MOS":
            return cls("MOS", delta_det_uev=250.0, sigma_ref_uev=40.0, ge_ref=0.3, field_v_per_nm=0.02, delta_det_per_field_uev_nm_per_v=1.0e4)
        return cls("SiGe")

    def sigma_uev(self, radius_nm: float, ge: float | None = None, interface_width_nm: float | None = None, step_factor: float = 1.0) -> float:
        """σ(ℓ, x̄, interface width): the stated scalings around the reference values."""
        ge = self.ge_ref if ge is None else ge
        width = self.interface_width_ref_nm if interface_width_nm is None else interface_width_nm
        alloy = np.sqrt(ge * (1 - ge) / (self.ge_ref * (1 - self.ge_ref)))
        overlap = np.sqrt(width / self.interface_width_ref_nm)
        return float(self.sigma_ref_uev * (self.radius_ref_nm / radius_nm) * alloy * overlap * step_factor)

    def delta_det_total_uev(self) -> float:
        return float(self.delta_det_uev + self.delta_det_per_field_uev_nm_per_v * self.field_v_per_nm)

    def mean_uev(self, radius_nm: float, **kwargs) -> float:
        """Mean of the Rice distribution (numerical), µeV."""
        return rice_mean(self.delta_det_total_uev(), self.sigma_uev(radius_nm, **kwargs))


def rice_mean(nu: float, sigma: float) -> float:
    """Mean of |ν + σ(n₁ + i n₂)|, ν ≥ 0, by quadrature; Rayleigh σ√(π/2) when ν = 0."""
    if sigma <= 0:
        return float(nu)
    r = np.linspace(0, nu + 12 * sigma, 20000)
    # Rice density (r/σ²) exp(-(r²+ν²)/2σ²) I0(rν/σ²), written with the exponentially scaled Bessel
    # function so that large ν/σ does not overflow
    pdf = (r / sigma**2) * np.exp(-((r - nu) ** 2) / (2 * sigma**2)) * i0e(r * nu / sigma**2)
    pdf /= np.trapezoid(pdf, r)
    return float(np.trapezoid(r * pdf, r))


def draw_valley_splitting_uev(params: ValleyParams, radius_nm: np.ndarray, rng: np.random.Generator, **kwargs) -> np.ndarray:
    """One splitting per dot radius in `radius_nm` (array), in µeV: |Δ_det + σ(n₁ + i n₂)|."""
    radii = np.asarray(radius_nm, dtype=float)
    sigma = np.array([params.sigma_uev(float(r), **kwargs) for r in radii.reshape(-1)]).reshape(radii.shape)
    nu = params.delta_det_total_uev()
    re = nu + sigma * rng.normal(0.0, 1.0, radii.shape)
    im = sigma * rng.normal(0.0, 1.0, radii.shape)
    return np.hypot(re, im)
