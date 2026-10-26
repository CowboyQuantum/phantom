"""Random fixed charges at an interface, as a perturbation of the electron energy landscape.

Catalog entry `interface-fixed-charge` (a `trait`): charged centres at the semiconductor–dielectric
interface with an areal density per material class. Each charge q = ±e sits at depth `h` below the
pinned surface; in the Davies–Larkin–Stern half-space (surface at zero potential, uniform dielectric,
ADR 0013) its potential at a point of the electron plane at depth `d` is that of the charge plus its
image above the surface:

    φ(ρ) = q / (4π ε₀ ε_r) · [ 1 / √(ρ² + (d − h)² + z₀²) − 1 / √(ρ² + (d + h)²) ]

with ρ the in-plane distance and z₀ a smoothing length standing in for the finite thickness of the
electron layer, which keeps a charge sitting in the plane itself (the MOS case, h = d) finite. The
electron energy perturbation is −e φ. Screening by the electrons is not included, as in the rest of
the closed-form electrostatics; the perturbation is therefore an upper bound, and the entry is
`calibration: none` until its density and sign fraction come from a confirmed reference.

Draws are seeded (ADR 0003) and scale with a `severity` knob that multiplies the density
(`severity = 0` draws no charges and leaves the landscape bit-for-bit unchanged).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phantom.constants import E_CHARGE_C, EPSILON_0_F_PER_M

MODEL_VERSION = 1
"""Bumped whenever any output of this module changes for the same parameters and seed (ADR 0003)."""


@dataclass(frozen=True)
class FixedChargeParams:
    """Areal density and sign statistics of interface charges. Defaults are the catalog draft's
    placeholders (MOS ~3e10 cm⁻², SiGe ~3e9 cm⁻²; mostly positive), `calibration: none`."""

    density_cm2: float = 3e9
    positive_fraction: float = 0.8
    """Fraction of charges that are positive (attractive for electrons)."""
    depth_nm: float = 7.0
    """Depth of the charged interface below the surface (SiGe: oxide/cap interface; MOS: the oxide/Si interface = the electron plane)."""
    smoothing_nm: float = 3.0
    """z₀: effective thickness of the electron layer, keeps in-plane charges finite."""

    def __post_init__(self) -> None:
        if self.density_cm2 < 0 or not 0 <= self.positive_fraction <= 1 or self.smoothing_nm <= 0:
            raise ValueError("density >= 0, 0 <= positive_fraction <= 1, smoothing > 0")


@dataclass(frozen=True)
class FixedChargeDraw:
    """One realisation: positions (nm) and signs (±1) of the charges inside a region."""

    x_nm: np.ndarray
    y_nm: np.ndarray
    sign: np.ndarray
    depth_nm: float
    smoothing_nm: float

    @property
    def count(self) -> int:
        return len(self.sign)

    def to_dict(self) -> dict:
        return {
            "x_nm": self.x_nm.tolist(),
            "y_nm": self.y_nm.tolist(),
            "sign": self.sign.tolist(),
            "depth_nm": self.depth_nm,
            "smoothing_nm": self.smoothing_nm,
        }

    @classmethod
    def from_dict(cls, d: dict) -> FixedChargeDraw:
        return cls(np.asarray(d["x_nm"], float), np.asarray(d["y_nm"], float), np.asarray(d["sign"], float), float(d["depth_nm"]), float(d["smoothing_nm"]))


def draw_fixed_charges(
    params: FixedChargeParams,
    region_nm: tuple[float, float, float, float],
    rng: np.random.Generator,
    severity: float = 1.0,
) -> FixedChargeDraw:
    """Poisson-distributed charges over `region_nm` = (x0, x1, y0, y1) at `severity × density`."""
    x0, x1, y0, y1 = region_nm
    area_cm2 = (x1 - x0) * (y1 - y0) * 1e-14
    n = int(rng.poisson(max(severity, 0.0) * params.density_cm2 * area_cm2)) if severity > 0 else 0
    x = rng.uniform(x0, x1, n)
    y = rng.uniform(y0, y1, n)
    sign = np.where(rng.random(n) < params.positive_fraction, 1.0, -1.0)
    return FixedChargeDraw(x, y, sign, params.depth_nm, params.smoothing_nm)


def charge_energy_ev(draw: FixedChargeDraw, x_nm: np.ndarray, y_nm: np.ndarray, electron_depth_nm: float, epsilon_r: float) -> np.ndarray:
    """Electron energy perturbation −eφ (eV) of all charges in `draw` on the grid (y, x)."""
    x = np.asarray(x_nm, float)[None, :]
    y = np.asarray(y_nm, float)[:, None]
    out = np.zeros(np.broadcast(x, y).shape)
    if draw.count == 0:
        return out
    prefactor_ev_nm = E_CHARGE_C / (4 * np.pi * EPSILON_0_F_PER_M * epsilon_r) * 1e9  # e/(4πε) in V·nm → eV·nm per e
    dz_direct2 = (electron_depth_nm - draw.depth_nm) ** 2 + draw.smoothing_nm**2
    dz_image2 = (electron_depth_nm + draw.depth_nm) ** 2
    for cx, cy, s in zip(draw.x_nm, draw.y_nm, draw.sign, strict=True):
        rho2 = (x - cx) ** 2 + (y - cy) ** 2
        phi = s * prefactor_ev_nm * (1.0 / np.sqrt(rho2 + dz_direct2) - 1.0 / np.sqrt(rho2 + dz_image2))
        out -= phi  # electron energy is −eφ
    return out
