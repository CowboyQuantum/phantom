"""Dots and barriers as features of the potential-energy landscape.

Given a gate basis and a gate vector, the electron potential energy U(x, y) has minima where dots form
and saddles where barriers sit. This module finds them and turns them into the numbers the charge
layer consumes (ADR 0013):

- dot positions: the deepest local minima of U inside the layout's channel region (and sensor region);
- dot size: the harmonic radius ℓ = (ħ² / (m* k))^{1/4} from the curvature k of U at the minimum
  (the ground-state extent of a 2D harmonic oscillator), which `charge.Dot` turns into a self-capacitance;
- lever arms: the basis values at the minimum, dimensionless potential per volt, equal to the
  derivative of the dot's energy with respect to each gate voltage;
- barrier height between two dots: the maximum of U along the straight line between their minima,
  measured from the higher of the two minima.

All of this is compile-tier (a grid search); the runtime consumes the resulting parameters.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phantom.electrostatics.closed_form import GateBasis

MODEL_VERSION = 1
"""Bumped whenever any output of this module changes for the same inputs (ADR 0003)."""

HBAR2_OVER_2ME_EV_NM2 = 0.0380998212
"""ħ² / (2 m_e) in eV·nm² (CODATA 2018)."""


@dataclass(frozen=True)
class DerivedDot:
    name: str
    x_nm: float
    y_nm: float
    energy_ev: float
    """U at the minimum (eV)."""
    curvature_ev_per_nm2: float
    """Mean of the two principal curvatures of U at the minimum."""
    radius_nm: float
    """Harmonic ground-state radius ℓ."""
    lever_arms: np.ndarray
    """Dimensionless potential per volt of every gate at the minimum (gates,)."""


@dataclass(frozen=True)
class DerivedBarrier:
    between: tuple[str, str]
    height_ev: float
    """Barrier top minus the higher of the two dot minima; negative when no barrier separates them."""
    x_nm: float
    y_nm: float


def _region_mask(basis: GateBasis, region: tuple[float, float, float, float]) -> np.ndarray:
    x0, x1, y0, y1 = region
    return ((basis.x_nm >= x0) & (basis.x_nm <= x1))[None, :] & ((basis.y_nm >= y0) & (basis.y_nm <= y1))[:, None]


def local_minima(u: np.ndarray, mask: np.ndarray) -> list[tuple[int, int]]:
    """Grid points inside `mask` lower than their eight neighbours, deepest first, as (iy, ix)."""
    c = u[1:-1, 1:-1]
    is_min = np.ones_like(c, dtype=bool)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy == 0 and dx == 0:
                continue
            is_min &= c < u[1 + dy : u.shape[0] - 1 + dy, 1 + dx : u.shape[1] - 1 + dx]
    is_min &= mask[1:-1, 1:-1]
    iy, ix = np.nonzero(is_min)
    order = np.argsort(c[iy, ix])
    return [(int(iy[k]) + 1, int(ix[k]) + 1) for k in order]


def harmonic_radius_nm(curvature_ev_per_nm2: float, effective_mass_me: float) -> float:
    """ℓ = (ħ² / (m* k))^{1/4} with k the curvature of U (eV/nm²)."""
    k = max(curvature_ev_per_nm2, 1e-12)
    return float((2 * HBAR2_OVER_2ME_EV_NM2 / (effective_mass_me * k)) ** 0.25)


def find_dots(
    basis: GateBasis,
    voltages_v: np.ndarray,
    region: tuple[float, float, float, float],
    n_dots: int,
    effective_mass_me: float,
    names: tuple[str, ...] | None = None,
    min_separation_nm: float = 20.0,
) -> list[DerivedDot]:
    """The `n_dots` deepest, mutually separated minima of U inside `region`, ordered by x."""
    u = basis.energy_ev(voltages_v)
    dx = float(basis.x_nm[1] - basis.x_nm[0])
    dy = float(basis.y_nm[1] - basis.y_nm[0])
    picked: list[tuple[int, int]] = []
    for iy, ix in local_minima(u, _region_mask(basis, region)):
        if all(np.hypot((ix - jx) * dx, (iy - jy) * dy) >= min_separation_nm for jy, jx in picked):
            picked.append((iy, ix))
        if len(picked) == n_dots:
            break
    if len(picked) < n_dots:
        raise ValueError(f"found {len(picked)} minima in {region}, need {n_dots}; adjust the gate voltages")
    picked.sort(key=lambda p: p[1])
    dots = []
    for k, (iy, ix) in enumerate(picked):
        kxx = (u[iy, ix + 1] - 2 * u[iy, ix] + u[iy, ix - 1]) / dx**2
        kyy = (u[iy + 1, ix] - 2 * u[iy, ix] + u[iy - 1, ix]) / dy**2
        curv = float((kxx + kyy) / 2)
        name = names[k] if names else f"dot{k + 1}"
        dots.append(
            DerivedDot(
                name,
                float(basis.x_nm[ix]),
                float(basis.y_nm[iy]),
                float(u[iy, ix]),
                curv,
                harmonic_radius_nm(curv, effective_mass_me),
                basis.at(float(basis.x_nm[ix]), float(basis.y_nm[iy])),
            )
        )
    return dots


def barrier_between(basis: GateBasis, voltages_v: np.ndarray, a: DerivedDot, b: DerivedDot, samples: int = 201) -> DerivedBarrier:
    """Highest point of U on the straight line from dot `a` to dot `b`."""
    u = basis.energy_ev(voltages_v)
    xs = np.linspace(a.x_nm, b.x_nm, samples)
    ys = np.linspace(a.y_nm, b.y_nm, samples)
    ix = np.clip(np.rint((xs - basis.x_nm[0]) / (basis.x_nm[1] - basis.x_nm[0])).astype(int), 0, len(basis.x_nm) - 1)
    iy = np.clip(np.rint((ys - basis.y_nm[0]) / (basis.y_nm[1] - basis.y_nm[0])).astype(int), 0, len(basis.y_nm) - 1)
    line = u[iy, ix]
    k = int(np.argmax(line))
    return DerivedBarrier((a.name, b.name), float(line[k] - max(a.energy_ev, b.energy_ev)), float(xs[k]), float(ys[k]))


def line_cut(basis: GateBasis, voltages_v: np.ndarray, x0: float, y0: float, x1: float, y1: float, samples: int = 301) -> tuple[np.ndarray, np.ndarray]:
    """U (eV) along the straight line from (x0, y0) to (x1, y1); returns (distance_nm, u_ev)."""
    u = basis.energy_ev(voltages_v)
    xs = np.linspace(x0, x1, samples)
    ys = np.linspace(y0, y1, samples)
    ix = np.clip(np.rint((xs - basis.x_nm[0]) / (basis.x_nm[1] - basis.x_nm[0])).astype(int), 0, len(basis.x_nm) - 1)
    iy = np.clip(np.rint((ys - basis.y_nm[0]) / (basis.y_nm[1] - basis.y_nm[0])).astype(int), 0, len(basis.y_nm) - 1)
    return np.hypot(xs - x0, ys - y0), u[iy, ix]
