"""Closed-form potential of polygonal surface gates over a two-dimensional electron gas.

Davies, Larkin and Stern (J. Appl. Phys. 77, 4504, 1995; reference entry `davies1995-patterned-2deg`)
solve Laplace's equation in a half-space of uniform dielectric with the surface held at the gate
voltages and at zero between gates (Fermi-level pinning). The potential at depth d below a gate of
polygon P at voltage V is then V times the solid angle the polygon subtends at the point, over 2π:

    φ(x, y, d) = V · Ω_P(x, y, d) / 2π

Their rectangle formula is this solid angle written out in arctangents; the solid angle form holds
for any planar polygon and is what is implemented here, as a fan of triangles with the Van Oosterom
and Strackee expression per triangle. Superposition over gates gives the full potential; each gate's
Ω/2π on a grid is its *basis*, dimensionless potential per volt, which is also the electrostatic lever
arm of that gate at every point. Screening by the electrons themselves is not included (the formula
is for the depleted regions that define dots and barriers); Step 7's solver adds it (ADR 0001, 0013).

Electron potential energy is U = −e φ, so in eV it is numerically −φ in volts.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phantom.layout.geometry import Layout

MODEL_VERSION = 1
"""Bumped whenever any output of this module changes for the same inputs (ADR 0003)."""


def triangle_solid_angle(r1: np.ndarray, r2: np.ndarray, r3: np.ndarray) -> np.ndarray:
    """Signed solid angle (sr) of the triangle with corner vectors r1, r2, r3 (..., 3) from the origin
    (Van Oosterom and Strackee, IEEE Trans. Biomed. Eng. 30, 125, 1983)."""
    n1, n2, n3 = (np.linalg.norm(r, axis=-1) for r in (r1, r2, r3))
    num = np.einsum("...i,...i->...", r1, np.cross(r2, r3))
    den = n1 * n2 * n3 + np.einsum("...i,...i->...", r1, r2) * n3 + np.einsum("...i,...i->...", r2, r3) * n1 + np.einsum("...i,...i->...", r3, r1) * n2
    return 2.0 * np.arctan2(num, den)


def polygon_solid_angle(vertices_nm: np.ndarray, x_nm: np.ndarray, y_nm: np.ndarray, depth_nm: float) -> np.ndarray:
    """Solid angle (sr, non-negative) of a planar polygon at z = 0 as seen from the points
    (x, y, −depth); `x_nm` and `y_nm` broadcast against each other."""
    x = np.asarray(x_nm, dtype=float)
    y = np.asarray(y_nm, dtype=float)
    x, y = np.broadcast_arrays(x, y)
    v = np.asarray(vertices_nm, dtype=float)
    pts = np.stack([x, y, np.full(x.shape, -float(depth_nm))], axis=-1)  # (..., 3)
    corners = np.concatenate([v, np.zeros((len(v), 1))], axis=1)  # (n, 3) in the gate plane
    total = np.zeros(x.shape)
    r0 = corners[0] - pts
    for i in range(1, len(v) - 1):
        total = total + triangle_solid_angle(r0, corners[i] - pts, corners[i + 1] - pts)
    return np.abs(total)


def rectangle_potential_dls(x: np.ndarray, y: np.ndarray, depth: float, left: float, right: float, bottom: float, top: float) -> np.ndarray:
    """The rectangle result in closed form, potential per volt of gate: the solid angle of a rectangle
    is the signed sum over its four corners of the quadrant solid angle arctan(u v / (d √(u² + v² + d²))),
    with u and v the signed distances from the point to the two edges meeting at that corner. This is
    the Davies–Larkin–Stern rectangle formula (their eq. for a rectangular gate; the human confirms the
    printed form when the reference entry is confirmed) and an independent check of the general
    polygon code."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    def quadrant(u, v):
        return np.arctan(u * v / (depth * np.sqrt(u**2 + v**2 + depth**2)))

    total = quadrant(x - left, y - bottom) + quadrant(x - left, top - y) + quadrant(right - x, y - bottom) + quadrant(right - x, top - y)
    return total / (2 * np.pi)


@dataclass(frozen=True)
class GateBasis:
    """Per-gate potential per volt on a grid: `basis[g, iy, ix] = Ω_g / 2π` at (x[ix], y[iy], depth)."""

    gate_names: tuple[str, ...]
    x_nm: np.ndarray
    y_nm: np.ndarray
    depth_nm: float
    basis: np.ndarray  # (gates, ny, nx)
    offset_ev: np.ndarray | None = None
    """Gate-independent electron energy added to `energy_ev` (ny, nx): disorder (Step 3). None = clean."""

    @classmethod
    def from_layout(cls, layout: Layout, pitch_nm: float = 2.0, grid: tuple[np.ndarray, np.ndarray] | None = None) -> GateBasis:
        """Evaluate every gate's solid angle on a grid covering the layout's bounding box at the
        stack's effective depth, or on the given `grid` = (x_nm, y_nm) so that perturbed layouts share
        the clean one's grid. Compile-tier work: milliseconds for a few gates on a 300×300 grid."""
        if grid is not None:
            x, y = (np.asarray(grid[0], dtype=float), np.asarray(grid[1], dtype=float))
        else:
            x0, x1, y0, y1 = layout.bounding_box_nm()
            x = np.arange(x0, x1 + pitch_nm / 2, pitch_nm)
            y = np.arange(y0, y1 + pitch_nm / 2, pitch_nm)
        depth = layout.stack.effective_depth_nm
        basis = np.stack([polygon_solid_angle(g.vertices, x[None, :], y[:, None], depth) / (2 * np.pi) for g in layout.gates])
        return cls(layout.gate_names, x, y, depth, basis)

    def potential_v(self, voltages_v: np.ndarray) -> np.ndarray:
        """φ(x, y) in volts for gate voltages (gates,) or (..., gates); returns (..., ny, nx)."""
        return np.tensordot(np.asarray(voltages_v, dtype=float), self.basis, axes=([-1], [0]))

    def energy_ev(self, voltages_v: np.ndarray) -> np.ndarray:
        """Electron potential energy U = −eφ in eV plus any disorder offset, same shape as `potential_v`."""
        u = -self.potential_v(voltages_v)
        return u if self.offset_ev is None else u + self.offset_ev

    def at(self, x_nm: float, y_nm: float) -> np.ndarray:
        """Basis of every gate at one point by bilinear interpolation: the lever-arm vector there (gates,)."""
        ix = np.clip(np.searchsorted(self.x_nm, x_nm) - 1, 0, len(self.x_nm) - 2)
        iy = np.clip(np.searchsorted(self.y_nm, y_nm) - 1, 0, len(self.y_nm) - 2)
        fx = (x_nm - self.x_nm[ix]) / (self.x_nm[ix + 1] - self.x_nm[ix])
        fy = (y_nm - self.y_nm[iy]) / (self.y_nm[iy + 1] - self.y_nm[iy])
        b = self.basis
        return (
            b[:, iy, ix] * (1 - fx) * (1 - fy)
            + b[:, iy, ix + 1] * fx * (1 - fy)
            + b[:, iy + 1, ix] * (1 - fx) * fy
            + b[:, iy + 1, ix + 1] * fx * fy
        )

    def index_of(self, name: str) -> int:
        return self.gate_names.index(name)
