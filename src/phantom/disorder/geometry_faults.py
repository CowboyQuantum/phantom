"""Gate misalignment and width error as transforms of the layout polygons.

Catalog entry `gate-misalignment` (a `defect`): a lithographic layer is offset from its design by a
common shift, and each gate's edges move by an independent width error. Because the closed-form
electrostatics takes polygons as input, these faults are transforms of the layout, and the perturbed
lever-arm basis follows without any new physics: the dot is not under the gate the design says.

Draws are seeded and scale with `severity` (multiplies both standard deviations); severity 0 returns
the layout unchanged, so the clean device is reproduced bit for bit.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phantom.layout.geometry import Gate, Layout

MODEL_VERSION = 1


@dataclass(frozen=True)
class MisalignmentParams:
    """Placeholders until the catalog entry is written (`calibration: none`): a few nm of overlay error
    and of edge (critical-dimension) error, typical of electron-beam lithography at this scale."""

    shift_sigma_nm: float = 5.0
    """Standard deviation of the common offset of the whole gate layer, each axis."""
    width_sigma_nm: float = 2.0
    """Standard deviation of each gate's edge displacement (outward positive), per gate."""

    def __post_init__(self) -> None:
        if self.shift_sigma_nm < 0 or self.width_sigma_nm < 0:
            raise ValueError("sigmas must be non-negative")


@dataclass(frozen=True)
class MisalignmentDraw:
    shift_nm: tuple[float, float]
    edge_nm: dict[str, float]
    """Outward displacement of every edge of each gate (a width error of twice this)."""

    def to_dict(self) -> dict:
        return {"shift_nm": list(self.shift_nm), "edge_nm": dict(self.edge_nm)}

    @classmethod
    def from_dict(cls, d: dict) -> MisalignmentDraw:
        return cls((float(d["shift_nm"][0]), float(d["shift_nm"][1])), {k: float(v) for k, v in d["edge_nm"].items()})


def draw_misalignment(params: MisalignmentParams, layout: Layout, rng: np.random.Generator, severity: float = 1.0) -> MisalignmentDraw:
    if severity <= 0:
        return MisalignmentDraw((0.0, 0.0), {g.name: 0.0 for g in layout.gates})
    shift = rng.normal(0.0, severity * params.shift_sigma_nm, 2)
    edges = rng.normal(0.0, severity * params.width_sigma_nm, len(layout.gates))
    return MisalignmentDraw((float(shift[0]), float(shift[1])), {g.name: float(e) for g, e in zip(layout.gates, edges, strict=True)})


def _dilate(vertices: np.ndarray, by_nm: float) -> np.ndarray:
    """Move every vertex outward by `by_nm` along each axis according to its side of the centroid:
    each edge of an axis-aligned rectangle moves out by exactly `by_nm`; for other convex polygons it is
    the same per-axis edge displacement."""
    if by_nm == 0.0:
        return vertices.copy()
    c = vertices.mean(axis=0)
    return vertices + by_nm * np.sign(vertices - c)


def apply_misalignment(layout: Layout, draw: MisalignmentDraw) -> Layout:
    """The layout with the common shift and per-gate edge displacements applied. Severity 0 (all zeros)
    returns an equal layout."""
    if draw.shift_nm == (0.0, 0.0) and all(v == 0.0 for v in draw.edge_nm.values()):
        return layout
    gates = []
    for g in layout.gates:
        v = _dilate(g.vertices, draw.edge_nm.get(g.name, 0.0)) + np.array(draw.shift_nm)
        gates.append(Gate(g.name, tuple((float(x), float(y)) for x, y in v), g.role))
    return Layout(tuple(gates), layout.stack, layout.material, layout.channel_nm, layout.sensor_nm, layout.margin_nm)
