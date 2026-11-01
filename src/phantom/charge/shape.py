"""The double dot's shape along the barrier gates, compiled on a grid and interpolated at runtime (#59, ADR 0030).

The closed-form potential is linear in the gate voltages, but the dots it forms are not: a barrier gate
moves the minima, changes their curvature and so their size, the lever arms and the disorder offset at
them, the capacitances, and the actions of the barriers around them. Compile re-derives the double dot
on a grid of nodes along each barrier gate, the other gates at the operating point, and records at every
node the change of each quantity from the operating point. At runtime each gate's changes are
interpolated (monotone cubic, PCHIP) at its carried voltage and the gates' changes are added: the
effects of two barrier gates are taken as independent. Beyond a gate's last node the changes are held
there. Where the grid ends because the potential no longer forms two dots (a dot merged with its
neighbour or with its lead), a point beyond that end is flagged `merged`. One quantity is read
differently: the interdot barrier's top is the potential on the cut between the held dots, which the
gates keep moving, so `quantity_along` continues it linearly past the ends (`charge.barrier`, #70).

The quantities (`QUANTITIES`, per dot where indexed):

- `alpha` (2, gates): the dimensionless lever arms at the dot, dφ/dV;
- `offset_ev` (2,): the disorder's energy offset at the dot, eV;
- `c11_af`, `c22_af`, `cm_af`: the CI capacitance matrix [[C₁, −C_m], [−C_m, C₂]], aF;
- `x_nm`, `y_nm` (2,): the dot's position, nm, refined between the grid points of the basis;
- `radius_nm` (2,): the harmonic radius, nm;
- `log_rate` (3,): the exponents −2S of the barriers to the left lead, to the right lead and between
  the dots, at the tunnelling levels compile holds fixed (ADR 0027);
- `barrier_top_ev` (1,): the interdot barrier's top above its tunnelling level, eV (#66).

At a dot the potential is φ_i = α_i·V − o_i (U = −eφ + o). With the shape, the CI model's gate-induced
charge is n_g = C(V)·(α(V)·V − Δo(V)) / E_UNIT, so the honeycomb follows the moving dots. By the envelope
theorem this is the potential at the true minimum to first order in the plungers, which stay linear.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from scipy.interpolate import PchipInterpolator

MODEL_VERSION = 3
"""2: the interdot barrier's top on the grid (#66). 3: `quantity_along`, a quantity weighted per gate
and continued linearly past the grid's ends (#70); `change` and `merged` are unchanged."""


def quantity_layout(n_gates: int) -> dict[str, slice]:
    """Where each quantity sits in a node's row of length `row_length(n_gates)`."""
    sizes = (("alpha", 2 * n_gates), ("offset_ev", 2), ("c11_af", 1), ("c22_af", 1), ("cm_af", 1), ("x_nm", 2), ("y_nm", 2), ("radius_nm", 2), ("log_rate", 3), ("barrier_top_ev", 1))
    out, start = {}, 0
    for name, size in sizes:
        out[name] = slice(start, start + size)
        start += size
    return out


def row_length(n_gates: int) -> int:
    return max(s.stop for s in quantity_layout(n_gates).values())


QUANTITIES = tuple(quantity_layout(1))


@dataclass(frozen=True)
class ShapeGrid:
    """Per barrier gate, the nodes (carried voltages) and each node's changes from the operating point."""

    gate_indices: tuple[int, ...]
    """The barrier gates' positions in the gate vector."""
    n_gates: int
    """Length of the gate vector."""
    nodes_v: tuple[tuple[float, ...], ...]
    """Per barrier gate, its carried voltage at each node, ascending, V."""
    changes: tuple[tuple[tuple[float, ...], ...], ...]
    """Per barrier gate, (nodes, `row_length`) changes from the operating point; zero at its node there."""
    merged_below: tuple[bool, ...]
    """Per barrier gate, whether the grid ends below because the potential stops forming two dots."""
    merged_above: tuple[bool, ...]
    """Likewise above."""

    def __post_init__(self) -> None:
        if not (len(self.gate_indices) == len(self.nodes_v) == len(self.changes) == len(self.merged_below) == len(self.merged_above)):
            raise ValueError("one entry per barrier gate")
        for nodes, rows in zip(self.nodes_v, self.changes, strict=True):
            if len(nodes) < 2 or len(rows) != len(nodes) or any(len(r) != row_length(self.n_gates) for r in rows):
                raise ValueError("each gate needs two or more nodes and one full row per node")
            if not np.all(np.diff(nodes) > 0):
                raise ValueError("nodes must ascend")

    @cached_property
    def _interpolants(self) -> tuple[PchipInterpolator, ...]:
        return tuple(PchipInterpolator(np.asarray(n), np.asarray(c), axis=0, extrapolate=False) for n, c in zip(self.nodes_v, self.changes, strict=True))

    @cached_property
    def _end_slopes(self) -> tuple[tuple[np.ndarray, np.ndarray], ...]:
        """Per barrier gate, each interpolant's derivative row at its first and its last node."""
        out = []
        for nodes, f in zip(self.nodes_v, self._interpolants, strict=True):
            d = f.derivative()
            out.append((np.asarray(d(nodes[0])), np.asarray(d(nodes[-1]))))
        return tuple(out)

    @property
    def layout(self) -> dict[str, slice]:
        return quantity_layout(self.n_gates)

    def change(self, v: np.ndarray) -> np.ndarray:
        """(..., `row_length`) the summed changes at carried voltages `v` (..., gates), each gate's
        voltage held at its end nodes."""
        v = np.asarray(v, dtype=float)
        out = np.zeros(v.shape[:-1] + (row_length(self.n_gates),))
        for j, nodes, f in zip(self.gate_indices, self.nodes_v, self._interpolants, strict=True):
            out = out + f(np.clip(v[..., j], nodes[0], nodes[-1]))
        return out

    def quantity(self, change: np.ndarray, name: str) -> np.ndarray:
        """One quantity's changes (..., size) from a `change` row; `alpha` as (..., 2, gates)."""
        x = change[..., self.layout[name]]
        return x.reshape(x.shape[:-1] + (2, self.n_gates)) if name == "alpha" else x

    def quantity_along(self, v: np.ndarray, name: str, weights: np.ndarray | None = None, beyond: str = "hold") -> np.ndarray:
        """(..., size) one quantity's summed changes at carried voltages `v` (..., gates), each barrier
        gate's change multiplied by `weights[j]` (indexed as the gate vector; None is 1). Beyond a gate's
        end nodes its change is held (`beyond="hold"`, as `change`) or continued linearly with the
        interpolant's slope at that end (`"linear"`): for a quantity the gates keep moving where the
        dots are held, the potential at a point of the held dots' cut such as the interdot barrier's
        top (#70, ADR 0034). Positions, capacitances and lever arms stay held."""
        if beyond not in ("hold", "linear"):
            raise ValueError("beyond is 'hold' or 'linear'")
        v = np.asarray(v, dtype=float)
        sl = self.layout[name]
        out = np.zeros(v.shape[:-1] + (sl.stop - sl.start,))
        for j, nodes, f, (low, high) in zip(self.gate_indices, self.nodes_v, self._interpolants, self._end_slopes, strict=True):
            x = v[..., j]
            clipped = np.clip(x, nodes[0], nodes[-1])
            part = f(clipped)[..., sl]
            if beyond == "linear":
                part = part + np.where((x < nodes[0])[..., None], (x - nodes[0])[..., None] * low[sl], 0.0) + np.where((x > nodes[-1])[..., None], (x - nodes[-1])[..., None] * high[sl], 0.0)
            out = out + (part if weights is None else weights[j] * part)
        return out

    def merged(self, v: np.ndarray) -> np.ndarray:
        """(...) True where a barrier gate lies beyond an end of its grid at which the potential
        stops forming two dots."""
        v = np.asarray(v, dtype=float)
        out = np.zeros(v.shape[:-1], dtype=bool)
        for j, nodes, lo, hi in zip(self.gate_indices, self.nodes_v, self.merged_below, self.merged_above, strict=True):
            if lo:
                out |= v[..., j] < nodes[0]
            if hi:
                out |= v[..., j] > nodes[-1]
        return out

    def to_dict(self) -> dict:
        return {
            "gate_indices": list(self.gate_indices),
            "n_gates": self.n_gates,
            "nodes_v": [list(n) for n in self.nodes_v],
            "changes": [[list(r) for r in rows] for rows in self.changes],
            "merged_below": list(self.merged_below),
            "merged_above": list(self.merged_above),
        }

    @classmethod
    def from_dict(cls, d: dict) -> ShapeGrid:
        return cls(
            tuple(int(j) for j in d["gate_indices"]),
            int(d["n_gates"]),
            tuple(tuple(float(x) for x in n) for n in d["nodes_v"]),
            tuple(tuple(tuple(float(x) for x in r) for r in rows) for rows in d["changes"]),
            tuple(bool(x) for x in d["merged_below"]),
            tuple(bool(x) for x in d["merged_above"]),
        )


__all__ = ["MODEL_VERSION", "QUANTITIES", "ShapeGrid", "quantity_layout", "row_length"]
