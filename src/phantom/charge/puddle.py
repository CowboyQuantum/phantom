"""The dots grow with their occupancy: each dot is a classical puddle whose radius scales as N^(1/3), so its capacitances follow the gate-induced charge and the interdot barrier is seen from its edge (ADR 0035).

The closed-form dot of ADR 0013 is a harmonic minimum with its one-electron radius, and the constant-
interaction model adds electrons without changing anything geometric. In a real dot the electrons form
a puddle. For a two-dimensional puddle of charge in a parabolic well the classical equilibrium density
is semi-elliptic, σ(r) = σ₀√(1 − r²/R²), because that density produces a parabolic in-plane potential
that cancels the well inside the puddle; the total charge is then N = (2π/3)σ₀R² with σ₀ ∝ kR/e², so

    N ∝ R³,   R ∝ N^(1/3),   e²/C ∝ 1/R ∝ N^(−1/3),   N ∝ φ^(3/2)

(the last from N = Cφ/e with C ∝ R). The classical density of a parabolic dot is Shikin, Nazin, Heitmann
and Demel's (Phys. Rev. B 43, 11903, 1991; `proposed` entry `shikin1991-parabolic-dot-classical-density`).
phantom uses the scaling, anchored at the compiled dot:

- **Scale.** At a gate point the gate-induced charge of the geometric model is n (the compiled lever arms
  and capacitances, the shape grid). Each dot's puddle has s_i = √(max(n_i, n_floor)/N_ref,i), with N_ref
  its gate-induced charge at the operating point: C ∝ R ∝ N^(1/3) and N = Cφ/e give C ∝ √n. The induced
  charge with the scaled capacitances is then the puddle's own N = s·n ∝ n^(3/2), the classical relation,
  and it equals n at the operating point.
- **Capacitances.** Each dot's self-capacitance (its total minus the mutual one, ADR 0013) is multiplied
  by s_i; the mutual capacitance, the sphere formula's 4πε r₁r₂/d, by s₁s₂. So the addition energy falls
  as N^(−1/3) and the geometric fractional splitting f_C = 2C_m/(C_t + C_m) grows about as s: cells are
  large and nearly a grid near (1,1), smaller hexagons at high occupancy.
- **Radius.** Each puddle's radius at the operating point is the bare potential's turning point at the
  interdot tunnelling level on the cut toward the other dot: the puddle of N_ref electrons fills the
  bare allowed region. Elsewhere R_i = s_i·R_ref,i. (`State.dot_radii_nm` keeps the harmonic radius of
  the shape grid; `State.puddle_radii_nm` is this one.)
- **The barrier seen from the edge.** An electron at the Fermi level sees the barrier from the puddle's
  edge. Along the compiled cut (the bare potential at the operating point, above the level), a puddle
  smaller than the bare allowed region sees the profile from further in, where it is lower, so the
  barrier top above the edge is higher by u(R_ref) − u(R): that is a rigid shift of the profile outside
  the puddle, as the flat self-consistent potential of a puddle meets the level at its edge. A larger
  puddle reaches toward the saddle and lowers it; at the saddle the barrier is gone. The interdot
  barrier (`charge.barrier`) adds the mean of the two dots' shifts to its top, so t, Γ_m and the
  conductance fall toward (1,1) and rise toward high occupancy at fixed gates. The shift is zero at the
  operating point and depends on the occupancy alone; the gates move the top as in ADR 0034.

Placeholders (`calibration: none`): the N^(1/3) scaling is classical and ignores the kinetic energy and
the gates' screening of the puddle; the compiled cut is held where the shape grid moves the dots; the
compiled dot's radius is the one-electron harmonic one (ADR 0013), here read as the puddle of N_ref
electrons.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np

MODEL_VERSION = 1


@dataclass(frozen=True)
class PuddleGrowth:
    """How the two dots' puddles scale with their gate-induced charge, anchored at the operating point."""

    n_ref: tuple[float, float]
    """Each dot's gate-induced charge at the operating point, electrons."""
    turning_nm: tuple[float, float]
    """Each puddle's radius at the operating point: the bare potential's turning point at the interdot
    tunnelling level on the compiled cut, measured from the dot's own centre toward the other, nm."""
    cut_nm: tuple[float, ...]
    """Distance from dot 1 along the compiled cut between the dots at the operating point, nm, ascending."""
    cut_ev: tuple[float, ...]
    """The bare potential above the interdot tunnelling level along the cut, eV (negative inside the dots)."""
    n_floor: float = 0.5
    """The smallest gate-induced charge the scaling is taken at, electrons: below it a puddle stops shrinking."""

    def __post_init__(self) -> None:
        if len(self.cut_nm) != len(self.cut_ev) or len(self.cut_nm) < 3:
            raise ValueError("the cut needs three or more samples with one energy each")
        if not np.all(np.diff(self.cut_nm) > 0):
            raise ValueError("the cut's distances must ascend")
        if not (all(n > 0 for n in self.n_ref) and all(r > 0 for r in self.turning_nm) and self.n_floor > 0):
            raise ValueError("n_ref, turning_nm and n_floor must be positive")

    @cached_property
    def _arrays(self) -> tuple[np.ndarray, np.ndarray]:
        return np.asarray(self.cut_nm, dtype=float), np.asarray(self.cut_ev, dtype=float)

    @property
    def length_nm(self) -> float:
        return float(self.cut_nm[-1])

    @cached_property
    def saddle_nm(self) -> float:
        """Where the cut is highest, from dot 1, nm."""
        return float(self.cut_nm[int(np.argmax(self._arrays[1]))])

    @cached_property
    def _reference_edges_ev(self) -> tuple[float, float]:
        """The cut's energy at each dot's turning point: the edge the shift is measured from."""
        return (float(self._edge_ev(0, np.asarray(self.turning_nm[0]))), float(self._edge_ev(1, np.asarray(self.turning_nm[1]))))

    def scale(self, n_fixed: np.ndarray) -> np.ndarray:
        """s (..., 2) = √(max(n, n_floor)/N_ref) from the geometric model's gate-induced charge (..., 2)."""
        n = np.maximum(np.asarray(n_fixed, dtype=float), self.n_floor)
        return np.sqrt(n / np.asarray(self.n_ref))

    def radii_nm(self, scale: np.ndarray) -> np.ndarray:
        """Each puddle's radius (..., 2), nm."""
        return np.asarray(scale, dtype=float) * np.asarray(self.turning_nm)

    def capacitances_af(self, c: np.ndarray, scale: np.ndarray) -> np.ndarray:
        """The CI matrix (..., 2, 2) with each dot's self-capacitance (its total minus the mutual one)
        multiplied by s_i and the mutual capacitance by s₁s₂."""
        c = np.asarray(c, dtype=float)
        s = np.asarray(scale, dtype=float)
        c_m = -c[..., 0, 1]
        both = s[..., 0] * s[..., 1]
        # written as changes from `c`, so that a scale of exactly 1 returns `c` bit for bit
        c11 = c[..., 0, 0] + (c[..., 0, 0] - c_m) * (s[..., 0] - 1.0) + c_m * (both - 1.0)
        c22 = c[..., 1, 1] + (c[..., 1, 1] - c_m) * (s[..., 1] - 1.0) + c_m * (both - 1.0)
        c_m = c_m * both
        return np.stack([np.stack([c11, -c_m], axis=-1), np.stack([-c_m, c22], axis=-1)], axis=-2)

    def _edge_ev(self, dot: int, radius_nm: np.ndarray) -> np.ndarray:
        """The cut's energy at `radius_nm` from dot `dot` toward the other, the radius clipped at the saddle."""
        xs, ys = self._arrays
        saddle = self.saddle_nm
        if dot == 0:
            r = np.clip(radius_nm, 0.0, saddle)
            return np.interp(r, xs, ys)
        r = np.clip(radius_nm, 0.0, self.length_nm - saddle)
        return np.interp(self.length_nm - r, xs, ys)

    def barrier_shift_ev(self, scale: np.ndarray) -> np.ndarray:
        """(...) how much higher the interdot barrier's top is above the puddles' edges than above the
        compiled turning points: the mean over the two dots of u(R_ref) − u(R), zero at the operating
        point, positive for smaller puddles, −ΔU_ref once both reach the saddle."""
        radii = self.radii_nm(scale)
        ref = self._reference_edges_ev
        return 0.5 * ((ref[0] - self._edge_ev(0, radii[..., 0])) + (ref[1] - self._edge_ev(1, radii[..., 1])))

    def to_dict(self) -> dict:
        return {
            "n_ref": list(self.n_ref),
            "turning_nm": list(self.turning_nm),
            "cut_nm": list(self.cut_nm),
            "cut_ev": list(self.cut_ev),
            "n_floor": self.n_floor,
        }

    @classmethod
    def from_dict(cls, d: dict) -> PuddleGrowth:
        return cls(
            (float(d["n_ref"][0]), float(d["n_ref"][1])),
            (float(d["turning_nm"][0]), float(d["turning_nm"][1])),
            tuple(float(x) for x in d["cut_nm"]),
            tuple(float(x) for x in d["cut_ev"]),
            float(d.get("n_floor", 0.5)),
        )


__all__ = ["MODEL_VERSION", "PuddleGrowth"]
