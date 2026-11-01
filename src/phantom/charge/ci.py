"""Constant-interaction (CI) model of a few-electron double quantum dot.

Seeded from qdt-sim's `DoubleDotModel` (Step 1); the electrostatics only, no sensor and no noise.
Energies in eV, capacitances in aF, voltages in V, charges in units of the electron charge e.

For gate voltages `V` the gate-induced charge on the dots is `n_g = C_gate · V / E_UNIT`; the
electrostatic energy of the integer occupation `N` is

    U(N) = ½ · E_UNIT · (N − n_g)ᵀ · C_dot⁻¹ · (N − n_g)      [eV]

(gate-only terms dropped). At electron temperature `T` the occupation is Boltzmann-distributed over a
window of integer states around `n_g`. Everything is a pure function of the gate vector.

Derived quantities (all exact for this model, used by the P0 tests):
- addition energy of dot i: `E_UNIT · (C_dot⁻¹)_ii`;
- lever arm of gate j on dot i, in eV per V: `(C_dot⁻¹ · C_gate)_ij`;
- transition period along gate j for dot i: `E_UNIT / C_gate[i, j]` volts per electron;
- the two triple points of the (0,0)/(1,0)/(0,1)/(1,1) honeycomb vertex, in n_g, and their separation,
  which vanishes when the mutual capacitance vanishes.

With an interdot tunnel `coupling` (#58, ADR 0029) the window's charge states of one total charge
hybridize (`charge.hybrid`). The occupation, the ground state and the draws then come from the
eigenstates, and the energies and the derived quantities above stay those of the charge states.

With a `shape` (#59, ADR 0030) the capacitances, the lever arms and the offsets at the dots follow the
barrier gates (`charge.shape`): at each gate point C_dot(V) and n_g = C_dot(V)·(α(V)·V − Δo(V))/E_UNIT. The
fields `c_dot_af` and `c_gate_af` and the derived quantities above are those at the operating point.

With a `splitting` (#66, ADR 0032) the interdot conductance adds an effective interdot capacitance
ΔC_m(V)·[[1, −1], [−1, 1]] to the energies (`charge.splitting`), so the triple points separate from a
grid towards a superdot as the barrier opens; the gate-induced charge keeps the geometric capacitances.

With a `puddle` (#72, ADR 0035) each dot's capacitances scale with its gate-induced charge
(`charge.puddle`): the self-capacitance by s_i = √(n_i/N_ref,i), the mutual one by s₁s₂, and the induced
charge is that of the scaled capacitances, so the addition energy falls as N^(−1/3) and the cells shrink
towards high occupancy. The same scale is handed to the interdot barrier, whose top the puddles' edges
see higher at low occupancy and lower at high.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np

from phantom.charge.hybrid import (
    STATIC_WINDOW,
    InterdotCoupling,
    WindowEigenstates,
    boltzmann,
    window_eigenstates,
)
from phantom.charge.puddle import PuddleGrowth
from phantom.charge.shape import ShapeGrid
from phantom.charge.splitting import InterdotSplitting
from phantom.constants import E_UNIT_V_PER_E_AF, K_B_EV_PER_K

MODEL_VERSION = 3
"""Bumped whenever any output of this module changes for the same parameters (ADR 0003). 2: duplicated
clipped states near the empty corner carry no weight (qdt-sim over-counted them). 3: with a `puddle`
the capacitances and the induced charge scale with the occupancy, and the occupancy scale reaches the
coupling and the splitting (#72, ADR 0035). Without a `puddle` the output is that of version 2."""


@dataclass(frozen=True)
class ConstantInteractionParams:
    """Capacitance parameters of a two-dot constant-interaction model.

    Defaults are a plausible Si/SiGe-class device inherited from qdt-sim; they are not a measured device.
    Ranges: dot self-capacitances 20–200 aF, mutual 5–40 aF, gate-to-dot 2–30 aF, cross-capacitances a
    few aF, electron temperature 0.05–0.5 K.
    """

    gate_names: tuple[str, ...] = ("P1", "P2", "SG")
    """Gate voltages the model takes, in this order."""
    c_dot_af: np.ndarray = field(default_factory=lambda: np.array([[85.0, -20.0], [-20.0, 90.0]]))
    """Dot capacitance matrix [[C1, −Cm], [−Cm, C2]] in aF; Cm is the interdot mutual capacitance."""
    c_gate_af: np.ndarray = field(
        default_factory=lambda: np.array([[12.0, 2.5, 0.3], [2.5, 11.5, 0.3]])
    )
    """Gate-to-dot capacitances [dot, gate] in aF; off-diagonals are cross-capacitances."""
    temperature_k: float = 0.15
    """Electron temperature of the dots in K."""
    window: int = 2
    """Integer occupations considered: round(n_g) ± window per dot."""
    coupling: InterdotCoupling | None = None
    """The interdot tunnel coupling at the gate voltages; None keeps the charge states classical."""
    shape: ShapeGrid | None = None
    """The capacitances, lever arms and offsets along the barrier gates; None keeps them fixed."""
    splitting: InterdotSplitting | None = None
    """The interdot conductance and the triple-point splitting it sets (#66): an effective interdot
    capacitance added to the energies; None keeps the geometric capacitances."""
    puddle: PuddleGrowth | None = None
    """The dots' growth with their occupancy (#72): the capacitances scale with the gate-induced charge
    and the interdot barrier is seen from the puddles' edges; None keeps the dots' size fixed."""

    def __post_init__(self) -> None:
        n_g = len(self.gate_names)
        if self.c_gate_af.shape != (2, n_g):
            raise ValueError(f"c_gate_af must be (2, {n_g}) for gates {self.gate_names}")
        if self.c_dot_af.shape != (2, 2):
            raise ValueError("c_dot_af must be (2, 2)")
        if self.window < 1:
            raise ValueError("window must be >= 1")

    # --- electrostatics ------------------------------------------------------------------------

    def induced_charge(self, v: np.ndarray) -> np.ndarray:
        """Gate-induced charge on (dot 1, dot 2) in electrons for gate voltages `v` (last axis = gates)."""
        if self.shape is None:
            return np.asarray(v, dtype=float) @ self.c_gate_af.T / E_UNIT_V_PER_E_AF
        return self._geometric(np.asarray(v, dtype=float))[1]

    def local(self, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(C_dot (..., 2, 2) in aF, n_g (..., 2) in electrons) at gate voltages `v` (..., gates): with a
        `shape`, the capacitances and the induced charge of the dots as the barrier gates have shaped them;
        with a `puddle`, scaled with the occupancy; with a `splitting`, the effective interdot capacitance
        added."""
        v = np.asarray(v, dtype=float)
        c, n_g, s = self._geometry(v)
        if self.splitting is not None:
            _, extra = self.splitting.splitting(v, c, s)
            c = c + extra[..., None, None] * np.array([[1.0, -1.0], [-1.0, 1.0]])
        return c, n_g

    def _bare(self, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(C_dot (..., 2, 2) in aF, φ (..., 2) in V) of the compiled dots, the potential at each dot from
        the lever arms; with a `shape`, as the barrier gates have shaped them."""
        if self.shape is None:
            alpha = np.linalg.solve(self.c_dot_af, self.c_gate_af)
            return np.broadcast_to(self.c_dot_af, v.shape[:-1] + (2, 2)), v @ alpha.T
        g = self.shape
        change = g.change(v)
        c11, c22, cm = (g.quantity(change, k)[..., 0] for k in ("c11_af", "c22_af", "cm_af"))
        c = self.c_dot_af + np.stack([np.stack([c11, -cm], axis=-1), np.stack([-cm, c22], axis=-1)], axis=-2)
        alpha = np.linalg.solve(self.c_dot_af, self.c_gate_af) + g.quantity(change, "alpha")
        return c, np.einsum("...ig,...g->...i", alpha, v) - g.quantity(change, "offset_ev")

    def _geometry(self, v: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
        """(C_dot (..., 2, 2), n_g (..., 2), the puddles' scale (..., 2) or None) from the geometry alone:
        the shape grid and the puddles, no splitting."""
        if self.shape is None and self.puddle is None:
            return np.broadcast_to(self.c_dot_af, v.shape[:-1] + (2, 2)), v @ self.c_gate_af.T / E_UNIT_V_PER_E_AF, None
        c, phi = self._bare(v)
        n_g = np.einsum("...ij,...j->...i", c, phi) / E_UNIT_V_PER_E_AF
        if self.puddle is None:
            return c, n_g, None
        s = self.puddle.scale(n_g)
        c = self.puddle.capacitances_af(c, s)
        return c, np.einsum("...ij,...j->...i", c, phi) / E_UNIT_V_PER_E_AF, s

    def _geometric(self, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(C_dot (..., 2, 2), n_g (..., 2)) from the geometry alone: the shape grid and the puddles, no splitting."""
        c, n_g, _ = self._geometry(v)
        return c, n_g

    def occupancy_scale(self, v: np.ndarray) -> np.ndarray:
        """s (..., 2) at gate voltages `v`: each puddle's √(n/N_ref) (`charge.puddle`), the scale the
        interdot barrier takes; ones without a `puddle`."""
        v = np.asarray(v, dtype=float)
        if self.puddle is None:
            return np.ones(v.shape[:-1] + (2,))
        return self._geometry(v)[2]

    @property
    def varies(self) -> bool:
        """Whether the capacitances depend on the gate voltages (a `shape`, a `puddle` or a `splitting`)."""
        return self.shape is not None or self.puddle is not None or self.splitting is not None

    def _c_inv(self, v: np.ndarray) -> np.ndarray:
        """C_dot⁻¹: (2, 2) when the capacitances are fixed, (..., 2, 2) at gate voltages `v` otherwise."""
        if not self.varies:
            return np.linalg.inv(self.c_dot_af)
        return np.linalg.inv(self.local(v)[0])

    def interdot_splitting(self, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(g, f) (...) at gate voltages `v`: the interdot conductance per channel and the fractional
        triple-point splitting; without a `splitting`, (0, the geometric f_C = 2C_m/(C_t + C_m))."""
        v = np.asarray(v, dtype=float)
        c, _, s = self._geometry(v)
        if self.splitting is None:
            c_t, c_m = 0.5 * (c[..., 0, 0] + c[..., 1, 1]), -c[..., 0, 1]
            return np.zeros(v.shape[:-1]), 2 * c_m / (c_t + c_m)
        f, _ = self.splitting.splitting(v, c, s)
        return self.splitting.conductance(v, s), f

    def addition_energies_at_ev(self, v: np.ndarray) -> np.ndarray:
        """(..., 2) each dot's addition energy E_UNIT·(C_dot⁻¹)_ii at gate voltages `v`, eV."""
        v = np.asarray(v, dtype=float)
        c_inv = np.broadcast_to(self._c_inv(v), v.shape[:-1] + (2, 2))
        return E_UNIT_V_PER_E_AF * np.stack([c_inv[..., 0, 0], c_inv[..., 1, 1]], axis=-1)

    def _quadratic(self, d: np.ndarray, c_inv: np.ndarray | None) -> np.ndarray:
        """½ E_UNIT dᵀ C⁻¹ d over the last axis of `d` (..., k, 2), C⁻¹ (2, 2) or one per point (..., 2, 2)."""
        if c_inv is None or c_inv.ndim == 2:
            c_inv = np.linalg.inv(self.c_dot_af) if c_inv is None else c_inv
            return 0.5 * E_UNIT_V_PER_E_AF * np.einsum("...ki,ij,...kj->...k", d, c_inv, d)
        return 0.5 * E_UNIT_V_PER_E_AF * np.einsum("...ki,...ij,...kj->...k", d, c_inv, d)

    def energies(self, n_g: np.ndarray, c_inv: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
        """Occupations `N` (..., k, 2) in the window around `n_g` (..., 2) and their energies (..., k) in eV;
        `c_inv` (..., 2, 2) is C_dot⁻¹ at each point when the capacitances vary (the operating point's when None)."""
        n_g = np.asarray(n_g, dtype=float)
        center = np.rint(n_g).astype(int)
        offsets = np.array(
            list(itertools.product(range(-self.window, self.window + 1), repeat=2)), dtype=int
        )
        states = np.maximum(center[..., None, :] + offsets, 0)
        d = states - n_g[..., None, :]
        u = self._quadratic(d, c_inv)
        # Clipping at zero electrons makes some window entries identical near the empty corner. Only the
        # first occurrence of a state keeps its energy; duplicates get +inf so they carry no weight
        # (MODEL_VERSION 2: qdt-sim counted them several times).
        key = states[..., 0] * (4 * self.window + 3) + states[..., 1]
        k = key.shape[-1]
        dup = np.zeros(key.shape, dtype=bool)
        for i in range(1, k):
            dup[..., i] = (key[..., :i] == key[..., i : i + 1]).any(axis=-1)
        u = np.where(dup, np.inf, u)
        return states, u

    def state_energies_ev(self, v: np.ndarray, states: np.ndarray) -> np.ndarray:
        """Energies (k,) of the given integer occupations `states` (k, 2) at gate voltages `v` (gates,),
        for any states, not only the window: U(N) = ½ E_UNIT (N − n_g)ᵀ C_dot⁻¹ (N − n_g)."""
        v = np.asarray(v, dtype=float)
        n_g = self.induced_charge(v)
        d = np.asarray(states, dtype=float) - n_g
        c_inv = self._c_inv(v)
        return 0.5 * E_UNIT_V_PER_E_AF * np.einsum("ki,ij,kj->k", d, c_inv, d)

    def occupation_probabilities(
        self, v: np.ndarray, kt_ev: float | np.ndarray | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """States (..., k, 2) and Boltzmann probabilities (..., k) at the model temperature, or at the
        energy width `kt_ev` (scalar or one per gate point) when given, e.g. a lifetime-broadened width.
        With a `coupling` the probabilities are those of the charge states in the thermal state of the
        eigenstates (the density matrix's diagonal), so their mean is the thermal ⟨N⟩."""
        if self.coupling is not None:
            eig = self.eigenstates(v)
            return eig.states, eig.populations(boltzmann(eig.energies_ev, self._kt(kt_ev)))
        states, u = self.energies(self.induced_charge(v), None if not self.varies else self._c_inv(v))
        if kt_ev is None:
            kt = K_B_EV_PER_K * max(self.temperature_k, 1e-6)
        else:
            kt = np.maximum(np.asarray(kt_ev, dtype=float), 1e-12)[..., None]
        w = np.exp(-(u - u.min(axis=-1, keepdims=True)) / kt)
        return states, w / w.sum(axis=-1, keepdims=True)

    def mean_occupation(self, v: np.ndarray, kt_ev: float | np.ndarray | None = None) -> np.ndarray:
        """Thermal-average occupation ⟨N⟩ (..., 2) in electrons for gate voltages `v`."""
        if self.coupling is not None:
            eig = self.eigenstates(v)
            return eig.mean_charge(boltzmann(eig.energies_ev, self._kt(kt_ev)))
        states, p = self.occupation_probabilities(v, kt_ev)
        return np.einsum("...k,...ki->...i", p, states)

    def ground_state_occupation(self, v: np.ndarray) -> np.ndarray:
        """The lowest-energy integer occupation (..., 2) for gate voltages `v` (the T → 0 charge state);
        with a `coupling`, the charge state with the largest weight in the ground eigenstate."""
        if self.coupling is not None:
            eig = self.eigenstates(v)
            return eig.dominant_state(np.argmin(eig.energies_ev, axis=-1))
        states, u = self.energies(self.induced_charge(v), None if not self.varies else self._c_inv(v))
        idx = np.argmin(u, axis=-1)
        return np.take_along_axis(states, idx[..., None, None], axis=-2)[..., 0, :]

    def sample_occupation(self, v: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        """One occupation (2,) drawn from the Boltzmann distribution at gate voltages `v` (1-D); with a
        `coupling`, one eigenstate drawn and its charge expectation (electrons, not an integer)."""
        if self.coupling is not None:
            eig = self.eigenstates(np.asarray(v, dtype=float))
            p = boltzmann(eig.energies_ev, self._kt(None))
            return eig.charge()[rng.choice(len(p), p=p)]
        states, p = self.occupation_probabilities(np.asarray(v, dtype=float))
        return states[rng.choice(len(p), p=p)]

    # --- interdot tunnelling (#58) ------------------------------------------------------------------

    def tunnel_coupling_ev(self, v: np.ndarray) -> np.ndarray:
        """Interdot t (...) in eV at gate voltages `v` (..., gates); zero without a `coupling`."""
        v = np.asarray(v, dtype=float)
        if self.coupling is None:
            return np.zeros(v.shape[:-1])
        return self.coupling.tunnel_coupling_ev(v, None if self.puddle is None else self.occupancy_scale(v))

    def eigenstates(self, v: np.ndarray, level_shift_ev: np.ndarray | None = None) -> WindowEigenstates:
        """The window's hybridized eigenstates at gate voltages `v` (..., gates) (`charge.hybrid`). The last
        result is kept, so the occupation, the ground state and the sensor at one batch of points
        diagonalise it once."""
        v = np.asarray(v, dtype=float)
        if level_shift_ev is None:
            last = self.__dict__.get("_last_eigenstates")
            if last is not None and last[0].shape == v.shape and np.array_equal(last[0], v):
                return last[1]
            eig = self._eigenstates(v, None)
            object.__setattr__(self, "_last_eigenstates", (v.copy(), eig))
            return eig
        return self._eigenstates(v, level_shift_ev)

    def _eigenstates(self, v: np.ndarray, level_shift_ev: np.ndarray | None) -> WindowEigenstates:
        c_inv = None if not self.varies else self._c_inv(v)
        window = max(self.window, STATIC_WINDOW)
        return window_eigenstates(self.induced_charge(v), window, lambda states, n_g: self._grid_energies(states, n_g, c_inv), self.tunnel_coupling_ev(v), level_shift_ev)

    def _grid_energies(self, states: np.ndarray, n_g: np.ndarray, c_inv: np.ndarray | None = None) -> np.ndarray:
        """U(N) (..., K) of the states (..., K, 2) around `n_g` (..., 2), +inf for a negative occupation."""
        d = states - n_g[..., None, :]
        u = self._quadratic(d, c_inv)
        return np.where((states >= 0).all(axis=-1), u, np.inf)

    def _kt(self, kt_ev: float | np.ndarray | None) -> float | np.ndarray:
        if kt_ev is None:
            return K_B_EV_PER_K * max(self.temperature_k, 1e-6)
        return np.maximum(np.asarray(kt_ev, dtype=float), 1e-12)[..., None]

    # --- exact derived quantities --------------------------------------------------------------

    def addition_energies_ev(self) -> np.ndarray:
        """Addition energy of each dot in eV: `E_UNIT · (C_dot⁻¹)_ii`."""
        return E_UNIT_V_PER_E_AF * np.diag(np.linalg.inv(self.c_dot_af))

    def mutual_charging_energy_ev(self) -> float:
        """Interdot electrostatic coupling in eV: `E_UNIT · (C_dot⁻¹)_12`; zero when Cm = 0."""
        return float(E_UNIT_V_PER_E_AF * np.linalg.inv(self.c_dot_af)[0, 1])

    def lever_arms_ev_per_v(self) -> np.ndarray:
        """Lever arm of each gate on each dot, [dot, gate], in eV per V: `C_dot⁻¹ · C_gate`."""
        return np.linalg.inv(self.c_dot_af) @ self.c_gate_af

    def transition_period_v(self, dot: int, gate: str) -> float:
        """Volts of `gate` per electron added to `dot` along that gate alone: `E_UNIT / C_gate[dot, gate]`."""
        j = self.gate_names.index(gate)
        return float(E_UNIT_V_PER_E_AF / self.c_gate_af[dot, j])

    def charge_lattice(self, x_gate: str = "P1", y_gate: str = "P2") -> np.ndarray:
        """Electrons per volt on (dot 1, dot 2) from (x_gate, y_gate): rows of C_gate / e, shape (2, 2).

        The honeycomb is periodic under one-electron shifts of the gate-induced charge, so these rows are
        its reciprocal lattice vectors in the (x, y) gate plane, away from the empty-dot corner.
        """
        i, j = self.gate_names.index(x_gate), self.gate_names.index(y_gate)
        return self.c_gate_af[:, [i, j]] / E_UNIT_V_PER_E_AF

    def triple_points_ng(self) -> tuple[np.ndarray, np.ndarray]:
        """The two triple points of the (0,0)/(1,0)/(0,1)/(1,1) vertex in gate-induced charge n_g (electrons).

        With A = C_dot⁻¹: the first solves A11 x + A12 y = A11/2 and A12 x + A22 y = A22/2 (states (0,0),
        (1,0), (0,1) degenerate); the second is the first plus A12 · (C1 − Cm, C2 − Cm) (states (1,0),
        (0,1), (1,1) degenerate). Their separation is the interdot transition and vanishes with Cm.
        """
        a = np.linalg.inv(self.c_dot_af)
        tp1 = np.linalg.solve(a, np.array([a[0, 0] / 2, a[1, 1] / 2]))
        tp2 = tp1 + a[0, 1] * (self.c_dot_af @ np.ones(2))
        return tp1, tp2

    def gate_voltages_for_ng(self, n_g: np.ndarray, x_gate: str = "P1", y_gate: str = "P2") -> np.ndarray:
        """Voltages of (x_gate, y_gate), other gates at 0 V, that give gate-induced charge `n_g` (2,)."""
        lattice = self.charge_lattice(x_gate, y_gate)
        return np.linalg.solve(lattice, np.asarray(n_g, dtype=float))
