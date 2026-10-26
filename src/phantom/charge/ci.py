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
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np

MODEL_VERSION = 1
"""Bumped whenever any output of this module changes for the same parameters (ADR 0003)."""

from phantom.constants import E_UNIT_V_PER_E_AF, K_B_EV_PER_K


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
        return np.asarray(v, dtype=float) @ self.c_gate_af.T / E_UNIT_V_PER_E_AF

    def energies(self, n_g: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Occupations `N` (..., k, 2) in the window around `n_g` (..., 2) and their energies (..., k) in eV."""
        n_g = np.asarray(n_g, dtype=float)
        center = np.rint(n_g).astype(int)
        offsets = np.array(
            list(itertools.product(range(-self.window, self.window + 1), repeat=2)), dtype=int
        )
        states = np.maximum(center[..., None, :] + offsets, 0)
        d = states - n_g[..., None, :]
        c_inv = np.linalg.inv(self.c_dot_af)
        u = 0.5 * E_UNIT_V_PER_E_AF * np.einsum("...ki,ij,...kj->...k", d, c_inv, d)
        return states, u

    def occupation_probabilities(
        self, v: np.ndarray, kt_ev: float | np.ndarray | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """States (..., k, 2) and Boltzmann probabilities (..., k) at the model temperature, or at the
        energy width `kt_ev` (scalar or one per gate point) when given, e.g. a lifetime-broadened width."""
        states, u = self.energies(self.induced_charge(v))
        if kt_ev is None:
            kt = K_B_EV_PER_K * max(self.temperature_k, 1e-6)
        else:
            kt = np.maximum(np.asarray(kt_ev, dtype=float), 1e-12)[..., None]
        w = np.exp(-(u - u.min(axis=-1, keepdims=True)) / kt)
        return states, w / w.sum(axis=-1, keepdims=True)

    def mean_occupation(self, v: np.ndarray, kt_ev: float | np.ndarray | None = None) -> np.ndarray:
        """Thermal-average occupation ⟨N⟩ (..., 2) in electrons for gate voltages `v`."""
        states, p = self.occupation_probabilities(v, kt_ev)
        return np.einsum("...k,...ki->...i", p, states)

    def ground_state_occupation(self, v: np.ndarray) -> np.ndarray:
        """The lowest-energy integer occupation (..., 2) for gate voltages `v` (the T → 0 charge state)."""
        states, u = self.energies(self.induced_charge(v))
        idx = np.argmin(u, axis=-1)
        return np.take_along_axis(states, idx[..., None, None], axis=-2)[..., 0, :]

    def sample_occupation(self, v: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        """One occupation (2,) drawn from the Boltzmann distribution at gate voltages `v` (1-D)."""
        states, p = self.occupation_probabilities(np.asarray(v, dtype=float))
        return states[rng.choice(len(p), p=p)]

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
