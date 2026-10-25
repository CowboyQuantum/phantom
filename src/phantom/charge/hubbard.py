"""Two-site extended Hubbard model of a few-electron double dot.

Spinful fermions on two sites with at most two electrons per site (Fock space of 16 states):

    H = Σ_i ε_i(V) n_i + Σ_i U_i n_i↑ n_i↓ + V_12 n_1 n_2 − t(V_B) Σ_σ (c†_1σ c_2σ + h.c.)

The on-site energies are linear in the gate voltages through the lever arms, the tunnel coupling
`t(V_B)` is a `Barrier`, and the thermal state at the electron temperature gives the mean occupation.
The lead chemical potential is zero, so the Fock space is compared across electron numbers exactly as
the constant-interaction (CI) model compares its occupation window.

In the limit t → 0 the eigenstates are the charge states and, with `HubbardParams.from_ci`, the
energies equal the CI energies for N_i ≤ 2:  U_i = E_UNIT·(C⁻¹)_ii, V_12 = E_UNIT·(C⁻¹)_12 and
ε_i(V) = ½ E_UNIT (C⁻¹)_ii − (C⁻¹ C_gate V)_i. That identity is a P0 test.

Runtime cost: one batched 16×16 Hermitian diagonalisation per gate point.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phantom.charge.ci import E_UNIT_V_PER_E_AF, K_B_EV_PER_K, ConstantInteractionParams
from phantom.charge.dot import Barrier

MODEL_VERSION = 1
"""Bumped whenever any output of this module changes for the same parameters (ADR 0003)."""

# Fock basis: bit 0 = dot 1 spin up, bit 1 = dot 1 spin down, bit 2 = dot 2 up, bit 3 = dot 2 down.
_BASIS = np.arange(16)
_OCC = (_BASIS[:, None] >> np.arange(4)) & 1  # (16, 4)
N1 = (_OCC[:, 0] + _OCC[:, 1]).astype(float)
N2 = (_OCC[:, 2] + _OCC[:, 3]).astype(float)
D1 = (_OCC[:, 0] * _OCC[:, 1]).astype(float)
D2 = (_OCC[:, 2] * _OCC[:, 3]).astype(float)


def _hopping_matrix() -> np.ndarray:
    """Σ_σ (c†_1σ c_2σ + h.c.) in the Fock basis with Jordan–Wigner signs for the mode order 0,1,2,3."""
    h = np.zeros((16, 16))
    for a, b in ((0, 2), (1, 3)):  # (dot 1 mode, dot 2 mode) for spin up and spin down
        for s in _BASIS:
            if _OCC[s, b] and not _OCC[s, a]:
                between = int(_OCC[s, a + 1 : b].sum())
                sign = -1.0 if between % 2 else 1.0
                target = s ^ (1 << a) ^ (1 << b)
                h[target, s] += sign
                h[s, target] += sign
    return h


HOP = _hopping_matrix()


@dataclass(frozen=True)
class HubbardParams:
    """Parameters of the two-site Hubbard model; energies in eV, voltages in V."""

    gate_names: tuple[str, ...]
    """Gate voltages the model takes, in this order; includes the barrier gate."""
    lever_arms_ev_per_v: np.ndarray
    """On-site energy shift per volt, [dot, gate], in eV/V (a positive entry lowers ε_i as V rises)."""
    epsilon_0_ev: np.ndarray
    """On-site energies at zero gate voltage (2,)."""
    onsite_u_ev: np.ndarray
    """On-site charging energies U_i (2,)."""
    mutual_v_ev: float
    """Interdot electrostatic coupling V_12."""
    barrier: Barrier
    """t(V_B) for the gate named `barrier_gate`."""
    barrier_gate: str
    temperature_k: float = 0.15

    def __post_init__(self) -> None:
        n_g = len(self.gate_names)
        if self.lever_arms_ev_per_v.shape != (2, n_g):
            raise ValueError(f"lever_arms_ev_per_v must be (2, {n_g})")
        if self.epsilon_0_ev.shape != (2,) or self.onsite_u_ev.shape != (2,):
            raise ValueError("epsilon_0_ev and onsite_u_ev must be (2,)")
        if self.barrier_gate not in self.gate_names:
            raise ValueError(f"barrier_gate {self.barrier_gate!r} not in {self.gate_names}")

    @classmethod
    def from_ci(
        cls,
        ci: ConstantInteractionParams,
        barrier: Barrier,
        barrier_gate: str = "B",
        barrier_lever_arms_ev_per_v: tuple[float, float] = (0.0, 0.0),
    ) -> HubbardParams:
        """The Hubbard model whose t → 0 limit is the CI model `ci`, with `barrier` on a gate appended to
        the CI gate list (or an existing CI gate if `barrier_gate` already names one)."""
        a = np.linalg.inv(ci.c_dot_af)
        lever = a @ ci.c_gate_af
        if barrier_gate in ci.gate_names:
            names = ci.gate_names
        else:
            names = ci.gate_names + (barrier_gate,)
            lever = np.concatenate([lever, np.array(barrier_lever_arms_ev_per_v)[:, None]], axis=1)
        return cls(
            gate_names=names,
            lever_arms_ev_per_v=lever,
            epsilon_0_ev=0.5 * E_UNIT_V_PER_E_AF * np.diag(a),
            onsite_u_ev=E_UNIT_V_PER_E_AF * np.diag(a),
            mutual_v_ev=float(E_UNIT_V_PER_E_AF * a[0, 1]),
            barrier=barrier,
            barrier_gate=barrier_gate,
            temperature_k=ci.temperature_k,
        )

    # --- Hamiltonian and observables ----------------------------------------------------------

    def hamiltonian(self, v: np.ndarray) -> np.ndarray:
        """H (..., 16, 16) in eV for gate voltages `v` (..., gates)."""
        v = np.asarray(v, dtype=float)
        eps = self.epsilon_0_ev - v @ self.lever_arms_ev_per_v.T  # (..., 2)
        t = self.barrier.tunnel_coupling_ev(v[..., self.gate_names.index(self.barrier_gate)])
        diag = (
            eps[..., 0, None] * N1
            + eps[..., 1, None] * N2
            + self.onsite_u_ev[0] * D1
            + self.onsite_u_ev[1] * D2
            + self.mutual_v_ev * N1 * N2
        )  # (..., 16)
        h = np.zeros(diag.shape + (16,))
        idx = np.arange(16)
        h[..., idx, idx] = diag
        return h - t[..., None, None] * HOP

    def spectrum(self, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Eigenvalues (..., 16) ascending and eigenvectors (..., 16, 16) as columns."""
        return np.linalg.eigh(self.hamiltonian(v))

    def mean_occupation(self, v: np.ndarray) -> np.ndarray:
        """Thermal-average occupation ⟨N⟩ (..., 2) in electrons at the model temperature."""
        e, psi = self.spectrum(v)
        kt = K_B_EV_PER_K * max(self.temperature_k, 1e-6)
        w = np.exp(-(e - e[..., :1]) / kt)
        p = w / w.sum(axis=-1, keepdims=True)  # (..., n)
        weight = np.abs(psi) ** 2  # (..., basis, n)
        n1 = np.einsum("...n,b,...bn->...", p, N1, weight)
        n2 = np.einsum("...n,b,...bn->...", p, N2, weight)
        return np.stack([n1, n2], axis=-1)

    def ground_state_occupation(self, v: np.ndarray) -> np.ndarray:
        """⟨N⟩ (..., 2) in the ground state; not integer near an interdot anticrossing."""
        _, psi = self.spectrum(v)
        weight = np.abs(psi[..., :, 0]) ** 2
        return np.stack([weight @ N1, weight @ N2], axis=-1)
