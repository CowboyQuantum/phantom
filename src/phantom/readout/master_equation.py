"""A master equation over the few-electron charge states of the double dot.

The states are the constant-interaction occupation window (N₁, N₂) around the gate-induced charge,
extended during a sweep by any state that still carries probability, so a charge that cannot follow
the gates is never pushed around by the moving window. Three kinds of transition connect them, each
with a rate that respects detailed balance against the Boltzmann distribution of the CI energies at
the electron temperature (ADR 0017):

- dot 1 exchanges an electron with the left lead at rate Γ_L: loading N₁ → N₁+1 at Γ_L · f(ΔU/kT) and
  unloading at Γ_L · (1 − f(ΔU/kT)), f the Fermi function and ΔU = U(N₁+1, N₂) − U(N₁, N₂);
- dot 2 likewise with the right lead at Γ_R;
- an electron moves between the dots at Γ_m: forward and backward rates Γ_m · f(±ΔU/kT), so their
  ratio is exp(−ΔU/kT).

With equal lead chemical potentials the stationary distribution is exactly the Boltzmann one, which the
CI model already uses for its thermal occupation (P0 test). The rates are constants, or, with
`gate_rates` (a compiled device, #53), they follow the gate voltages of each rate matrix through their
compiled slopes and saturate at the attempt rate (`readout.rates.GateRates`): closing a barrier gate slows
its lead by orders of magnitude while the equilibrium, which the prefactors do not enter, stays exact.
Propagating the distribution along a sweep
with the matrix exponential per sample gives the charge actually present: when a rate is small compared
with the inverse dwell time the distribution lags, which is latching, and the lag depends on the sweep
direction. With sequential tunnelling only (Γ_m = 0) a cold double dot can stay in a metastable charge
state when both lead paths out of it are uphill by many kT; that is real physics and a P0 test. This
replaces the Step 1 placeholder in `readout.loading` for everything that needs dynamics.

With an interdot coupling t in the CI model (#58, ADR 0029) the master equation runs over the
eigenstates of the hybridized window (`charge.hybrid.eigenbasis`), a Pauli master equation in the
secular approximation:

- a lead adds an electron to dot i, taking eigenstate k (N electrons) to l (N + 1), at
  Γ · |⟨l|c†_i|k⟩|² · f((E_l − E_k)/kT), and removes one at Γ · |⟨l|c†_i|k⟩|² · (1 − f);
- within one total charge, eigenstate k relaxes to l at Γ_m · g_kl · f((E_l − E_k)/kT), with
  g_kl = Σ_b (|c_ik|² + |c_jk|²)(|c_il|² + |c_jl|²) over the bonds b = (i, j) between charge states one
  moved electron apart. This is the sum over bonds of |⟨l|X_b|k⟩|² + |⟨l|Z_b|k⟩|², the bond's transfer
  and its dipole. At t = 0 it is Γ_m between neighbouring charge states and zero otherwise, and for an
  isolated pair it is Γ_m at every detuning.

Detailed balance holds in the eigenenergies, so the stationary state is the Boltzmann distribution of
the eigenstates, the CI model's thermal state with the coupling (P0). When the gates move between
samples, the distribution over the old eigenstates is projected onto the new ones,
q_l = Σ_k p_k |⟨l|k⟩|²: coherences are dropped (a sudden change followed by dephasing). At t = 0 this is
the old placement of each charge state. Without a coupling, every output is that of version 1, bit for
bit.

Cost: one k×k matrix exponential per sample (k ≈ 25 for the default window), tens of microseconds; with
a coupling, a few small diagonalisations more.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm

from phantom.charge.ci import ConstantInteractionParams
from phantom.charge.hybrid import dynamics_window, eigenbasis
from phantom.constants import K_B_EV_PER_K
from phantom.readout.rates import GateRates

MODEL_VERSION = 3
"""2: with an interdot coupling, the rate equation runs over the hybridized eigenstates (#58). 3: its
window follows t (`dynamics_window`), and with a shape its energies use C(V) (#59)."""

PROBABILITY_FLOOR = 1e-12
"""States below this probability are dropped from the active set between samples."""


def fermi(x: np.ndarray) -> np.ndarray:
    """1 / (1 + eˣ) with overflow guarded."""
    x = np.clip(np.asarray(x, dtype=float), -500, 500)
    return 1.0 / (1.0 + np.exp(x))


def _unique_states(states: np.ndarray) -> np.ndarray:
    """Distinct integer occupations, sorted, as (k, 2)."""
    return np.unique(np.asarray(states, dtype=int), axis=0)


def _move_pairs(states: np.ndarray, move: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    """(sources, destinations) index arrays of the states one `move` apart."""
    index = {tuple(s): i for i, s in enumerate(np.asarray(states).tolist())}
    pairs = [(i, index[(n1 + move[0], n2 + move[1])]) for i, (n1, n2) in enumerate(np.asarray(states).tolist()) if (n1 + move[0], n2 + move[1]) in index]
    if not pairs:
        return np.zeros(0, dtype=int), np.zeros(0, dtype=int)
    src, dst = np.array(pairs, dtype=int).T
    return src, dst


def carried_states(states: np.ndarray, vectors: np.ndarray, p: np.ndarray) -> np.ndarray:
    """The charge states (k', 2) that still carry probability above `PROBABILITY_FLOOR`, given
    probabilities `p` (m,) of eigenstates `vectors` (k, m) over `states` (k, 2); eigenstates at or below
    the floor count as empty."""
    p = np.where(p > PROBABILITY_FLOOR, p, 0.0)
    return np.asarray(states)[(vectors * vectors) @ p > PROBABILITY_FLOOR]


def project(states: np.ndarray, vectors: np.ndarray, p: np.ndarray, new_states: np.ndarray, new_vectors: np.ndarray) -> np.ndarray:
    """q (k',) normalised: the distribution `p` (m,) over eigenstates `vectors` (k, m) on `states`,
    projected onto the eigenstates `new_vectors` (k', k') on `new_states`, q_l = Σ_k p_k |⟨l|k⟩|².
    Eigenstates at or below `PROBABILITY_FLOOR` and charge states missing from `new_states` are dropped."""
    index = {tuple(s): i for i, s in enumerate(np.asarray(new_states).tolist())}
    embedded = np.zeros((len(new_states), vectors.shape[1]))
    for i, s in enumerate(np.asarray(states).tolist()):
        j = index.get(tuple(s))
        if j is not None:
            embedded[j] = vectors[i]
    overlap = new_vectors.T @ embedded
    q = (overlap * overlap) @ np.where(p > PROBABILITY_FLOOR, p, 0.0)
    return q / q.sum()


@dataclass(frozen=True)
class RateEquation:
    """The double dot's charge states and their transition rates."""

    charge: ConstantInteractionParams
    gamma_left_hz: float
    """Dot 1 to the left lead."""
    gamma_right_hz: float
    """Dot 2 to the right lead."""
    gamma_interdot_hz: float = 0.0
    """Electron transfer between the dots (inelastic, phenomenological)."""
    gate_rates: GateRates | None = None
    """When given, the (left, right, interdot) rates at each gate point come from it (#53), and the
    three constants above are its rates at the reference point; None keeps them constant."""

    def rates_hz(self, v: np.ndarray) -> np.ndarray:
        """(Γ_L, Γ_R, Γ_m) (..., 3), Hz, at carried voltages `v` (..., gates)."""
        if self.gate_rates is not None:
            return self.gate_rates.rates_hz(v)
        v = np.asarray(v, dtype=float)
        return np.broadcast_to(np.array([self.gamma_left_hz, self.gamma_right_hz, self.gamma_interdot_hz]), v.shape[:-1] + (3,)).copy()

    @property
    def kt_ev(self) -> float:
        return K_B_EV_PER_K * max(self.charge.temperature_k, 1e-6)

    def window_states(self, v: np.ndarray) -> np.ndarray:
        """The CI occupation window at `v`, without duplicates; with a coupling, as wide as t there needs
        (`charge.hybrid.dynamics_window`) for the truncated chains to stay exact (ADR 0029)."""
        v = np.asarray(v, dtype=float)
        ci = self.charge
        if ci.coupling is not None:
            ratio = float(ci.tunnel_coupling_ev(v)) / float(np.min(ci.addition_energies_ev()))
            ci = dataclasses.replace(ci, window=max(ci.window, dynamics_window(ratio)))
        states, _ = ci.energies(ci.induced_charge(v), None if ci.shape is None else ci._c_inv(v))
        return _unique_states(states)

    def rate_matrix(self, v: np.ndarray, states: np.ndarray, level_shift_ev: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
        """(energies (k,), W (k, k)) over the given states at `v`; W[j, i] is the rate from i to j and
        W[i, i] = −Σ_j W[j, i]. `level_shift_ev` (2,) adds that energy per electron on each dot (the
        charge sensor's back-action, ADR 0028); the rates keep detailed balance in the shifted energies."""
        v = np.asarray(v, dtype=float)
        states = np.asarray(states, dtype=int)
        u = self.charge.state_energies_ev(v, states)
        if level_shift_ev is not None:
            u = u + states @ np.asarray(level_shift_ev, dtype=float)
        k = len(u)
        w = np.zeros((k, k))
        index = {tuple(s): i for i, s in enumerate(states.tolist())}
        kt = self.kt_ev
        if self.gate_rates is None:
            g_l, g_r, g_m = self.gamma_left_hz, self.gamma_right_hz, self.gamma_interdot_hz
        else:
            g_l, g_r, g_m = (float(x) for x in self.gate_rates.rates_hz(v))
        moves = (((1, 0), g_l), ((0, 1), g_r), ((1, -1), g_m))
        listed = states.tolist()
        for (dn1, dn2), gamma in moves:
            if gamma <= 0:
                continue
            pairs = [(i, index[(n1 + dn1, n2 + dn2)]) for i, (n1, n2) in enumerate(listed) if (n1 + dn1, n2 + dn2) in index]
            if not pairs:
                continue
            src, dst = np.array(pairs, dtype=int).T
            du = u[dst] - u[src]
            # each entry is written once (no move is the reverse of another), so this equals the
            # element-by-element loop bit for bit
            w[dst, src] += gamma * fermi(du / kt)  # i → j
            w[src, dst] += gamma * fermi(-du / kt)  # j → i
        w[np.diag_indices(k)] = -w.sum(axis=0)
        return u, w

    def eigen_rate_matrix(self, v: np.ndarray, states: np.ndarray, level_shift_ev: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(energies (k,), eigenvectors (k, k) as columns over `states`, W (k, k) over the eigenstates)
        at `v`, with the interdot coupling of the CI model there (#58). Without a coupling the
        eigenvectors are the identity and W is `rate_matrix`'s, bit for bit."""
        v = np.asarray(v, dtype=float)
        states = np.asarray(states, dtype=int)
        t = float(self.charge.tunnel_coupling_ev(v))
        if t == 0.0:
            u, w = self.rate_matrix(v, states, level_shift_ev)
            return u, np.eye(len(u)), w
        u = self.charge.state_energies_ev(v, states)
        if level_shift_ev is not None:
            u = u + states @ np.asarray(level_shift_ev, dtype=float)
        e, c = eigenbasis(states, u, t)
        kt = self.kt_ev
        de = (e[:, None] - e[None, :]) / kt  # [l, k]: E_l − E_k
        up, down = fermi(de), fermi(-de)
        g_l, g_r, g_m = (float(x) for x in self.rates_hz(v))
        w = np.zeros((len(e), len(e)))
        for move, gamma in (((1, 0), g_l), ((0, 1), g_r)):
            src, dst = _move_pairs(states, move)
            if gamma <= 0 or not len(src):
                continue
            m = c[dst].T @ c[src]  # ⟨l|c†|k⟩
            weight = gamma * m * m
            w += weight * up + (weight * down).T
        src, dst = _move_pairs(states, (1, -1))
        if g_m > 0 and len(src):
            bond = c[src] ** 2 + c[dst] ** 2  # (bonds, k)
            g = bond.T @ bond
            np.fill_diagonal(g, 0.0)
            w += g_m * g * up
        w[np.diag_indices(len(e))] = 0.0
        w[np.diag_indices(len(e))] = -w.sum(axis=0)
        return e, c, w

    def states_and_rates(self, v: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(states, energies, W) over the CI window at `v`."""
        states = self.window_states(v)
        u, w = self.rate_matrix(v, states)
        return states, u, w

    def stationary(self, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """States of the window and the stationary probability of each (null vector of W, normalised;
        with a coupling, the charge states' populations in the stationary state of the eigenstates)."""
        states = self.window_states(v)
        _, c, w = self.eigen_rate_matrix(v, states)
        vals, vecs = np.linalg.eig(w)
        p = np.abs(np.real(vecs[:, np.argmin(np.abs(vals))]))
        return states, (c * c) @ (p / p.sum())

    def sweep(
        self,
        v_sequence: np.ndarray,
        dwell_s: float,
        p0: np.ndarray | None = None,
        states0: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Propagate the distribution along `v_sequence` (n, gates), one dwell per sample, starting from
        the stationary distribution at the first point (or `p0` over `states0`). Returns
        (⟨N⟩ (n, 2), final distribution, final states): each sample's distribution is the one *after* its
        dwell. The active state set at each sample is the CI window there plus every state still carrying
        probability, so nothing is lost or displaced as the window moves. With a coupling the
        distribution lives on the eigenstates and is projected between samples (`project`); the
        returned distribution is then the final charge states' populations, and `p0` a distribution
        over the charge states `states0`."""
        v_sequence = np.asarray(v_sequence, dtype=float)
        n = len(v_sequence)
        if p0 is None:
            states, p = self.stationary(v_sequence[0])
        else:
            states, p = np.asarray(states0, dtype=int), np.asarray(p0, dtype=float)
        vectors = np.eye(len(p))
        mean = np.empty((n, 2))
        for t in range(n):
            active = _unique_states(np.concatenate([self.window_states(v_sequence[t]), carried_states(states, vectors, p)], axis=0))
            _, c, w = self.eigen_rate_matrix(v_sequence[t], active)
            q = project(states, vectors, p, active, c)
            q = expm(w * dwell_s) @ q
            q = np.maximum(q, 0.0)
            q /= q.sum()
            mean[t] = ((c * c) @ q) @ active
            states, vectors, p = active, c, q
        return mean, (vectors * vectors) @ p, states
