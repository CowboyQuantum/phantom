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

Cost: one k×k matrix exponential per sample (k ≈ 25 for the default window), tens of microseconds.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm

from phantom.charge.ci import ConstantInteractionParams
from phantom.constants import K_B_EV_PER_K
from phantom.readout.rates import GateRates

MODEL_VERSION = 1

PROBABILITY_FLOOR = 1e-12
"""States below this probability are dropped from the active set between samples."""


def fermi(x: np.ndarray) -> np.ndarray:
    """1 / (1 + eˣ) with overflow guarded."""
    x = np.clip(np.asarray(x, dtype=float), -500, 500)
    return 1.0 / (1.0 + np.exp(x))


def _unique_states(states: np.ndarray) -> np.ndarray:
    """Distinct integer occupations, sorted, as (k, 2)."""
    return np.unique(np.asarray(states, dtype=int), axis=0)


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
        """The CI occupation window at `v`, without duplicates."""
        states, _ = self.charge.energies(self.charge.induced_charge(np.asarray(v, dtype=float)))
        return _unique_states(states)

    def rate_matrix(self, v: np.ndarray, states: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(energies (k,), W (k, k)) over the given states at `v`; W[j, i] is the rate from i to j and
        W[i, i] = −Σ_j W[j, i]."""
        v = np.asarray(v, dtype=float)
        states = np.asarray(states, dtype=int)
        u = self.charge.state_energies_ev(v, states)
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

    def states_and_rates(self, v: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(states, energies, W) over the CI window at `v`."""
        states = self.window_states(v)
        u, w = self.rate_matrix(v, states)
        return states, u, w

    def stationary(self, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """States of the window and the stationary distribution (null vector of W, normalised)."""
        states, _, w = self.states_and_rates(v)
        vals, vecs = np.linalg.eig(w)
        p = np.abs(np.real(vecs[:, np.argmin(np.abs(vals))]))
        return states, p / p.sum()

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
        probability, so nothing is lost or displaced as the window moves."""
        v_sequence = np.asarray(v_sequence, dtype=float)
        n = len(v_sequence)
        if p0 is None:
            states, p = self.stationary(v_sequence[0])
        else:
            states, p = np.asarray(states0, dtype=int), np.asarray(p0, dtype=float)
        mean = np.empty((n, 2))
        for t in range(n):
            keep = states[p > PROBABILITY_FLOOR]
            active = _unique_states(np.concatenate([self.window_states(v_sequence[t]), keep], axis=0))
            q = np.zeros(len(active))
            index = {tuple(s): i for i, s in enumerate(active.tolist())}
            for s, prob in zip(states.tolist(), p, strict=True):
                if prob > PROBABILITY_FLOOR:
                    q[index[tuple(s)]] += prob
            q /= q.sum()
            _, w = self.rate_matrix(v_sequence[t], active)
            q = expm(w * dwell_s) @ q
            q = np.maximum(q, 0.0)
            q /= q.sum()
            mean[t] = q @ active
            states, p = active, q
        return mean, p, states
