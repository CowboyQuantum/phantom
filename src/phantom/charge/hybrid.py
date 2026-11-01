"""Interdot tunnelling in the constant-interaction window: the charge states hybridize (#58, ADR 0029).

The CI model's charge states (N₁, N₂) are classical. An interdot tunnel coupling t mixes the states of
one total charge N = N₁ + N₂ that differ by one electron moved between the dots:

    H_N = Σ U(N₁, N₂) |N₁ N₂⟩⟨N₁ N₂| − t Σ (|N₁+1, N₂−1⟩⟨N₁ N₂| + h.c.),

a chain over the window's states of that charge, with U the CI energies. This is spinless and has one
level per dot: the polarization of a charge qubit, with its anticrossing gap 2t. Its eigenstates take
the place of the charge states. The thermal state is the Boltzmann distribution over their energies,
and the charge the sensor reads is the expectation Σ_i |c_ik|² N_i of each eigenstate. Across an
interdot line at detuning ε, the polarization of the two states then follows
½(1 − ε/Ω · tanh(Ω/2kT)), Ω = √(ε² + 4t²), which is wider than the thermal 2kT once t ≳ kT. At t = 0 the
eigenstates are the charge states and everything is the CI model (P0).

`InterdotCoupling` gives t at the carried gate voltages. With an `InterdotBarrier` (#70, ADR 0034) it
is t = t_max·√T, T the barrier's one transmission (Kemble's, of its screened top continued past the
merge), shared with the interdot rate and the conductance. Without one the exponent moves from its
compiled value with half the slopes of ln Γ_m, since the interdot rate and t² come from the same
transmission (Γ_m ∝ T, t ∝ √T, ADR 0027), and saturates at t_max. t_max = ħω₀/π is compile's estimate
t = ħω₀√T/π at full transmission (the barrier between the dots gone).

The window is a (2w+1)² grid around the gate-induced charge, its centre held at or above zero
electrons. A state with a negative occupation is left out: its energy is +inf and it couples to nothing.
Each total charge is one block, a chain of at most 2w + 1 states, so a gate point costs a few small
symmetric diagonalisations, done in batches. Truncating a chain leaves out admixtures of order
(t/E_add)^(2w); the centre rint(n_g) changes at half-integer n_g, the middle of every interdot segment, so
the truncation must stay far below the cross-platform tolerance there (ADR 0014). The statics use
`STATIC_WINDOW` = 4 (below 10⁻¹² up to t's ceiling), and the dynamics `dynamics_window(t / E_add)`, the
smallest window whose truncation is below 10⁻¹³.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import expit

from phantom.charge.barrier import InterdotBarrier
from phantom.charge.shape import ShapeGrid

MODEL_VERSION = 4
"""2: the statics' window is ±4 and the dynamics' follows t (#59, ADR 0029): the truncation no longer
shows at the middle of an interdot segment. 3: with an `InterdotBarrier`, t = t_max·√T of the barrier's
one transmission (#70, ADR 0034). 4: the puddles' occupancy scale reaches the barrier (#72, ADR 0035)."""

STATIC_WINDOW = 4
"""Half-width of the hybridized window in the statics: chains of up to nine states."""


def dynamics_window(t_over_e_add: float) -> int:
    """The smallest window half-width whose truncation, about (t/E_add)^(2w), stays below 10⁻¹³ at the
    ratio of t to the smaller addition energy: 2 below 10⁻³, 3 below 1.5·10⁻², 4 above (measured on the
    reference device, where the error at t/E_add = 0.1 is 2·10⁻⁵, 4·10⁻⁹ and 2·10⁻¹³ for w = 2, 3, 4)."""
    if t_over_e_add < 1e-3:
        return 2
    return 3 if t_over_e_add < 1.5e-2 else 4


EXCLUDED_GAP_EV = 1.0
"""How far above the block's highest allowed level a left-out state is put before the diagonalisation,
eV. It couples to nothing, so the value only has to stay clear of the allowed levels."""


@dataclass(frozen=True)
class InterdotCoupling:
    """The interdot tunnel coupling t at carried gate voltages: compiled at a reference point and moved
    with the slopes of its exponent, saturating at `t_max_ev`."""

    v_ref_v: tuple[float, ...]
    """Carried gate voltages t is compiled at, one per gate, V."""
    t_ref_ev: float
    """t at `v_ref_v`, eV."""
    t_max_ev: float
    """The ceiling of t, eV: ħω₀/π, the coupling through a transmission of one."""
    log_slopes_per_v: tuple[float, ...]
    """Slopes of the exponent of t against every gate, 1/V: half those of ln Γ_m."""
    shape: ShapeGrid | None = None
    """When given, the exponent along its barrier gates is half the change of the interdot action on
    its grid (#59), and the slopes count for the other gates only."""
    barrier: InterdotBarrier | None = None
    """When given (#70), t = t_max·√T with T the barrier's transmission at the gates, and the slopes and
    the shape grid's action are not used; `t_ref_ev` is then t_max·√T at the reference point."""

    def __post_init__(self) -> None:
        if not (self.t_ref_ev > 0 and self.t_max_ev > 0):
            raise ValueError("t_ref_ev and t_max_ev must be positive")
        if len(self.log_slopes_per_v) != len(self.v_ref_v):
            raise ValueError("slopes need one entry per gate")

    def tunnel_coupling_ev(self, v: np.ndarray, scale: np.ndarray | None = None) -> np.ndarray:
        """t (...) in eV at carried voltages `v` (..., gates). With a `barrier`, t_max·√T(V), the barrier
        seen from the puddles' edges when their `scale` (..., 2) is given (#72). Otherwise
        t_max · σ(ln(t_ref/(t_max − t_ref)) + s·ΔV), exactly t_ref at the reference point (the form of
        `readout.rates.saturating_rates_hz`)."""
        v = np.asarray(v, dtype=float)
        if self.barrier is not None:
            return self.t_max_ev * np.sqrt(self.barrier.transmission(v, scale))
        dv = v - np.asarray(self.v_ref_v)
        anchor = np.log(self.t_ref_ev) - np.log(max(self.t_max_ev - self.t_ref_ev, 1e-12 * self.t_max_ev))
        if self.shape is None:
            return self.t_max_ev * expit(anchor + dv @ np.asarray(self.log_slopes_per_v))
        slopes = np.array(self.log_slopes_per_v)
        slopes[list(self.shape.gate_indices)] = 0.0
        along = 0.5 * self.shape.quantity(self.shape.change(v), "log_rate")[..., 2]
        return self.t_max_ev * expit(anchor + dv @ slopes + along)


def grid_offsets(window: int) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
    """The window's (2w+1)² occupation offsets (K, 2), ordered block by block (total offset a + b
    ascending, then a ascending), and each block's indices into them. Consecutive states of a block
    differ by one electron moved from dot 2 to dot 1."""
    offsets, blocks, start = [], [], 0
    for s in range(-2 * window, 2 * window + 1):
        chain = [(a, s - a) for a in range(-window, window + 1) if -window <= s - a <= window]
        offsets.extend(chain)
        blocks.append(np.arange(start, start + len(chain)))
        start += len(chain)
    return np.array(offsets, dtype=int), tuple(blocks)


@dataclass(frozen=True)
class WindowEigenstates:
    """The hybridized window at a batch of gate points: K states per point."""

    center: np.ndarray
    """(..., 2) the window's centre, electrons: charges are summed relative to it, so their digits do
    not depend on the size of the occupation."""
    states: np.ndarray
    """(..., K, 2) charge states of the window, block by block (`grid_offsets`); entries with a
    negative occupation are left out (energy +inf)."""
    energies_ev: np.ndarray
    """(..., K) eigenstate energies, eV; +inf for the left-out states. Eigenstate k belongs to the
    block of window state k."""
    weights: tuple[np.ndarray, ...]
    """Per block, (..., m, m): |c_ik|², the weight of eigenstate k on charge state i of the block."""
    blocks: tuple[np.ndarray, ...]
    """Per block, its indices into the K states and eigenstates."""

    def charge(self) -> np.ndarray:
        """(..., K, 2) the charge expectation Σ_i |c_ik|² N_i of each eigenstate, electrons."""
        out = np.zeros(self.states.shape, dtype=float)
        rel = (self.states - self.center[..., None, :]).astype(float)
        for idx, w in zip(self.blocks, self.weights, strict=True):
            out[..., idx, :] = np.einsum("...ik,...ix->...kx", w, rel[..., idx, :])
        return self.center[..., None, :] + out

    def mean_charge(self, p: np.ndarray) -> np.ndarray:
        """(..., 2) ⟨N⟩ given probabilities `p` (..., K) of the eigenstates, summed relative to the
        centre."""
        rel = (self.states - self.center[..., None, :]).astype(float)
        return self.center + np.einsum("...k,...kx->...x", self.populations(p), rel)

    def populations(self, p: np.ndarray) -> np.ndarray:
        """(..., K) the probability Σ_k p_k |c_ik|² of each charge state, given probabilities `p`
        (..., K) of the eigenstates: the diagonal of the density matrix in the charge basis."""
        out = np.zeros(p.shape, dtype=float)
        for idx, w in zip(self.blocks, self.weights, strict=True):
            out[..., idx] = np.einsum("...ik,...k->...i", w, p[..., idx])
        return out

    def dominant_state(self, k: np.ndarray) -> np.ndarray:
        """(..., 2) the charge state with the largest weight in eigenstate `k` (...)."""
        k = np.asarray(k)
        best = np.zeros(k.shape, dtype=int)
        for idx, w in zip(self.blocks, self.weights, strict=True):
            inside = (k >= idx[0]) & (k <= idx[-1])
            col = np.clip(k - idx[0], 0, len(idx) - 1)
            column = np.take_along_axis(w, col[..., None, None], axis=-1)[..., 0]  # (..., m)
            best = np.where(inside, idx[0] + np.argmax(column, axis=-1), best)
        return np.take_along_axis(self.states, best[..., None, None], axis=-2)[..., 0, :]


def window_eigenstates(n_g: np.ndarray, window: int, energy_fn, t_ev: np.ndarray, level_shift_ev: np.ndarray | None = None) -> WindowEigenstates:
    """Diagonalise the window around gate-induced charge `n_g` (..., 2). `energy_fn(states, n_g)` gives
    the CI energies (..., K) of the states (..., K, 2), `t_ev` (...) the coupling at each point, and
    `level_shift_ev` (2,), when given, adds that energy per electron on each dot."""
    n_g = np.asarray(n_g, dtype=float)
    offsets, blocks = grid_offsets(window)
    center = np.maximum(np.rint(n_g).astype(int), 0)
    states = center[..., None, :] + offsets
    allowed = (states >= 0).all(axis=-1)
    u = energy_fn(states, n_g)
    if level_shift_ev is not None:
        u = u + states @ np.asarray(level_shift_ev, dtype=float)
    t = np.asarray(t_ev, dtype=float)[..., None]
    energies = np.empty(u.shape)
    weights = []
    for idx in blocks:
        m = len(idx)
        ub, ok = u[..., idx], allowed[..., idx]
        # the block's lowest allowed level is taken out before the diagonalisation, so the matrix's norm
        # is the spread of its levels, not their size, and the eigenvalues keep their digits (ADR 0014)
        low = np.min(np.where(ok, ub, np.inf), axis=-1, keepdims=True)
        low = np.where(np.isfinite(low), low, 0.0)
        top = np.max(np.where(ok, ub - low, -np.inf), axis=-1, keepdims=True)
        top = np.where(np.isfinite(top), top, 0.0)
        h = np.zeros(u.shape[:-1] + (m, m))
        h[..., np.arange(m), np.arange(m)] = np.where(ok, ub - low, top + EXCLUDED_GAP_EV)
        if m > 1:
            hop = np.where(ok[..., :-1] & ok[..., 1:], -t, 0.0)
            h[..., np.arange(m - 1), np.arange(1, m)] = hop
            h[..., np.arange(1, m), np.arange(m - 1)] = hop
        e, c = np.linalg.eigh(h)
        w = c * c
        # a left-out state couples to nothing, so its eigenvector is its own unit vector
        out = np.einsum("...ik,...i->...k", w, (~ok).astype(float)) > 0.5
        energies[..., idx] = np.where(out, np.inf, e + low)
        weights.append(w)
    return WindowEigenstates(center.astype(float), states, energies, tuple(weights), blocks)


def boltzmann(energies_ev: np.ndarray, kt_ev: float | np.ndarray) -> np.ndarray:
    """Normalised Boltzmann probabilities over the last axis of `energies_ev`, +inf entries at zero."""
    e = np.asarray(energies_ev, dtype=float)
    w = np.exp(-(e - np.min(e, axis=-1, keepdims=True)) / kt_ev)
    return w / w.sum(axis=-1, keepdims=True)


def eigenbasis(states: np.ndarray, u: np.ndarray, t_ev: float) -> tuple[np.ndarray, np.ndarray]:
    """(energies (k,), eigenvectors (k, k) as columns) of the window Hamiltonian over arbitrary charge
    states `states` (k, 2) with energies `u` (k,): one chain per total charge, coupled by −t between
    states one moved electron apart. Columns that belong to one total charge are nonzero only on it.
    With t = 0 the result is (u, identity) exactly."""
    states = np.asarray(states, dtype=int)
    u = np.asarray(u, dtype=float)
    k = len(u)
    if t_ev == 0.0:
        return u.copy(), np.eye(k)
    energies, vectors = np.empty(k), np.zeros((k, k))
    index = {tuple(s): i for i, s in enumerate(states.tolist())}
    total = states.sum(axis=1)
    for n in np.unique(total):
        members = np.flatnonzero(total == n)
        order = members[np.argsort(states[members, 0], kind="stable")]
        m = len(order)
        low = float(np.min(u[order]))
        h = np.diag(u[order] - low)  # the norm is the spread of the levels (ADR 0014)
        for a in range(m):
            n1, n2 = states[order[a]]
            j = index.get((int(n1) + 1, int(n2) - 1))
            if j is not None:
                b = int(np.flatnonzero(order == j)[0])
                h[a, b] = h[b, a] = -t_ev
        e, c = np.linalg.eigh(h)
        energies[order] = e + low
        vectors[np.ix_(order, order)] = c
    return energies, vectors


__all__ = ["EXCLUDED_GAP_EV", "MODEL_VERSION", "STATIC_WINDOW", "InterdotCoupling", "WindowEigenstates", "boltzmann", "dynamics_window", "eigenbasis", "grid_offsets", "window_eigenstates"]
