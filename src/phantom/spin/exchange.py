"""Exchange J(ε, t_c) from the two-site Hubbard model, and the detuning axis in gate space (ADR 0021).

Two electrons on two dots form a spin singlet and three triplets. With the tunnel coupling t the (1,1)
singlet mixes with the doubly occupied singlets (0,2) and (2,0) and is pushed down; the triplets are
not (Pauli exclusion). The exchange energy is the splitting

    J(V) = E_T(1,1) − E_S(1,1) ≥ 0,

read directly from the two-electron sector of the Hubbard Hamiltonian: the lowest triplet is the
|↑↑⟩ state, which is diagonal, and the lowest singlet is the ground state of the six-state sector.
Far from the anticrossings the second-order result is J ≈ 2t²/Δ_R + 2t²/Δ_L with Δ_R = E(0,2) − E(1,1)
and Δ_L = E(2,0) − E(1,1) (the (1,1)S–(0,2)S matrix element is √2 t); that limit is a P0 test.

Detuning ε is the on-site energy difference ε = ε₁ − ε₂ measured from the symmetric point, so ε > 0
pushes the pair toward (0,2). `detuning_voltages` moves the two plungers by the inverse of their
lever-arm block to realise a chosen ε at fixed ε₁ + ε₂. `exchange_sensitivity_ev_per_ev` is dJ/dε, the
number that turns charge noise into exchange dephasing.
"""

from __future__ import annotations

import numpy as np

from phantom.charge.hubbard import N1, N2, HubbardParams

MODEL_VERSION = 1

_TWO = np.flatnonzero(N1 + N2 == 2)  # the six two-electron Fock states
_UPUP = 0b0101  # dot 1 up, dot 2 up: a pure triplet, diagonal in H


def exchange_j_ev(hub: HubbardParams, v: np.ndarray) -> np.ndarray:
    """J (..., ) in eV for gate voltages `v` (..., gates): E(|↑↑⟩) − E(lowest two-electron singlet)."""
    h = hub.hamiltonian(v)
    block = h[..., _TWO[:, None], _TWO[None, :]]
    e_s = np.linalg.eigvalsh(block)[..., 0]
    e_t = h[..., _UPUP, _UPUP]
    return e_t - e_s


def exchange_second_order_ev(hub: HubbardParams, v: np.ndarray) -> np.ndarray:
    """The perturbative J = 2t²/Δ_R + 2t²/Δ_L for comparison; valid where t ≪ Δ."""
    v = np.asarray(v, dtype=float)
    eps = hub.epsilon_0_ev - v @ hub.lever_arms_ev_per_v.T
    t = hub.barrier.tunnel_coupling_ev(v[..., hub.gate_names.index(hub.barrier_gate)])
    e11 = eps[..., 0] + eps[..., 1] + hub.mutual_v_ev
    d_r = 2 * eps[..., 1] + hub.onsite_u_ev[1] - e11
    d_l = 2 * eps[..., 0] + hub.onsite_u_ev[0] - e11
    return 2 * t**2 / d_r + 2 * t**2 / d_l


def symmetric_point_v(hub: HubbardParams, fixed: dict[str, float] | None = None, plungers: tuple[str, str] = ("P1", "P2")) -> np.ndarray:
    """Gate voltages (gates,) at the centre of the (1,1) cell: on-site energies ε_i = −U_i/2 − V₁₂, which
    put the pair midway between the one- and three-electron boundaries and between the (2,0) and
    (0,2) anticrossings. Gates other than the two plungers take the values in `fixed` (default 0 V)."""
    fixed = fixed or {}
    v = np.array([float(fixed.get(g, 0.0)) for g in hub.gate_names])
    ip = [hub.gate_names.index(g) for g in plungers]
    target = -hub.onsite_u_ev / 2 - hub.mutual_v_ev
    others = [i for i in range(len(hub.gate_names)) if i not in ip]
    rhs = hub.epsilon_0_ev - hub.lever_arms_ev_per_v[:, others] @ v[others] - target
    v[ip] = np.linalg.solve(hub.lever_arms_ev_per_v[:, ip], rhs)
    return v


def detuning_voltages(hub: HubbardParams, v_center: np.ndarray, eps_ev: np.ndarray, plungers: tuple[str, str] = ("P1", "P2")) -> np.ndarray:
    """Voltages (..., gates) that shift the on-site energies by (+ε/2, −ε/2) from `v_center` using the
    two plungers only (ε₁ + ε₂ fixed): ε > 0 lowers dot 2 and pushes the pair toward (0,2)."""
    ip = [hub.gate_names.index(g) for g in plungers]
    eps = np.asarray(eps_ev, dtype=float)
    d_eps = np.stack([eps / 2, -eps / 2], axis=-1)  # (..., 2)
    dv = -np.linalg.solve(hub.lever_arms_ev_per_v[:, ip], d_eps[..., None])[..., 0]
    vc = np.asarray(v_center, dtype=float)
    shape = np.broadcast_shapes(eps.shape, vc.shape[:-1]) + (len(hub.gate_names),)
    v = np.broadcast_to(vc, shape).copy()
    v[..., ip] += np.broadcast_to(dv, shape[:-1] + (2,))
    return v


def detuning_noise_sigma_ev(hub: HubbardParams, sigma_gate_v: float, gates: tuple[str, ...] = ("P1", "P2")) -> float:
    """rms detuning noise from independent rms voltage noise `sigma_gate_v` on each named gate:
    σ_ε² = Σ_j (α₁ⱼ − α₂ⱼ)² σ_V²."""
    ig = [hub.gate_names.index(g) for g in gates]
    diff = hub.lever_arms_ev_per_v[0, ig] - hub.lever_arms_ev_per_v[1, ig]
    return float(sigma_gate_v * np.sqrt(np.sum(diff**2)))


def exchange_sensitivity_ev_per_ev(hub: HubbardParams, v_center: np.ndarray, eps_ev: np.ndarray, plungers: tuple[str, str] = ("P1", "P2")) -> np.ndarray:
    """dJ/dε by the Hellmann–Feynman theorem. Along the detuning axis the Hamiltonian changes as
    dH/dε = ½(N₁ − N₂) (the plungers move so that ε₁ gains +ε/2 and ε₂ gains −ε/2; the barrier gate
    does not move). The triplet |↑↑⟩ has N₁ = N₂, so dJ/dε = −⟨S|½(N₁ − N₂)|S⟩ = ½⟨N₂ − N₁⟩_S: half the
    charge polarisation of the singlet ground state. Exact, and free of the digit cancellation of a
    finite difference of eigenvalues (which differed between platforms beyond ADR 0014's tolerance)."""
    v = detuning_voltages(hub, v_center, eps_ev, plungers)
    block = hub.hamiltonian(v)[..., _TWO[:, None], _TWO[None, :]]
    _, psi = np.linalg.eigh(block)
    weight = np.abs(psi[..., :, 0]) ** 2  # (..., 6) over the two-electron Fock states
    polarisation = weight @ (N2[_TWO] - N1[_TWO])
    return 0.5 * polarisation


def quasi_static_average(f, sigma: float, n_nodes: int = 31) -> np.ndarray:
    """⟨f(δ)⟩ over δ ~ N(0, σ²) by Gauss–Hermite quadrature: Σ w_k f(√2 σ x_k)/√π. `f` takes an array of
    offsets and returns an array whose last axis is the offset axis."""
    x, w = np.polynomial.hermite.hermgauss(n_nodes)
    vals = f(np.sqrt(2.0) * sigma * x)
    return np.tensordot(vals, w, axes=([-1], [0])) / np.sqrt(np.pi)
