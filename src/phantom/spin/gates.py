"""Fidelity of an exchange gate on the pair with the device's noise (ADR 0022).

A SWAP-type gate switches J on for t = πħ/J at an operating detuning. Its errors here are the
quasi-static detuning noise of the draw, which shifts J through dJ/dε during the gate, the Zeeman
difference of the g draw (which the ideal target includes, so it is not counted as an error unless it
is left out of the target), and the Markovian T₁ and T_φ. The entanglement fidelity of the realised
channel E against the target unitary U is

    F_e = (1/d²) Σ_{ij} ⟨i| U† E(|i⟩⟨j|) U |j⟩,    F_avg = (d F_e + 1) / (d + 1),

evaluated with the Lindblad superoperator (a linear map on vectorised operators) and averaged over the
quasi-static offset by quadrature. d = 4 for the pair.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import expm

from phantom.charge.hubbard import HubbardParams
from phantom.constants import HBAR_EV_S
from phantom.spin.dynamics import Lindblad, SpinParams, rotating_frame_hamiltonian_ev
from phantom.spin.exchange import detuning_voltages, exchange_j_ev

MODEL_VERSION = 1


def entanglement_fidelity(superoperator_t: np.ndarray, u_target: np.ndarray) -> float:
    """F_e of the channel given by the propagated superoperator (d² × d², row-major vec) against U."""
    d = u_target.shape[0]
    total = 0.0 + 0.0j
    for i in range(d):
        for j in range(d):
            op = np.zeros((d, d), dtype=complex)
            op[i, j] = 1.0
            out = (superoperator_t @ op.reshape(-1)).reshape(d, d)
            total += (u_target.conj().T @ out @ u_target)[i, j]
    return float(np.real(total) / d**2)


def average_gate_fidelity(f_e: float, d: int = 4) -> float:
    return (d * f_e + 1.0) / (d + 1.0)


def swap_gate_fidelity(
    hub: HubbardParams,
    v_center: np.ndarray,
    eps_ev: float,
    sigma_eps_ev: float,
    params: SpinParams,
    delta_zeeman_ev: float = 0.0,
    n_nodes: int = 15,
) -> tuple[float, float]:
    """(F_avg, gate time) of the SWAP gate J t = πħ at detuning `eps_ev`, averaged over Gaussian detuning
    offsets of rms `sigma_eps_ev`, with the Markovian T₁ and T_φ of `params`."""
    j_nom = float(exchange_j_ev(hub, detuning_voltages(hub, v_center, np.array([eps_ev])))[0])
    t_gate = np.pi * HBAR_EV_S / j_nom
    delta = np.array([delta_zeeman_ev / 2, -delta_zeeman_ev / 2])
    u_target = expm(-1j * rotating_frame_hamiltonian_ev(delta, np.zeros(2), j_nom) * t_gate / HBAR_EV_S)
    x, w = np.polynomial.hermite.hermgauss(n_nodes)
    offsets = np.sqrt(2.0) * sigma_eps_ev * x if sigma_eps_ev > 0 else np.zeros(1)
    weights = w / np.sqrt(np.pi) if sigma_eps_ev > 0 else np.ones(1)
    j_k = exchange_j_ev(hub, detuning_voltages(hub, v_center, eps_ev + offsets))
    f_e = 0.0
    for jk, wk in zip(j_k, weights):
        lb = Lindblad.for_spins(rotating_frame_hamiltonian_ev(delta, np.zeros(2), float(jk)), t1_s=params.t1_s, t_phi_s=params.t_phi_markov_s)
        f_e += wk * entanglement_fidelity(expm(lb.superoperator() * t_gate), u_target)
    return average_gate_fidelity(f_e), float(t_gate)
