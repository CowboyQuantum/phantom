"""Few-spin dynamics: Hamiltonians in the rotating frame, a Lindblad propagator, and the analytic
Rabi and Ramsey signals that are its noiseless limits (ADR 0021).

Energies are in eV and times in seconds; a frequency is E/h with h = 2πħ. One qubit is a spin-½ with
basis (|↑⟩, |↓⟩); n qubits are the Kronecker product in that order. In the frame rotating at the drive
frequency and with the rotating-wave approximation, spin i sees

    H_i = ½ δ_i σ_z + ½ Ω_i σ_x,      δ_i = g_i μ_B B − h f_drive,   Ω_i = h f_Rabi,

and a pair sees in addition the exchange H_J = J (S₁·S₂ − ¼), which is −J on the singlet and 0 on the
triplets. Decoherence enters as Lindblad operators: relaxation √(1/T₁) σ₋ and Markovian dephasing
√(1/(2 T_φ)) σ_z per spin. The slow part of the charge noise and the hyperfine field are not Markovian;
they are quasi-static offsets averaged over (`exchange.quasi_static_average`), which gives the
Gaussian envelope exp(−(τ/T₂*)²) with T₂* = √2 ħ / σ_E for an rms energy offset σ_E.

The propagator exponentiates the Liouvillian once per step with `scipy.linalg.expm` and applies it
repeatedly, which keeps the output reproducible across platforms (ADR 0010).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import expm

from phantom.constants import HBAR_EV_S
from phantom.valley.signatures import MU_B_UEV_PER_T

MODEL_VERSION = 1

H_EV_S = 2 * np.pi * HBAR_EV_S
"""Planck constant in eV·s."""
MU_B_EV_PER_T = MU_B_UEV_PER_T * 1e-6

SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
SM = np.array([[0, 0], [1, 0]], dtype=complex)
"""σ₋ = |↓⟩⟨↑|: relaxation from the upper state |↑⟩ to |↓⟩."""
I2 = np.eye(2, dtype=complex)


def embed(op: np.ndarray, i: int, n: int) -> np.ndarray:
    """`op` acting on qubit i of n."""
    out = np.array([[1.0 + 0j]])
    for k in range(n):
        out = np.kron(out, op if k == i else I2)
    return out


@dataclass(frozen=True)
class SpinParams:
    """Per-dot spin parameters for a pair. Placeholders (`calibration: none`) until the Step 6 reference
    entries are confirmed; the g-factor spread and dg/dV are the MOS/SiGe order of magnitude."""

    g_mean: float = 2.0
    g_sigma: float = 2e-3
    """Dot-to-dot spread of the g-factor (standard deviation)."""
    dg_dv_per_v: float = 1e-3
    """Gate dependence of g (per volt): turns gate-referred charge noise into frequency noise."""
    hyperfine_sigma_ev: float = 1e-9
    """rms Overhauser energy per dot: ~1 neV for natural silicon (T₂* ≈ 1 µs), ~0.05 neV enriched."""
    t1_s: float = 1.0
    """Spin relaxation time away from hot spots."""
    t_phi_markov_s: float = 100e-6
    """Markovian (fast, white) dephasing time."""

    def __post_init__(self) -> None:
        if self.g_sigma < 0 or self.hyperfine_sigma_ev < 0 or self.t1_s <= 0 or self.t_phi_markov_s <= 0:
            raise ValueError("spreads >= 0, times > 0")


@dataclass(frozen=True)
class SpinDraw:
    """One device's spin draw: g-factors per dot."""

    g_factors: tuple[float, ...]
    seed: int

    def zeeman_ev(self, b_t: float | np.ndarray) -> np.ndarray:
        """E_Z = g μ_B B per dot, (..., dots)."""
        return np.asarray(self.g_factors)[None, :] * MU_B_EV_PER_T * np.atleast_1d(np.asarray(b_t, dtype=float))[:, None]


def draw_spin(params: SpinParams, rng: np.random.Generator, seed: int = 0, n_dots: int = 2) -> SpinDraw:
    g = params.g_mean + params.g_sigma * rng.normal(0.0, 1.0, n_dots)
    return SpinDraw(tuple(float(x) for x in g), seed)


def quasi_static_sigma_ev(params: SpinParams, b_t: float, sigma_gate_v: float) -> float:
    """rms quasi-static energy offset of one spin: g-noise through dg/dV plus the hyperfine field."""
    return float(np.hypot(params.dg_dv_per_v * sigma_gate_v * MU_B_EV_PER_T * b_t, params.hyperfine_sigma_ev))


def t2_star_s(sigma_e_ev: float) -> float:
    """T₂* = √2 ħ / σ_E of the Gaussian envelope exp(−(τ/T₂*)²)."""
    return float(np.sqrt(2.0) * HBAR_EV_S / sigma_e_ev) if sigma_e_ev > 0 else float("inf")


# --- Hamiltonians --------------------------------------------------------------------------------


def rotating_frame_hamiltonian_ev(delta_ev: np.ndarray, omega_ev: np.ndarray, j_ev: float = 0.0) -> np.ndarray:
    """n-qubit H (2ⁿ × 2ⁿ) in eV: Σ ½(δ_i σ_z + Ω_i σ_x) and, for n = 2, J (S₁·S₂ − ¼)."""
    delta = np.atleast_1d(np.asarray(delta_ev, dtype=float))
    omega = np.atleast_1d(np.asarray(omega_ev, dtype=float))
    n = len(delta)
    h = np.zeros((2**n, 2**n), dtype=complex)
    for i in range(n):
        h += 0.5 * delta[i] * embed(SZ, i, n) + 0.5 * omega[i] * embed(SX, i, n)
    if n == 2 and j_ev != 0.0:
        s1s2 = sum(0.25 * embed(s, 0, 2) @ embed(s, 1, 2) for s in (SX, SY, SZ))
        h += j_ev * (s1s2 - 0.25 * np.eye(4))
    return h


def singlet_state() -> np.ndarray:
    """(|↑↓⟩ − |↓↑⟩)/√2 in the two-qubit basis."""
    psi = np.zeros(4, dtype=complex)
    psi[1], psi[2] = 1 / np.sqrt(2), -1 / np.sqrt(2)
    return psi


# --- Lindblad ------------------------------------------------------------------------------------


@dataclass
class Lindblad:
    """dρ/dt = −(i/ħ)[H, ρ] + Σ_k (C_k ρ C_k† − ½{C_k†C_k, ρ}) as a superoperator on the row-major
    vectorised ρ: vec(AρB) = (A ⊗ Bᵀ) vec(ρ)."""

    h_ev: np.ndarray
    collapse: list[np.ndarray] = field(default_factory=list)

    @classmethod
    def for_spins(cls, h_ev: np.ndarray, t1_s: float | np.ndarray | None = None, t_phi_s: float | np.ndarray | None = None) -> Lindblad:
        n = h_ev.shape[0].bit_length() - 1  # 2ⁿ states
        ops: list[np.ndarray] = []
        for i in range(n):
            if t1_s is not None:
                t1 = float(np.broadcast_to(t1_s, (n,))[i])
                if np.isfinite(t1):
                    ops.append(np.sqrt(1.0 / t1) * embed(SM, i, n))
            if t_phi_s is not None:
                tp = float(np.broadcast_to(t_phi_s, (n,))[i])
                if np.isfinite(tp):
                    ops.append(np.sqrt(1.0 / (2.0 * tp)) * embed(SZ, i, n))
        return cls(h_ev, ops)

    def superoperator(self) -> np.ndarray:
        d = self.h_ev.shape[0]
        eye = np.eye(d, dtype=complex)
        lio = -1j / HBAR_EV_S * (np.kron(self.h_ev, eye) - np.kron(eye, self.h_ev.T))
        for c in self.collapse:
            cdc = c.conj().T @ c
            lio += np.kron(c, c.conj()) - 0.5 * np.kron(cdc, eye) - 0.5 * np.kron(eye, cdc.T)
        return lio

    def evolve(self, rho0: np.ndarray, dt_s: float, n_steps: int) -> np.ndarray:
        """ρ at times k·dt for k = 0..n_steps: (n_steps + 1, d, d)."""
        d = rho0.shape[0]
        u = expm(self.superoperator() * dt_s)
        out = np.empty((n_steps + 1, d, d), dtype=complex)
        r = np.asarray(rho0, dtype=complex).reshape(-1)
        out[0] = rho0
        for k in range(1, n_steps + 1):
            r = u @ r
            out[k] = r.reshape(d, d)
        return out


def probability(rho: np.ndarray, state: np.ndarray) -> np.ndarray:
    """⟨state|ρ|state⟩ for ρ (..., d, d), real."""
    return np.real(np.einsum("i,...ij,j->...", state.conj(), rho, state))


# --- Analytic limits -----------------------------------------------------------------------------


def rabi_probability(tau_s: np.ndarray, delta_ev: np.ndarray, omega_ev: float) -> np.ndarray:
    """P(flip) for a drive of strength Ω at detuning δ for time τ, starting from |↑⟩:
    Ω²/(Ω² + δ²) · sin²(√(Ω² + δ²) τ / 2ħ)."""
    tau = np.asarray(tau_s, dtype=float)
    delta = np.asarray(delta_ev, dtype=float)
    w = np.sqrt(omega_ev**2 + delta**2)
    return omega_ev**2 / w**2 * np.sin(w * tau / (2 * HBAR_EV_S)) ** 2


def ramsey_probability(tau_s: np.ndarray, delta_ev: float, sigma_e_ev: float = 0.0, t2_markov_s: float = float("inf")) -> np.ndarray:
    """P(|↓⟩) after π/2 – wait τ – π/2 at detuning δ: ½(1 + cos(δτ/ħ) e^{−(τ/T₂*)²} e^{−τ/T₂})."""
    tau = np.asarray(tau_s, dtype=float)
    env = np.exp(-((tau / t2_star_s(sigma_e_ev)) ** 2)) * np.exp(-tau / t2_markov_s)
    return 0.5 * (1 + np.cos(delta_ev * tau / HBAR_EV_S) * env)


def singlet_return_probability(tau_s: np.ndarray, j_ev: np.ndarray) -> np.ndarray:
    """A separated singlet evolving under J for τ, projected back on the singlet: ½(1 + cos(Jτ/ħ)) in
    the usual fingerprint regime, where the Zeeman difference between the dots is much larger than J
    and the dynamics is a rotation between S and T₀."""
    return 0.5 * (1 + np.cos(np.asarray(j_ev, dtype=float) * np.asarray(tau_s, dtype=float) / HBAR_EV_S))


# --- EDSR driving --------------------------------------------------------------------------------


def edsr_rabi_energy_ev(g: float, gradient_t_per_nm: float, displacement_nm: float) -> float:
    """Rabi energy ħΩ of electric-dipole spin resonance in a micromagnet gradient: the oscillating
    displacement x_ac of the electron in a field gradient dB/dx gives an oscillating transverse field
    ½ x_ac dB/dx, so ħΩ = ½ g μ_B (dB/dx) x_ac (`calibration: none`)."""
    return 0.5 * g * MU_B_EV_PER_T * gradient_t_per_nm * displacement_nm


def drive_displacement_nm(v_ac_v: float, lever_arm_ev_per_v: float, orbital_energy_ev: float, radius_nm: float) -> float:
    """Electron displacement under a gate drive of amplitude V_ac: the dipole estimate α V_ac / E_orb
    times the orbital radius ℓ (`calibration: none`)."""
    return lever_arm_ev_per_v * v_ac_v / orbital_energy_ev * radius_nm
