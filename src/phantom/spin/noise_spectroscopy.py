"""Charge-noise spectroscopy with the qubit as the spectrometer (ADR 0022).

A qubit whose frequency depends on a gate-referred voltage (through dJ/dε for an exchange qubit, or
through dg/dV for a single spin) dephases under the voltage noise S_V(f). With a dynamical-decoupling
sequence of N π pulses in a total time τ the coherence is

    W_N(τ) = exp(−χ_N(τ)),    χ_N(τ) = ½ ∫₀^∞ df S_ω(f) F_N(2π f τ) / (2π f)²,

where S_ω is the one-sided spectrum of the qubit's angular-frequency noise (∫ S_ω df = σ_ω²) and F_N
is the filter function of the sequence, F_N(x) = |Σ_{k=0}^{N} (−1)^k (e^{i x t_{k+1}} − e^{i x t_k})|²
with the pulses at t_k = (k − ½)/N (CPMG spacing, times as fractions of τ). N = 0 is a Ramsey
experiment: for noise slow on the scale of τ it gives χ = σ_ω² τ² / 2, the Gaussian envelope with
T₂* = √2 / σ_ω. N = 1 is a Hahn echo. The filter of a long CPMG sequence is a narrow band at
f_N = N / (2τ), so W_N read against f_N is the spectrum, the measurement a tune-up makes.

The spectrum is the noise draw's: the 1/f amplitude is defined as the rms over the draw's band
[f_min, f_max], S_V(f) = A²/f with A² = σ_V² / ln(f_max/f_min), and the same A extends above the
sensor's band to the qubit's (the join the roster asks for, `noise-psd` and this protocol sharing one
amplitude); the extension's upper cut-off is a parameter.
"""

from __future__ import annotations

import numpy as np

from phantom.constants import HBAR_EV_S

MODEL_VERSION = 1


def one_over_f_amplitude_sq(sigma: float, f_min_hz: float, f_max_hz: float) -> float:
    """A² of S(f) = A²/f whose rms over [f_min, f_max] is σ."""
    return sigma**2 / np.log(f_max_hz / f_min_hz)


def one_over_f_psd(f_hz: np.ndarray, amplitude_sq: float, f_low_hz: float, f_high_hz: float) -> np.ndarray:
    """One-sided S(f) = A²/f inside [f_low, f_high], zero outside."""
    f = np.asarray(f_hz, dtype=float)
    return np.where((f >= f_low_hz) & (f <= f_high_hz), amplitude_sq / np.maximum(f, 1e-300), 0.0)


def cpmg_filter(x: np.ndarray, n_pulses: int) -> np.ndarray:
    """F_N(x) with x = ωτ, for N ≥ 0 ideal π pulses at CPMG positions."""
    x = np.asarray(x, dtype=float)
    edges = np.concatenate([[0.0], (np.arange(1, n_pulses + 1) - 0.5) / max(n_pulses, 1), [1.0]]) if n_pulses > 0 else np.array([0.0, 1.0])
    total = np.zeros(x.shape, dtype=complex)
    for k in range(len(edges) - 1):
        total += (-1) ** k * (np.exp(1j * x * edges[k + 1]) - np.exp(1j * x * edges[k]))
    return np.abs(total) ** 2


def coherence(tau_s: float, n_pulses: int, f_hz: np.ndarray, s_omega: np.ndarray) -> float:
    """W_N(τ) = exp(−χ) for a spectrum sampled on the (log-spaced) grid `f_hz`; the integral is done in
    ln f so that decades carry equal weight."""
    f = np.asarray(f_hz, dtype=float)
    omega = 2 * np.pi * f
    integrand = s_omega * cpmg_filter(omega * tau_s, n_pulses) / omega**2
    chi = 0.5 * float(np.trapezoid(integrand * f, np.log(f)))
    return float(np.exp(-chi))


def omega_spectrum(s_v: np.ndarray, sensitivity_ev_per_v: float) -> np.ndarray:
    """S_ω from a voltage spectrum: (dE/dV / ħ)² S_V."""
    return (sensitivity_ev_per_v / HBAR_EV_S) ** 2 * np.asarray(s_v)


def filter_frequency_hz(tau_s: float, n_pulses: int) -> float:
    return n_pulses / (2 * tau_s)


def log_grid_hz(f_low_hz: float, f_high_hz: float, per_decade: int = 400) -> np.ndarray:
    n = int(np.ceil(np.log10(f_high_hz / f_low_hz) * per_decade)) + 1
    return np.logspace(np.log10(f_low_hz), np.log10(f_high_hz), n)
