"""What a valley splitting does to measurements that do not need the spin layer (ADR 0020).

- **Magnetospectroscopy of the second electron.** With one electron in the lower valley, adding a
  second makes a spin singlet in the same valley–orbital state while the Zeeman energy is below the
  valley splitting, and a spin triplet with one electron in the upper valley once it is above. The
  addition energy of the second electron therefore rises with field as +½ g μ_B B until g μ_B B = E_V
  and falls as −½ g μ_B B after: the transition position in gate voltage has a kink at B* = E_V / (g μ_B).
- **Relaxation hot spot.** Where the Zeeman energy crosses the valley splitting the spin and valley
  states mix and spin relaxation speeds up by orders of magnitude; T₁(B) has a dip at B*. Modelled as a
  Lorentzian in E_Z − E_V on top of a background rate; the background's own field dependence waits
  for the spin layer.
- **Pulsed-gate spectroscopy.** A square pulse of amplitude A on the plunger during a charge-sensing
  sweep splits a transition line in two, A apart; loading through the valley-excited state adds a
  weaker line displaced by E_V / α from the pulsed line (an idealised picture with fixed weights,
  `calibration: none`).
"""

from __future__ import annotations

import numpy as np

MU_B_UEV_PER_T = 57.883818060
"""Bohr magneton in µeV/T (CODATA 2018)."""


def second_electron_shift_uev(e_valley_uev: float, b_t: np.ndarray, g_factor: float = 2.0) -> np.ndarray:
    """Shift of the second electron's addition energy with field, µeV, relative to zero field."""
    ez = g_factor * MU_B_UEV_PER_T * np.asarray(b_t, dtype=float)
    return np.where(ez < e_valley_uev, 0.5 * ez, e_valley_uev - 0.5 * ez)


def kink_field_t(e_valley_uev: float, g_factor: float = 2.0) -> float:
    return float(e_valley_uev / (g_factor * MU_B_UEV_PER_T))


def relaxation_rate_hz(b_t: np.ndarray, e_valley_uev: float, background_hz: float = 1.0, peak_hz: float = 1e5, width_uev: float = 5.0, g_factor: float = 2.0) -> np.ndarray:
    """Spin relaxation rate 1/T₁ against field: a Lorentzian hot spot where E_Z = E_V."""
    ez = g_factor * MU_B_UEV_PER_T * np.asarray(b_t, dtype=float)
    return background_hz + peak_hz * width_uev**2 / ((ez - e_valley_uev) ** 2 + width_uev**2)


def pulsed_spectroscopy_occupation(
    v_plunger_v: np.ndarray,
    pulse_amplitude_v: np.ndarray,
    v_transition_v: float,
    e_valley_uev: float,
    lever_arm_ev_per_v: float,
    kt_uev: float = 13.0,
    excited_weight: float = 0.4,
) -> np.ndarray:
    """Time-averaged occupation of the second electron under a 50 % duty square pulse: half the time the
    plunger sits at V, half at V + A. Each half loads through the ground state (a thermally broadened
    step at the transition) and, in the pulsed half, additionally through the valley-excited state at
    E_V / α lower plunger voltage with weight `excited_weight`. The excited-state loading is a
    non-equilibrium effect of the fast pulse edge, so its weight is scaled by 1 − exp(−A / w) with w
    the thermal width in volts: it vanishes smoothly as the amplitude goes to zero, where the sweep
    is an equilibrium one. Broadcasts V (n,) against A (m,) → (m, n)."""
    v = np.asarray(v_plunger_v, dtype=float)[None, :]
    a = np.asarray(pulse_amplitude_v, dtype=float)[:, None]
    width_v = kt_uev * 1e-6 / lever_arm_ev_per_v
    shift_v = e_valley_uev * 1e-6 / lever_arm_ev_per_v

    def step(x):
        return 1.0 / (1.0 + np.exp(-(x) / width_v))

    low = step(v - v_transition_v)
    high_ground = step(v + a - v_transition_v)
    high_excited = step(v + a - v_transition_v + shift_v)
    weight = excited_weight * (1.0 - np.exp(-np.abs(a) / width_v))
    high = high_ground + weight * (high_excited - high_ground)
    return 0.5 * (low + high)
