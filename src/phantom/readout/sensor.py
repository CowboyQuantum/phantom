"""A charge sensor: a single-electron transistor (SET) read in sequential tunnelling.

Seeded from qdt-sim's `DoubleDotModel` sensor (Step 1). The sensor is itself a single quantum dot
with charging energy `E_c = charging_energy_ev` and lever arm `α_s = E_c / period_v` (eV per V of its
potential φ_s), biased by a source–drain voltage `V_sd`. Its levels are `μ_k = E_c·k − α_s·(φ_s − phase_v)`;
sequential tunnelling at the sensor's thermal energy `kT_s = α_s·width_v / 2` gives

    I = I_bg + I_0 · Σ_k [ f((μ_k − eV_sd/2)/kT_s) − f((μ_k + eV_sd/2)/kT_s) ],   f(x) = 1/(1+eˣ),

a periodic train of thermally broadened Coulomb peaks at small bias and Coulomb diamonds of height
`E_c/e` against `V_sd`. `I_0 = 4 kT_s · peak_a / (e · bias_v)` makes `peak_a` the linear-response peak
current at the default bias. Because the sensor *is* a single dot in transport, the same functions
serve the single-dot protocols (`coulomb-oscillations-single-dot`, `diamonds-single-dot`).

The sensor's potential is linear in the gate voltages and in the charge on the dots it watches:

    φ_s = α · V + β · ⟨N⟩ + δφ                                [V-equivalent]

so every dot charge transition shows as a step in sensor current.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

MODEL_VERSION = 1
"""Bumped whenever any output of this module changes for the same parameters (ADR 0003)."""


@dataclass(frozen=True)
class SensorSETParams:
    """Parameters of the SET charge sensor and its coupling to the dots and gates.

    Defaults inherited from qdt-sim; a plausible Si-class device, not a measured one. Ranges: charging
    energy 0.3–3 meV, peak period 0.05–0.5 V, thermal width 5–50 mV, peak current 0.1–10 nA.
    """

    alpha_v_per_v: np.ndarray = field(default_factory=lambda: np.array([0.01, 0.01, 1.0]))
    """Sensor lever arm per gate: dφ_s/dV_gate, one entry per gate of the charge model."""
    beta_v_per_e: np.ndarray = field(default_factory=lambda: np.array([-0.004, -0.003]))
    """Sensor potential shift per electron on each dot (V-equivalent)."""
    period_v: float = 0.25
    """Spacing of the sensor's Coulomb peaks in φ_s (V)."""
    width_v: float = 0.04
    """Thermal width of a sensor peak in φ_s (V)."""
    peak_a: float = 1.0e-9
    """Linear-response peak current (A) at `bias_v`."""
    background_a: float = 0.05e-9
    """Sensor current far from a peak (A)."""
    phase_v: float = -0.06
    """φ_s of one sensor peak (V); the rest follow at `period_v`."""
    charging_energy_ev: float = 1.0e-3
    """Charging energy of the sensor SET (eV): e times the bias at a Coulomb diamond's apex."""
    bias_v: float = 50e-6
    """Source–drain bias the sensor is read at unless one is given (V)."""

    def __post_init__(self) -> None:
        if self.beta_v_per_e.shape != (2,):
            raise ValueError("beta_v_per_e must be (2,)")
        if self.period_v <= 0 or self.width_v <= 0 or self.charging_energy_ev <= 0 or self.bias_v <= 0:
            raise ValueError("period_v, width_v, charging_energy_ev and bias_v must be positive")

    @property
    def lever_arm_ev_per_v(self) -> float:
        """α_s: sensor level shift per volt of φ_s (eV/V); one period in φ_s adds one electron."""
        return self.charging_energy_ev / self.period_v

    @property
    def kt_ev(self) -> float:
        """Thermal energy of the sensor SET (eV), set so its zero-bias peaks have width `width_v`."""
        return self.lever_arm_ev_per_v * self.width_v / 2

    def potential(self, v: np.ndarray, n: np.ndarray, delta_phi_v: float | np.ndarray = 0.0) -> np.ndarray:
        """φ_s (V-equivalent) for gate voltages `v` (..., gates) and dot occupation `n` (..., 2)."""
        return (
            np.asarray(v, dtype=float) @ self.alpha_v_per_v
            + np.asarray(n, dtype=float) @ self.beta_v_per_e
            + delta_phi_v
        )

    def levels_ev(self, phi: np.ndarray) -> np.ndarray:
        """μ_k (..., 5) in eV: the five sensor levels nearest the bias window at potential `phi`."""
        x = (np.asarray(phi, dtype=float) - self.phase_v) / self.period_v
        k = np.rint(x)[..., None] + np.arange(-2, 3)
        return self.charging_energy_ev * (k - x[..., None])

    def current_a(self, phi: np.ndarray, v_sd: float | np.ndarray | None = None) -> np.ndarray:
        """Sensor current (A) at potential `phi` and source–drain bias `v_sd` (V; default `bias_v`)."""
        v_sd = self.bias_v if v_sd is None else np.asarray(v_sd, dtype=float)
        kt = self.kt_ev
        mu = self.levels_ev(phi)
        half = np.asarray(v_sd)[..., None] / 2  # eV numerically, for V_sd in V
        with np.errstate(over="ignore"):
            window = 1 / (1 + np.exp((mu - half) / kt)) - 1 / (1 + np.exp((mu + half) / kt))
        i_0 = 4 * kt * self.peak_a / self.bias_v
        return self.background_a + i_0 * window.sum(axis=-1)

    def diamond_edges(self, k: int) -> tuple[tuple[float, float], tuple[float, float]]:
        """The two Coulomb-diamond edge lines through the k-th zero-bias degeneracy point, as
        (φ_s intercept at V_sd = 0, slope dV_sd/dφ_s): the rising edge and the falling one.
        An edge is where level k crosses a lead: μ_k = ±eV_sd/2 → V_sd = ±2 α_s (φ_s − φ_k)."""
        phi_k = self.phase_v + k * self.period_v
        s = 2 * self.lever_arm_ev_per_v
        return (phi_k, s), (phi_k, -s)

    def peak_positions_v(self, phi_min: float, phi_max: float) -> np.ndarray:
        """φ_s of every zero-bias Coulomb peak in [phi_min, phi_max]."""
        k0 = np.ceil((phi_min - self.phase_v) / self.period_v)
        k1 = np.floor((phi_max - self.phase_v) / self.period_v)
        return self.phase_v + np.arange(k0, k1 + 1) * self.period_v
