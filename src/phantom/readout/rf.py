"""Radio-frequency reflectometry: a lumped tank circuit with the charge sensor as its load.

The sensor's source–drain resistance R sits in parallel with the parasitic capacitance C_p at the end
of a series inductor L driven through a line of impedance Z₀ (ADR 0018):

    Z(f) = i ω L + 1 / (i ω C_p + 1 / R),     Γ(f) = (Z − Z₀) / (Z + Z₀).

Resonance is at f₀ = 1 / (2π √(L C_p)); at resonance the reflected amplitude vanishes when the sensor
resistance equals the matching value R_match = L / (C_p Z₀), and the phase of Γ flips by π as R crosses
it. A change of sensor conductance therefore changes the reflected amplitude (and phase) at the carrier
frequency, which is what a reflectometry readout demodulates. Cables, amplifier noise and the varactor
of a tunable match are not modelled; the sensor resistance is the small-bias ratio V_sd / I of the SET
model, so the RF flank mirrors the DC Coulomb peak.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MODEL_VERSION = 1


@dataclass(frozen=True)
class TankCircuit:
    """Series inductor into the sensor and its parasitic capacitance. Defaults: 820 nH, 0.5 pF, 50 Ω
    (f₀ ≈ 250 MHz, R_match ≈ 33 kΩ, a typical SET resistance). Ranges: L 100 nH–2 µH, C_p 0.2–2 pF."""

    inductance_h: float = 820e-9
    parasitic_capacitance_f: float = 0.5e-12
    z0_ohm: float = 50.0

    def __post_init__(self) -> None:
        if self.inductance_h <= 0 or self.parasitic_capacitance_f <= 0 or self.z0_ohm <= 0:
            raise ValueError("inductance, capacitance and line impedance must be positive")

    @property
    def resonance_hz(self) -> float:
        return float(1.0 / (2 * np.pi * np.sqrt(self.inductance_h * self.parasitic_capacitance_f)))

    @property
    def matched_resistance_ohm(self) -> float:
        """R at which the reflection vanishes at resonance: L / (C_p Z₀)."""
        return float(self.inductance_h / (self.parasitic_capacitance_f * self.z0_ohm))

    def impedance_ohm(self, f_hz: np.ndarray, r_sensor_ohm: np.ndarray) -> np.ndarray:
        """Complex input impedance; `f_hz` and `r_sensor_ohm` broadcast against each other."""
        w = 2 * np.pi * np.asarray(f_hz, dtype=float)
        r = np.asarray(r_sensor_ohm, dtype=float)
        y_load = 1j * w * self.parasitic_capacitance_f + 1.0 / r
        return 1j * w * self.inductance_h + 1.0 / y_load

    def reflection(self, f_hz: np.ndarray, r_sensor_ohm: np.ndarray) -> np.ndarray:
        """Complex reflection coefficient Γ."""
        z = self.impedance_ohm(f_hz, r_sensor_ohm)
        return (z - self.z0_ohm) / (z + self.z0_ohm)


def sensor_resistance_ohm(current_a: np.ndarray, bias_v: float, floor_ohm: float = 1e3, ceiling_ohm: float = 1e9) -> np.ndarray:
    """Small-bias resistance V_sd / I of the sensor, clipped to a plausible band."""
    i = np.maximum(np.abs(np.asarray(current_a, dtype=float)), 1e-18)
    return np.clip(bias_v / i, floor_ohm, ceiling_ohm)
