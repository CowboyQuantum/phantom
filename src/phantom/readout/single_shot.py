"""Single-shot charge readout: integrated sensor signals, their histograms, and signal-to-noise against
integration time.

At a readout point the double dot is in one of two charge states (the two the readout must tell apart);
the sensor gives a different current for each. A shot integrates the noisy current for `tau_s`. Shots
follow one another in time on the device's noise processes, so slow 1/f drift and telegraph switching
move the two distributions between shots as they do in a real histogram; white readout noise has a fixed
spectral density, so its per-sample standard deviation scales as 1/√dt and its contribution to a shot
falls as 1/√tau. The signal-to-noise ratio is the separation of the two means over their pooled
standard deviation; the readout fidelity is the fraction of shots on the right side of the midpoint
threshold (ADR 0018). Placeholder noise amplitudes (`calibration: none`).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phantom.readout.measure import SensedDoubleDot
from phantom.readout.noise import OneOverFNoise, TelegraphNoise, WhiteNoise

MODEL_VERSION = 1

REFERENCE_SAMPLE_S = 1e-3
"""The sampling interval at which the noise defaults (15 pA white per sample) are stated."""


@dataclass(frozen=True)
class ShotResult:
    signals_a: np.ndarray
    """Integrated current per shot (n_shots,)."""
    true_state: np.ndarray
    """0 or 1 per shot."""
    means_a: tuple[float, float]
    stds_a: tuple[float, float]

    @property
    def snr(self) -> float:
        pooled = np.sqrt(0.5 * (self.stds_a[0] ** 2 + self.stds_a[1] ** 2))
        return float(abs(self.means_a[1] - self.means_a[0]) / max(pooled, 1e-30))

    @property
    def fidelity(self) -> float:
        threshold = 0.5 * (self.means_a[0] + self.means_a[1])
        decided = (self.signals_a > threshold).astype(int) if self.means_a[1] > self.means_a[0] else (self.signals_a < threshold).astype(int)
        return float((decided == self.true_state).mean())


def noisy_device_for_sampling(dev: SensedDoubleDot, dt_s: float, seed: int) -> SensedDoubleDot:
    """The same device with white readout noise rescaled to the sampling interval `dt_s` (constant
    spectral density), fresh seeds for the shot sequence."""
    white = dev.readout_noise
    sigma_ref = white.sigma if isinstance(white, WhiteNoise) else 15e-12
    sigma = sigma_ref * np.sqrt(REFERENCE_SAMPLE_S / dt_s)
    gate = dev.gate_noise
    sensor = dev.sensor_noise
    gate_new = OneOverFNoise(gate.amplitude, gate.f_min_hz, gate.f_max_hz, gate.n_fluctuators, seed=seed) if isinstance(gate, OneOverFNoise) else gate
    sensor_new = TelegraphNoise(sensor.amplitude, sensor.rate_up_hz, sensor.rate_down_hz, seed=seed + 1) if isinstance(sensor, TelegraphNoise) else sensor
    return SensedDoubleDot(dev.charge, dev.sensor, gate_new, sensor_new, WhiteNoise(sigma, seed=seed + 2))


def single_shots(
    dev: SensedDoubleDot,
    v: np.ndarray,
    states: np.ndarray,
    tau_s: float,
    dt_s: float,
    n_shots: int,
    rng: np.random.Generator,
    p_state1: float = 0.5,
    t0_s: float = 0.0,
) -> ShotResult:
    """Integrate `n_shots` consecutive shots of length `tau_s` sampled every `dt_s` at gate point `v`
    (gates,), the charge in one of the two `states` (2, 2) drawn with probability `p_state1` for the
    second. `dev` should come from `noisy_device_for_sampling` for the chosen `dt_s`."""
    v = np.asarray(v, dtype=float)
    n_per = max(round(tau_s / dt_s), 1)
    true_state = (rng.random(n_shots) < p_state1).astype(int)
    times = t0_s + dt_s * np.arange(n_shots * n_per)
    g = dev.gate_noise.trace(times)
    s = dev.sensor_noise.trace(times)
    r = dev.readout_noise.trace(times)
    v_eff = v[None, :] + g[:, None]
    n = np.asarray(states, dtype=float)[np.repeat(true_state, n_per)]
    i = dev.sensor.current_a(dev.sensor.potential(v_eff, n, s), None, v_eff) + r
    signals = i.reshape(n_shots, n_per).mean(axis=1)
    means = tuple(float(signals[true_state == k].mean()) if np.any(true_state == k) else float("nan") for k in (0, 1))
    stds = tuple(float(signals[true_state == k].std()) if np.sum(true_state == k) > 1 else float("nan") for k in (0, 1))
    return ShotResult(signals, true_state, means, stds)
