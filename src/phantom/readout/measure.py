"""A double dot read by an SET sensor, with noise on an explicit time axis.

This is the Step 1 composition of `charge` and `readout`: the constant-interaction model gives the
dot occupation, the sensor turns it into a current, and three noise processes enter where they do in
a real measurement (as in qdt-sim's `SimulatedDevice`, minus any notion of channels or instruments):

- `gate`: a common-mode gate-voltage offset (V) added to every gate, slow charge noise that drifts
  the honeycomb;
- `sensor`: an offset of the sensor potential φ_s (V-equivalent), a trap near the sensor;
- `readout`: additive current noise (A).

Every point of a sweep has a time; the noise is sampled at those times in order. With no times given,
the sweep is noise-free and is the pure function of the gate vector that the P0 tests use.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from phantom.charge.ci import ConstantInteractionParams
from phantom.disorder.noise_draw import NoiseDraw
from phantom.readout.noise import (
    Composite,
    NoiseProcess,
    NoNoise,
    OneOverFNoise,
    TelegraphNoise,
    WhiteNoise,
)
from phantom.readout.sensor import SensorSETParams

MODEL_VERSION = 1
"""Bumped whenever any output of this module changes for the same parameters and seed (ADR 0003)."""


@dataclass
class SensedDoubleDot:
    """A CI double dot, its SET sensor and the noise processes, built for a given seed."""

    charge: ConstantInteractionParams = field(default_factory=ConstantInteractionParams)
    sensor: SensorSETParams = field(default_factory=SensorSETParams)
    gate_noise: NoiseProcess = field(default_factory=NoNoise)
    sensor_noise: NoiseProcess = field(default_factory=NoNoise)
    readout_noise: NoiseProcess = field(default_factory=NoNoise)

    def __post_init__(self) -> None:
        if self.sensor.alpha_v_per_v.shape != (len(self.charge.gate_names),):
            raise ValueError("sensor.alpha_v_per_v needs one entry per charge gate")

    @classmethod
    def default_noisy(cls, seed: int = 0, **kwargs) -> SensedDoubleDot:
        """qdt-sim's default noise: 150 µV rms 1/f on the gates between 10 mHz and 10 kHz, one 4 mV
        two-level fluctuator at the sensor switching at 0.5 Hz each way, 15 pA white readout noise."""
        return cls(
            gate_noise=OneOverFNoise(amplitude=150e-6, f_min_hz=1e-2, f_max_hz=1e4, seed=seed),
            sensor_noise=TelegraphNoise(amplitude=4e-3, rate_up_hz=0.5, rate_down_hz=0.5, seed=seed + 1),
            readout_noise=WhiteNoise(sigma=15e-12, seed=seed + 2),
            **kwargs,
        )

    @classmethod
    def from_noise_draw(cls, draw: NoiseDraw, sample_dt_s: float = 1e-3, switcher_scale: float = 1.0, **kwargs) -> SensedDoubleDot:
        """The processes a per-device noise draw (disorder.noise_draw) describes: a 1/f ensemble on the
        gates, the drawn switchers summed at the sensor with their amplitudes times `switcher_scale`
        (the draw states them in the typed sensor's units; a sensor dot converts them, ADR 0026), and the
        white floor rescaled to `sample_dt_s` (constant spectral density). Seeds derive from the draw's seed."""
        gate = OneOverFNoise(draw.one_over_f_amplitude_v, draw.f_min_hz, draw.f_max_hz, draw.n_fluctuators, seed=draw.seed) if draw.one_over_f_amplitude_v > 0 else NoNoise()
        switchers = tuple(
            TelegraphNoise(a * switcher_scale, up, down, seed=draw.seed + 10 + k)
            for k, (a, up, down) in enumerate(zip(draw.switcher_amplitudes_v, draw.switcher_rates_up_hz, draw.switcher_rates_down_hz, strict=True))
        )
        sensor = Composite(switchers) if switchers else NoNoise()
        white = WhiteNoise(draw.white_sigma_a_at_1ms * np.sqrt(1e-3 / sample_dt_s), seed=draw.seed + 2)
        return cls(gate_noise=gate, sensor_noise=sensor, readout_noise=white, **kwargs)

    # --- noise-free physics ----------------------------------------------------------------

    def current_a(
        self,
        v: np.ndarray,
        v_sd: float | np.ndarray | None = None,
        delta_phi_v: float | np.ndarray = 0.0,
    ) -> np.ndarray:
        """Sensor current (A) at gate voltages `v` (..., gates) with the thermal-average occupation."""
        n = self.charge.mean_occupation(v)
        return self.sensor.current_a(self.sensor.potential(v, n, delta_phi_v), v_sd, v)

    # --- with time and noise -----------------------------------------------------------------

    def sweep_a(
        self,
        v: np.ndarray,
        times_s: np.ndarray | None = None,
        v_sd: float | np.ndarray | None = None,
    ) -> np.ndarray:
        """Sensor current (A) for gate voltages `v` (..., gates) sampled at `times_s` (...) in the order
        the flattened points are listed; noise-free when `times_s` is None."""
        v = np.asarray(v, dtype=float)
        if times_s is None:
            return self.current_a(v, v_sd)
        shape = v.shape[:-1]
        times = np.asarray(times_s, dtype=float).reshape(-1)
        if times.shape != (int(np.prod(shape)),):
            raise ValueError("times_s must have one time per gate point")
        g = self.gate_noise.trace(times)
        s = self.sensor_noise.trace(times)
        r = self.readout_noise.trace(times)
        v_eff = v.reshape(-1, v.shape[-1]) + g[:, None]
        v_sd_flat = None if v_sd is None else np.broadcast_to(np.asarray(v_sd, dtype=float), shape).reshape(-1)
        i = self.current_a(v_eff, v_sd_flat, s) + r
        return i.reshape(shape)

    def trace_a(self, v: np.ndarray, times_s: np.ndarray, v_sd: float | None = None) -> np.ndarray:
        """Sensor current (A) at one fixed gate point `v` (gates,) over `times_s` (n,)."""
        v = np.asarray(v, dtype=float)
        return self.sweep_a(np.broadcast_to(v, (len(times_s), v.shape[-1])), times_s, v_sd)


def raster_times_s(shape: tuple[int, ...], dt_s: float, t0_s: float = 0.0) -> np.ndarray:
    """Times of a raster sweep of `shape` points, one every `dt_s` seconds in flattened order."""
    n = int(np.prod(shape))
    return (t0_s + dt_s * np.arange(n)).reshape(shape)
