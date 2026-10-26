"""Noise parameters drawn per device: the fluctuator ensemble behind 1/f, the dominant switchers, and
the white readout floor, from catalog placeholders (entries #1 oxide traps and fluctuators, #12 slow
drift; `calibration: none` until written).

A device's charge noise is a draw, not a constant (ADR 0002): the 1/f amplitude varies between devices
by a factor of a few, a few dominant two-level fluctuators near the sensor dominate the telegraph
signature, and their rates spread over decades. This module draws those numbers; `readout` turns them
into processes (`SensedDoubleDot.from_noise_draw`), keeping the layer order (readout imports disorder,
never the reverse). The 1/f band reaches down to 0.1 mHz so that cuts repeated over hours drift.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MODEL_VERSION = 1


@dataclass(frozen=True)
class NoiseParams:
    """Placeholder distributions per material class."""

    one_over_f_amplitude_v_median: float = 150e-6
    """Median rms gate-referred 1/f amplitude over the band (V); log-normal spread `log_sigma`."""
    one_over_f_log_sigma: float = 0.4
    f_min_hz: float = 1e-4
    f_max_hz: float = 1e4
    n_fluctuators: int = 32
    switchers_mean: float = 1.0
    """Poisson mean number of dominant two-level fluctuators near the sensor."""
    switcher_amplitude_v_median: float = 3e-3
    switcher_amplitude_log_sigma: float = 0.5
    switcher_rate_log10_hz_range: tuple[float, float] = (-1.0, 2.0)
    """log10 of each switcher's up and down rates, uniform in this range (0.1 Hz to 100 Hz)."""
    white_sigma_a_at_1ms: float = 15e-12

    @classmethod
    def for_material(cls, material: str) -> NoiseParams:
        if material == "MOS":
            return cls(one_over_f_amplitude_v_median=300e-6, switchers_mean=2.0)
        return cls()


@dataclass(frozen=True)
class NoiseDraw:
    one_over_f_amplitude_v: float
    f_min_hz: float
    f_max_hz: float
    n_fluctuators: int
    switcher_amplitudes_v: tuple[float, ...] = ()
    switcher_rates_up_hz: tuple[float, ...] = ()
    switcher_rates_down_hz: tuple[float, ...] = ()
    white_sigma_a_at_1ms: float = 15e-12
    seed: int = 0

    def to_dict(self) -> dict:
        return {
            "one_over_f_amplitude_v": self.one_over_f_amplitude_v,
            "f_min_hz": self.f_min_hz,
            "f_max_hz": self.f_max_hz,
            "n_fluctuators": self.n_fluctuators,
            "switcher_amplitudes_v": list(self.switcher_amplitudes_v),
            "switcher_rates_up_hz": list(self.switcher_rates_up_hz),
            "switcher_rates_down_hz": list(self.switcher_rates_down_hz),
            "white_sigma_a_at_1ms": self.white_sigma_a_at_1ms,
            "seed": self.seed,
        }

    @classmethod
    def from_dict(cls, d: dict) -> NoiseDraw:
        return cls(
            float(d["one_over_f_amplitude_v"]),
            float(d["f_min_hz"]),
            float(d["f_max_hz"]),
            int(d["n_fluctuators"]),
            tuple(float(x) for x in d.get("switcher_amplitudes_v", [])),
            tuple(float(x) for x in d.get("switcher_rates_up_hz", [])),
            tuple(float(x) for x in d.get("switcher_rates_down_hz", [])),
            float(d.get("white_sigma_a_at_1ms", 15e-12)),
            int(d.get("seed", 0)),
        )


def draw_noise(params: NoiseParams, rng: np.random.Generator, seed: int, severity: float = 1.0) -> NoiseDraw:
    """One device's noise parameters. Severity scales the 1/f amplitude and the switcher count;
    severity 0 gives no 1/f and no switchers, leaving only the white floor."""
    if severity <= 0:
        return NoiseDraw(0.0, params.f_min_hz, params.f_max_hz, params.n_fluctuators, white_sigma_a_at_1ms=params.white_sigma_a_at_1ms, seed=seed)
    amp = params.one_over_f_amplitude_v_median * np.exp(rng.normal(0.0, params.one_over_f_log_sigma)) * severity
    n_sw = int(rng.poisson(params.switchers_mean * severity))
    amps = params.switcher_amplitude_v_median * np.exp(rng.normal(0.0, params.switcher_amplitude_log_sigma, n_sw))
    lo, hi = params.switcher_rate_log10_hz_range
    up = 10 ** rng.uniform(lo, hi, n_sw)
    down = 10 ** rng.uniform(lo, hi, n_sw)
    return NoiseDraw(
        float(amp), params.f_min_hz, params.f_max_hz, params.n_fluctuators,
        tuple(float(a) for a in amps), tuple(float(r) for r in up), tuple(float(r) for r in down),
        params.white_sigma_a_at_1ms, seed,
    )
