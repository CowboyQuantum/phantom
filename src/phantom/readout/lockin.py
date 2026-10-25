"""The `sensor_lockin` channel: a hardware derivative of the sensor signal.

A square wave of peak-to-peak amplitude `amplitude_v` at `frequency_hz`, applied along a direction in
gate space (one plunger, both in common mode, or the two in opposition for detuning), and the sensor
current demodulated at that frequency: the output at a point is the average over `cycles` cycles of
(I(V + d·A/2) − I(V − d·A/2)) / A, in amperes per volt of modulation. It differs from a numerical
derivative of a map in two physical ways: the finite modulation amplitude smears features narrower than
A, and noise at and near the modulation frequency survives the demodulation while slow drift cancels.
Each point costs 2·cycles sensor evaluations.

`lockin_sweep_a_per_v` is the quasi-static channel: the charge follows the chop instantly.
`lockin_sweep_dynamic_a_per_v` adds the dot-to-lead tunnel rates (`LeadRateLandscape`): where the
rate is below the chop frequency the charge cannot follow and the derivative signal latches or
vanishes; where the rate is far above it the transition is lifetime-broadened and the derivative
peak is lower and wider. The charge present is carried from one half-cycle to the next and from one
point to the next in scan order, so latching also shows along the fast axis.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from phantom.readout.loading import LeadRateLandscape, latched_occupation
from phantom.readout.measure import SensedDoubleDot

MODEL_VERSION = 2
"""Bumped whenever any output of this module changes for the same parameters and seed (ADR 0003)."""


@dataclass(frozen=True)
class LockinParams:
    """Square-wave modulation along `direction` (gate name → weight; a single `gate` means weight 1)
    and demodulation of the sensor current. The modulation on gate g is weight_g · amplitude_v."""

    gate: str | None = None
    direction: dict[str, float] = field(default_factory=dict)
    amplitude_v: float = 2e-3
    frequency_hz: float = 1e3
    cycles: int = 4

    def __post_init__(self) -> None:
        if self.amplitude_v <= 0 or self.frequency_hz <= 0 or self.cycles < 1:
            raise ValueError("amplitude_v and frequency_hz must be positive, cycles >= 1")
        if not self.direction and self.gate is None:
            raise ValueError("give a gate or a direction")
        if not self.direction:
            object.__setattr__(self, "direction", {self.gate: 1.0})

    @property
    def evaluations_per_point(self) -> int:
        return 2 * self.cycles

    @property
    def dwell_s(self) -> float:
        """Time spent per point: `cycles` periods of the modulation."""
        return self.cycles / self.frequency_hz

    def vector_v(self, gate_names: tuple[str, ...]) -> np.ndarray:
        """Peak-to-peak modulation per gate (gates,) in volts."""
        d = np.zeros(len(gate_names))
        for g, w in self.direction.items():
            d[gate_names.index(g)] = w * self.amplitude_v
        return d


def _modulated(device: SensedDoubleDot, v: np.ndarray, lockin: LockinParams):
    v = np.asarray(v, dtype=float)
    shape, n_g = v.shape[:-1], v.shape[-1]
    n_pts = int(np.prod(shape))
    d = lockin.vector_v(device.charge.gate_names)
    half = np.array([+0.5, -0.5])
    v_mod = np.repeat(v.reshape(n_pts, 1, 1, n_g), lockin.cycles, axis=1).repeat(2, axis=2).copy()
    v_mod += half[None, None, :, None] * d[None, None, None, :]
    return shape, n_pts, v_mod


def _times(lockin: LockinParams, n_pts: int, t0_s: float) -> np.ndarray:
    period = 1.0 / lockin.frequency_hz
    return (
        t0_s
        + np.arange(n_pts)[:, None, None] * lockin.dwell_s
        + np.arange(lockin.cycles)[None, :, None] * period
        + np.array([0.0, 0.5])[None, None, :] * period
    )


def lockin_sweep_a_per_v(
    device: SensedDoubleDot,
    v: np.ndarray,
    lockin: LockinParams,
    t0_s: float | None = None,
    v_sd: float | None = None,
) -> np.ndarray:
    """Quasi-static demodulated sensor response (A/V) for gate voltages `v` (..., gates). Points are
    visited in flattened order starting at `t0_s`, each dwelling `lockin.dwell_s`; noise-free when
    `t0_s` is None."""
    shape, n_pts, v_mod = _modulated(device, v, lockin)
    times = None if t0_s is None else _times(lockin, n_pts, t0_s)
    i = device.sweep_a(v_mod, times, v_sd)  # (n_pts, cycles, 2)
    out = (i[..., 0] - i[..., 1]).mean(axis=1) / lockin.amplitude_v
    return out.reshape(shape)


def lockin_sweep_dynamic_a_per_v(
    device: SensedDoubleDot,
    v: np.ndarray,
    lockin: LockinParams,
    rates: LeadRateLandscape,
    rng: np.random.Generator,
    t0_s: float | None = None,
    v_sd: float | None = None,
) -> dict[str, np.ndarray]:
    """Demodulated sensor response with tunnel-rate dynamics. Returns `signal_a_per_v` (...,), the
    `log10_rate_hz` of the faster dot at each point (...,), and `n_present_mean` (..., 2), the charge
    actually present averaged over the point's half-cycles."""
    shape, n_pts, v_mod = _modulated(device, v, lockin)
    steps = n_pts * lockin.cycles * 2
    flat = v_mod.reshape(steps, -1)
    half_s = 0.5 / lockin.frequency_hz
    if t0_s is None:
        g = s = r = np.zeros(steps)
    else:
        times = _times(lockin, n_pts, t0_s).reshape(-1)
        g, s, r = device.gate_noise.trace(times), device.sensor_noise.trace(times), device.readout_noise.trace(times)
    v_eff = flat + g[:, None]
    gamma = rates.rate_hz(v_eff)  # (steps, 2)
    kt = rates.broadened_kt_ev(v_eff, device.charge.temperature_k)
    n_eq = device.charge.mean_occupation(v_eff, kt)
    n_present = latched_occupation(n_eq, gamma, half_s, rng, per_dot=True)
    i = device.sensor.current_a(device.sensor.potential(v_eff, n_present, s), v_sd) + r
    i = i.reshape(n_pts, lockin.cycles, 2)
    signal = (i[..., 0] - i[..., 1]).mean(axis=1) / lockin.amplitude_v
    log_rate = np.log10(gamma.max(axis=-1).reshape(n_pts, -1).mean(axis=1))
    n_mean = n_present.reshape(n_pts, -1, 2).mean(axis=1)
    return {
        "signal_a_per_v": signal.reshape(shape),
        "log10_rate_hz": log_rate.reshape(shape),
        "n_present_mean": n_mean.reshape(shape + (2,)),
    }
