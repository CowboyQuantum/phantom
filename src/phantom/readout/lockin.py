"""The `sensor_lockin` channel: a hardware derivative of the sensor signal.

A square wave of peak-to-peak amplitude `amplitude_v` at `frequency_hz` is added to one named gate and
the sensor current is demodulated at that frequency: the output at a point is the average over
`cycles` cycles of (I(V + A/2) − I(V − A/2)) / A, in amperes per volt. It differs from a numerical
derivative of a map in two physical ways: the finite modulation amplitude smears features narrower
than A, and noise at and near the modulation frequency survives the demodulation while slow drift
cancels. Each point costs 2·cycles sensor evaluations.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phantom.readout.measure import SensedDoubleDot

MODEL_VERSION = 1
"""Bumped whenever any output of this module changes for the same parameters and seed (ADR 0003)."""


@dataclass(frozen=True)
class LockinParams:
    """Square-wave modulation of `gate` and demodulation of the sensor current."""

    gate: str
    amplitude_v: float = 2e-3
    frequency_hz: float = 1e3
    cycles: int = 4

    def __post_init__(self) -> None:
        if self.amplitude_v <= 0 or self.frequency_hz <= 0 or self.cycles < 1:
            raise ValueError("amplitude_v and frequency_hz must be positive, cycles >= 1")

    @property
    def evaluations_per_point(self) -> int:
        return 2 * self.cycles

    @property
    def dwell_s(self) -> float:
        """Time spent per point: `cycles` periods of the modulation."""
        return self.cycles / self.frequency_hz


def lockin_sweep_a_per_v(
    device: SensedDoubleDot,
    v: np.ndarray,
    lockin: LockinParams,
    t0_s: float | None = None,
    v_sd: float | None = None,
) -> np.ndarray:
    """Demodulated sensor response (A/V) for gate voltages `v` (..., gates). Points are visited in
    flattened order starting at `t0_s`, each dwelling `lockin.dwell_s`; noise-free when `t0_s` is None."""
    v = np.asarray(v, dtype=float)
    shape, n_g = v.shape[:-1], v.shape[-1]
    j = device.charge.gate_names.index(lockin.gate)
    n_pts = int(np.prod(shape))
    half = np.array([+0.5, -0.5]) * lockin.amplitude_v
    v_mod = np.repeat(v.reshape(n_pts, 1, 1, n_g), lockin.cycles, axis=1).repeat(2, axis=2).copy()
    v_mod[..., j] += half[None, None, :]
    times = None
    if t0_s is not None:
        period = 1.0 / lockin.frequency_hz
        times = (
            t0_s
            + np.arange(n_pts)[:, None, None] * lockin.dwell_s
            + np.arange(lockin.cycles)[None, :, None] * period
            + np.array([0.0, 0.5])[None, None, :] * period
        )
    i = device.sweep_a(v_mod, times, v_sd)  # (n_pts, cycles, 2)
    out = (i[..., 0] - i[..., 1]).mean(axis=1) / lockin.amplitude_v
    return out.reshape(shape)
