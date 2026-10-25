"""Noise processes with an explicit time axis.

Seeded from qdt-sim's `noise.py` (Step 1); the sample sequences are bit-identical to it for the same
seed and times. Every process is sampled at non-decreasing times `t` in seconds and returns a value in
the unit of whatever it is added to: volts of gate-equivalent offset, or amperes of readout noise.
Processes are stateful and seeded, so a run is reproducible for a given seed (ADR 0003).

- `WhiteNoise`: memoryless Gaussian; readout noise per sample.
- `TelegraphNoise`: one two-level fluctuator with up and down rates; a single dominant charge trap.
  Advancing by Δt uses the exact two-state Markov propagator, so any Δt is sampled correctly.
- `OneOverFNoise`: an ensemble of symmetric fluctuators with log-uniform rates in [f_min, f_max],
  whose summed spectrum is ~1/f in that band; slow charge noise that drifts the honeycomb and the
  sensor working point.
- `Composite`: sum of processes. `NoNoise`: zero.

`sample(t)` returns one value; `trace(times)` returns the array the same `sample` calls would, faster.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable

import numpy as np

MODEL_VERSION = 1
"""Bumped whenever any output of this module changes for the same parameters and seed (ADR 0003)."""


@runtime_checkable
class NoiseProcess(Protocol):
    def sample(self, t: float) -> float:
        """Value at time `t` (s); `t` never decreases between calls."""
        ...

    def trace(self, times: np.ndarray) -> np.ndarray:
        """Values at non-decreasing `times` (s), identical to successive `sample` calls."""
        ...


@dataclass
class WhiteNoise:
    """Gaussian noise of standard deviation `sigma` (unit of the quantity it is added to)."""

    sigma: float
    seed: int = 0
    _rng: np.random.Generator = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self._rng = np.random.default_rng(self.seed)

    def sample(self, t: float) -> float:
        return float(self._rng.normal(0.0, self.sigma)) if self.sigma else 0.0

    def trace(self, times: np.ndarray) -> np.ndarray:
        n = len(times)
        return self._rng.normal(0.0, self.sigma, n) if self.sigma else np.zeros(n)


@dataclass
class TelegraphNoise:
    """Two-level fluctuator: value ±amplitude/2, switching with `rate_up_hz` (0→1) and `rate_down_hz` (1→0)."""

    amplitude: float
    rate_up_hz: float
    rate_down_hz: float
    seed: int = 0
    _rng: np.random.Generator = field(init=False, repr=False)
    _state: int = field(init=False, repr=False)
    _t: float | None = field(init=False, default=None, repr=False)

    def __post_init__(self) -> None:
        if self.rate_up_hz < 0 or self.rate_down_hz < 0 or self.rate_up_hz + self.rate_down_hz <= 0:
            raise ValueError("rates must be non-negative with a positive sum")
        self._rng = np.random.default_rng(self.seed)
        p_up = self.rate_up_hz / (self.rate_up_hz + self.rate_down_hz)
        self._state = int(self._rng.random() < p_up)

    def _advances(self, t: float) -> bool:
        """True when moving to `t` is a step forward that consumes one random number."""
        return self._t is not None and t > self._t

    def _advance(self, t: float, u: float) -> None:
        total = self.rate_up_hz + self.rate_down_hz
        p_inf = self.rate_up_hz / total
        decay = np.exp(-total * (t - self._t))
        self._state = int(u < p_inf + (self._state - p_inf) * decay)

    def sample(self, t: float) -> float:
        if self._advances(t):
            self._advance(t, self._rng.random())
        self._t = t
        return self.amplitude * (self._state - 0.5)

    def trace(self, times: np.ndarray) -> np.ndarray:
        """Identical to successive `sample` calls: one random number per forward step, drawn in a block."""
        times = np.asarray(times, dtype=float)
        n = len(times)
        if n == 0:
            return np.zeros(0)
        prev = np.concatenate([[np.nan if self._t is None else self._t], times[:-1]])
        steps = ~np.isnan(prev) & (times > prev)
        u = self._rng.random(int(steps.sum()))
        out = np.empty(n)
        k = 0
        for i in range(n):
            t = float(times[i])
            if steps[i]:
                self._advance(t, u[k])
                k += 1
            self._t = t
            out[i] = self._state
        return self.amplitude * (out - 0.5)


@dataclass
class OneOverFNoise:
    """Sum of `n_fluctuators` symmetric telegraph processes with rates log-uniform in [f_min_hz, f_max_hz],
    scaled so the total rms is `amplitude`. The spectrum is ~1/f between the two rates."""

    amplitude: float
    f_min_hz: float = 1e-2
    f_max_hz: float = 1e4
    n_fluctuators: int = 24
    seed: int = 0
    _fluctuators: list[TelegraphNoise] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not 0 < self.f_min_hz < self.f_max_hz:
            raise ValueError("need 0 < f_min_hz < f_max_hz")
        rng = np.random.default_rng(self.seed)
        rates = np.exp(rng.uniform(np.log(self.f_min_hz), np.log(self.f_max_hz), self.n_fluctuators))
        a = 2.0 * self.amplitude / np.sqrt(self.n_fluctuators)  # each ±a/2 → variance a²/4
        self._fluctuators = [
            TelegraphNoise(a, float(r) / 2, float(r) / 2, seed=int(rng.integers(2**31))) for r in rates
        ]

    def sample(self, t: float) -> float:
        return float(sum(f.sample(t) for f in self._fluctuators))

    def trace(self, times: np.ndarray) -> np.ndarray:
        times = np.asarray(times, dtype=float)
        total = np.zeros(len(times))
        for f in self._fluctuators:
            total += f.trace(times)
        return total


@dataclass
class Composite:
    """Sum of processes, sampled at the same time."""

    processes: tuple[NoiseProcess, ...]

    def sample(self, t: float) -> float:
        return float(sum(p.sample(t) for p in self.processes))

    def trace(self, times: np.ndarray) -> np.ndarray:
        total = np.zeros(len(times))
        for p in self.processes:
            total += p.trace(times)
        return total


class NoNoise:
    """Always zero."""

    def sample(self, t: float) -> float:
        return 0.0

    def trace(self, times: np.ndarray) -> np.ndarray:
        return np.zeros(len(times))
