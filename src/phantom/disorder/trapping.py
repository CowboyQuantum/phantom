"""Charge trapping with gate history: the hysteresis loop in a gate sweep.

Catalog entry `oxide-traps-fluctuators` (#1) at its slowest: traps in the dielectric that charge when
the gate is driven high and release only when it is driven low again. The trapped charge screens the
gate, so the same programmed voltage produces a different channel energy on the way up and on the way
down, and a turn-on sweep up then down opens a loop. This is the static, history-dependent partner of
the telegraph and 1/f processes in `readout.noise`.

Model: a trapped fraction f in [0, 1] follows the gate voltage with first-order kinetics per sample:
above `capture_v` it relaxes toward 1 at rate `capture_rate_per_sample`, below `release_v` toward 0 at
`release_rate_per_sample`, and between the two it holds. The effective gate voltage is V − f · shift_v.
Severity multiplies `shift_v`; severity 0 gives no shift and the two sweep directions coincide bit for bit.
Placeholder parameters, `calibration: none`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MODEL_VERSION = 1


@dataclass(frozen=True)
class TrappingParams:
    capture_v: float = 0.6
    release_v: float = 0.0
    shift_v: float = 0.15
    """Threshold shift when fully trapped (V of gate-equivalent)."""
    capture_rate_per_sample: float = 0.05
    release_rate_per_sample: float = 0.05

    def __post_init__(self) -> None:
        if self.release_v > self.capture_v or self.shift_v < 0:
            raise ValueError("release_v <= capture_v and shift_v >= 0")
        if not (0 <= self.capture_rate_per_sample <= 1 and 0 <= self.release_rate_per_sample <= 1):
            raise ValueError("rates per sample are fractions in [0, 1]")


def trapped_fraction(v_sequence: np.ndarray, params: TrappingParams, f0: float = 0.0) -> np.ndarray:
    """Trapped fraction after each sample of a gate-voltage sequence visited in order."""
    v = np.asarray(v_sequence, dtype=float)
    f = np.empty(len(v))
    cur = float(f0)
    for i, vi in enumerate(v):
        if vi >= params.capture_v:
            cur += params.capture_rate_per_sample * (1.0 - cur)
        elif vi <= params.release_v:
            cur -= params.release_rate_per_sample * cur
        f[i] = cur
    return f


def effective_gate_voltage(v_sequence: np.ndarray, params: TrappingParams, severity: float = 1.0, f0: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """(effective voltage, trapped fraction) along a sequence: V − severity · shift_v · f."""
    v = np.asarray(v_sequence, dtype=float)
    if severity <= 0:
        return v.copy(), np.zeros(len(v))
    f = trapped_fraction(v, params, f0)
    return v - severity * params.shift_v * f, f
