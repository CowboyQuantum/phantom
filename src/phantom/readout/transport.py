"""Channel turn-on and barrier pinch-off: the direct-current shapes of the electrical-health protocols.

Phenomenological closed forms for Step 1, so the `turn-on-single-gate` and `pinch-off-barrier`
protocols have a shape the owner can review. Step 2 replaces both with currents derived from the
layout's potential and marks these as superseded. Neither is a validated model: `calibration: none`.

- Turn-on: a channel under an accumulation gate conducts once the gate passes a threshold, with a
  soft knee of width `width_v`, and saturates: I = I_sat · softplus((V − V_th)/w) / softplus(V_span/w).
- Pinch-off: a barrier gate closes a channel; the current falls exponentially below the pinch-off
  voltage with the barrier's own `v_scale_v` and saturates above it: I = I_sat · t / (t + t_sat), with
  t the `Barrier` tunnel coupling.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phantom.charge.dot import Barrier

MODEL_VERSION = 1
"""Bumped whenever any output of this module changes for the same parameters (ADR 0003)."""


def _softplus(x: np.ndarray) -> np.ndarray:
    return np.logaddexp(0.0, x)


@dataclass(frozen=True)
class ChannelTurnOn:
    """A gated channel that turns on at `threshold_v` over `width_v` and saturates at `i_sat_a`
    once the gate is `span_v` past threshold. Range: thresholds 0.2–1.5 V, widths 10–100 mV."""

    threshold_v: float = 0.6
    width_v: float = 0.05
    span_v: float = 0.5
    i_sat_a: float = 5e-9

    def current_a(self, v_gate: float | np.ndarray) -> np.ndarray:
        x = (np.asarray(v_gate, dtype=float) - self.threshold_v) / self.width_v
        return self.i_sat_a * np.tanh(_softplus(x) / _softplus(self.span_v / self.width_v))


def pinch_off_current_a(barrier: Barrier, v_b: float | np.ndarray, i_sat_a: float = 2e-9) -> np.ndarray:
    """Current through a channel closed by `barrier`: exponential pinch-off below `barrier.v_ref_v`,
    saturation at `i_sat_a` well above it. The tunnel coupling at `v_ref_v` sets the half-current point."""
    t = barrier.tunnel_coupling_ev(v_b)
    return i_sat_a * t / (t + barrier.t0_ev)
