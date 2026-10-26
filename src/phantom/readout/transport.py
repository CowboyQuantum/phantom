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


# --- shapes derived from a layout's potential (Step 2) ------------------------------------------


def channel_transmission_current_a(
    u_barrier_ev: np.ndarray, threshold_ev: float = -0.05, width_ev: float = 0.004, i_sat_a: float = 2e-9
) -> np.ndarray:
    """Current through a channel whose highest point of electron energy along it is `u_barrier_ev`:
    open (saturated) when that energy lies below `threshold_ev`, closed above, with a soft step of width
    `width_ev`. The barrier energy comes from the layout's potential, so thresholds and slopes against
    every gate are derived; the step shape itself is phenomenological (calibration none)."""
    return i_sat_a / (1.0 + np.exp((np.asarray(u_barrier_ev, dtype=float) - threshold_ev) / width_ev))


def gate_leakage_current_a(
    v_gate: np.ndarray,
    area_nm2: float,
    thickness_nm: float,
    onset_field_v_per_nm: float = 0.08,
    steepness_v: float = 0.2,
    i_onset_a_per_um2: float = 1e-9,
) -> np.ndarray:
    """Leakage through the gate dielectric: antisymmetric in the voltage, proportional to the gate area,
    and steeply exponential (scale `steepness_v`) so that it reaches `i_onset_a_per_um2` per square
    micrometre at the onset voltage onset_field × thickness and is negligible well inside that window.
    Silicon dioxide breaks down near 1 V/nm; the onset field is set well below it. Phenomenological;
    calibration none."""
    v_onset = onset_field_v_per_nm * thickness_nm
    v = np.asarray(v_gate, dtype=float)
    return i_onset_a_per_um2 * (area_nm2 * 1e-6) * np.sinh(v / steepness_v) / np.sinh(v_onset / steepness_v)


def ohmic_current_a(v_sd: np.ndarray, series_resistance_ohm: float = 20e3) -> np.ndarray:
    """Source–drain current through an open channel with its contacts: linear, one series resistance
    (contacts plus the 2DEG squares). Phenomenological; calibration none."""
    return np.asarray(v_sd, dtype=float) / series_resistance_ohm
