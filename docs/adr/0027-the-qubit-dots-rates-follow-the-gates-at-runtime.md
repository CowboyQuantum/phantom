# 0027 — The qubit dots' rates follow the gates at runtime: compiled slopes of the WKB action at a fixed tunnelling level, capped at the attempt rate

**Status:** Accepted
**Date:** 2026-10-30
**Step:** 8 follow-up (issue #53, after ADR 0026)

## Context
Compile computes the double dot's three rates (dot 1 to the left lead, dot 2 to the right lead, and the
transfer between the dots) once, at the operating point. Each is the attempt rate times exp(−2S), with
S the WKB action (ADR 0017, 0023). At runtime they were constants. `State` reported the same Γ_L,
Γ_R and Γ_m at every gate point, and the master equation used them everywhere. On the reference
device they are 10⁹–10¹¹ Hz, so every dwell from a microsecond up saw equilibrium. No barrier gate
could slow a lead into the range where latching or telegraph switching shows, although the rate
equation and the Markov trajectories that would show them were already built.

ADR 0026 gave the sensor's own leads gate-dependent rates: slopes d ln Γ/dV against every gate at a
tunnelling level held fixed, and a logistic that caps the rate at the attempt rate. The owner asked
for the same for the qubit dots (2026-10-30: "gates to modulate those dynamics").

The closed-form potential has no screening (#48). On seed 2 the slopes of ln Γ_L reach 213/V on
P1, 173/V on R_L and 119/V on B_L. The slope of Γ_m on B_M is 232/V, which equals 2/v_scale of the
interdot coupling t. That is expected: t ∝ √T and Γ_m ∝ T come from the same cut.

## Decision
The dots' lead and interdot rates follow the gates at runtime. The rule is the one ADR 0026 set for
the sensor's leads:

- Compile places each barrier's tunnelling level once at the operating point, by
  `tunnelling_level_ev`, and holds it fixed.
- It then takes the central differences (±10 mV) of the action against every gate at that level.
- At runtime the exponent moves linearly with those slopes, and `saturating_rates_hz`, shared with
  the sensor's leads, anchors the rate at the compiled one and caps it at the dot's attempt rate.

Specifically:

- **Record.** `CompiledBarrier` gains `attempt_rate_hz`, `tunnelling_level_ev` and
  `log_rate_slopes_per_v`. The compiled rates are unchanged bit for bit. `DEVICE_FORMAT_VERSION` 3.
- **Runtime.** `readout.rates.GateRates` holds the reference point, the rates, the attempt rates and
  the slopes. `RateEquation.gate_rates` is optional: when it is given, each rate matrix takes its
  three rates at its own voltages, and without it the rates are the constants as before.
  `Device.runtime_rates()` always gives it.
- **Consumer surface.** `State.gamma_left_hz`, `gamma_right_hz` and `gamma_interdot_hz` are arrays
  with the input's leading axes, at the carried voltages. `Truth` keeps the compiled values at the
  operating point. A session's dynamics use the rates at each sample's voltages, gate noise included.
- **Equilibrium.** Only the prefactors move, so detailed balance and the Boltzmann stationary state
  are exact at every point (P0).
- **Release.** This is folded into the untagged 0.3.0, as the owner decided.
  Model versions: `compile.device` 4, `device` 4, `device.session` 3.

## Consequences
- Easier:
  - A consumer closes B_L or B_R and slows the lead into the kilohertz range. On seed 2, Γ_L at
    dot 1's transition is about 10⁴ Hz with B_L 0.1 V lower.
  - A sweep there latches. At 3 µs per point the transition opens into a loop about 3 mV wide.
  - A parked `sampled` session telegraphs from the dot's own lead, and closing B_M slows interdot
    transfer.
  - The gallery shows all of this (`session-closing-barrier`).
- Harder:
  - The rates inherit the unscreened sensitivity: a factor e every 4–10 mV, P1 moves Γ_L more than
    B_L does, and a reservoir gate reaches the lead barrier strongly (#48).
  - Across one P1 transition Γ_L changes by about e⁸. Within ±10 mV the runtime exponent stays within
    20 % of the action's change plus the grid's step (about 0.1 in ln Γ). Beyond that the linear
    exponent drifts from the compile tier.
  - The equilibrium stays exact, but the time scales far from the operating point are an
    extrapolation.
- Must now:
  - `State`'s three rates are arrays. A consumer that cached a scalar rate reads it per point.
  - Device files from earlier 0.3.0 commits do not load; recompile them from their seed.
  - The contract's hold check settles for the slower of its two points.
  - `session-hold-relaxation` overlays the total rate out of the held state at the read point. Its
    earlier overlay, exp(−Γ_L t), left out the interdot path that empties through dot 2's lead. That
    path made the mean relax twice as fast on the gallery device before this change too.

## Alternatives considered
- **The tunnelling level following the dot's ground state at each point.** Rejected for the same
  reason as in ADR 0026. At a charge transition the electron tunnels at the leads' electrochemical
  potential, which no gate moves. With a moving level, a plunger would close its own dot's lead
  instead of opening it, and the compile tier and the runtime would disagree.
- **Only the barrier gate's own slope (B_L on Γ_L, and so on), like t's single-gate `Barrier`.**
  Rejected: the action responds to every gate, the plunger more strongly than the barrier in this
  potential. Dropping that would make the rates follow a gate the potential barely weights.
- **Recomputing the action at runtime.** Rejected by ADR 0001: a line cut through the basis per
  sample is far slower than 1 ms per point.
- **The rates as an exponential without a ceiling.** Rejected: opening a barrier would raise the rate
  past the attempt frequency.
- **A new release, 0.4.0, for this change.** Rejected by the owner: 0.3.0 is untagged and nobody pins
  it, so one breaking release carries the sensor dot and the dots' dynamics.
