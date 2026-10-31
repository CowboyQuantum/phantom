# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** none in flight. Step 8's phantom side was released as `v0.2.0` (on `f72daac`). Steps 1–7
are built; their open items (sign-offs, catalog entries, P1 rows, the solver) are on #16, #17, #25, #28,
#32, #34 and #37. 0.3.0 is untagged and now carries #47, #53 and #54 (`CHANGELOG.md`).
**Branch:** `issue-54/sensor-charge-dynamics` (PR #56, against `main`). `main` is at `45e27a0`, the merge of
#55 (#53).
**Last session:** 2026-10-30. The owner asked for better charge dynamics for the sensor and for the
gates to modulate them. Two decisions were taken that day:
- **Back-action:** on, always with dynamics.
- **Release:** both changes fold into the untagged 0.3.0.

## Done this session
- **#53 ([ADR 0027](adr/0027-the-qubit-dots-rates-follow-the-gates-at-runtime.md)): the qubit dots'
  rates follow the gates.**
  - Each `CompiledBarrier` carries its attempt rate, a tunnelling level held fixed, and the slopes of
    ln Γ against every gate.
  - At runtime, `readout.rates.GateRates` and `saturating_rates_hz` (shared with the sensor's leads)
    give the rates, and `RateEquation.gate_rates` uses them.
  - `State`'s three rates are arrays. `DEVICE_FORMAT_VERSION` 3.
  - New protocol `session-closing-barrier`. With B_L 0.1 V lower, Γ_L is about 10⁴ Hz, a sweep at
    3 µs per point opens a 3.5 mV loop, and a parked session telegraphs.
  - `session-hold-relaxation` now overlays the total rate out of the held state. The interdot path
    had doubled the relaxation rate before as well.
  - Tests: `tests/test_dot_rates.py`, 9 P0.
- **#54 ([ADR 0028](adr/0028-the-sensor-dots-charge-is-a-process-in-time-with-back-action.md)): the
  sensor dot's charge is a process in time.**
  - With dynamics, the sensor's charge is a birth–death chain over its levels with broadened
    occupations, so its stationary current is exactly ADR 0026's.
  - A sample reads the electrons that left to the drain during its dwell: from exact trajectories at
    ≤ 200 jumps per dwell, otherwise a Gaussian with the full-counting-statistics variance.
  - Back-action is −β_i·(M − M̄), against the sensor's stationary mean.
  - New code: `readout.sensor_dynamics` and `Measurement.sensor_charge`.
  - New protocol `session-slow-sensor`.
  - Tests: `tests/test_sensor_dynamics.py`, 10 P0.
- **Totals:** 174 tests; 72 protocols within budget, all unsigned; contract OK; ledger re-recorded.

## In flight
- PR #56 (#54) waits for the owner's review of its preview at `/pr-56/`. Its first Windows gallery
  comparison failed by 1.5·10⁻⁹, from an exponential at ΓT ≈ 10⁸. The fix in `9a58df2` replaces it with
  the relaxed closed form beyond 50 relaxation times (ADR 0028).

## Next action
The owner reviews and merges #56. They also decide whether to tag `v0.3.0`
(ADR 0025) and when qdt moves its pin (CowboyQuantum/qdt #140, #142, and #169 for the run option).

## Open issues
#48 rates' unscreened sensitivity (now also the dots': P1 moves Γ_L more than B_L); #49 per-device
sensor-barrier tune; #50 gate 1/f placeholder; #51 cp1252 console; #41 spin driven by voltages.
Open blockers: none. Open decisions: plan §8 #1 (solver); release cadence against qdt's pins.

## Known gaps
- **Owner-dependent, never blocking:** sign-offs for all 72 protocols; catalog entries #1–#14;
  confirming the ten `proposed` reference entries.
- **Rates:** linearised exponents within about 20 % of the compile tier over ±10 mV, drifting beyond
  that. No screening until the FEM backend exists.
- **Sensor process:** the split is first order within a dwell. The noise is sequential with broadened
  occupations. There is no joint energy between the dots and the sensor (ADR 0028).
- **Placeholders** (`calibration: none`): as listed in the step files and ADRs 0019–0028.
- **Cross-platform (ADR 0014):** avoid digit-cancelling quantities; use closed forms, grid points and
  bisection.
- **Pages:** the source is `gh-pages` (root); per-PR previews are under `/pr-<n>/`.
