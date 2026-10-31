# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** none in flight. Released: `v0.2.0` (on `f72daac`, Step 8's phantom side) and `v0.3.0`
(on `6920f19`, carrying #47, #53 and #54). Steps 1–7 are built; their open items (sign-offs, catalog
entries, P1 rows, the solver) are on #16, #17, #25, #28, #32, #34 and #37.
**Branch:** `issue-54/release-v0.3.0` (docs only: the changelog heading and this file). `main` is at
`6920f19`, the merge of #56 (#54).
**Last session:** 2026-10-30 to 2026-10-31. The owner asked for better charge dynamics for the sensor and
for the gates to modulate them. Three decisions were taken:
- **Back-action:** on, always with dynamics.
- **Release:** both changes fold into 0.3.0.
- **Tag:** on 2026-10-31 the owner asked for `v0.3.0` and a GitHub Release made as for `v0.2.0`.

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
  - Tests: `tests/test_dot_rates.py`, 9 P0. Merged as #55.
- **#54 ([ADR 0028](adr/0028-the-sensor-dots-charge-is-a-process-in-time-with-back-action.md)): the
  sensor dot's charge is a process in time.**
  - A birth–death chain over the sensor's levels with broadened occupations; a sample reads the
    electrons counted to the drain (exact trajectories at ≤ 200 jumps, otherwise Gaussian with the
    full-counting-statistics variance; the relaxed closed form beyond 50 relaxation times).
  - Back-action −β_i·(M − M̄), against the sensor's stationary mean.
  - New: `readout.sensor_dynamics`, `Measurement.sensor_charge`, protocol `session-slow-sensor`.
  - Tests: `tests/test_sensor_dynamics.py`, 10 P0. Merged as #56.
- **Release `v0.3.0`:** annotated tag on `6920f19` after CI passed there (Windows gallery comparison
  included); GitHub Release "phantom 0.3.0" with the changelog entry as its notes, as for `v0.2.0`.
  Noted on #53 and #54.
- **Totals:** 174 tests; 72 protocols within budget, all unsigned; contract OK.

## In flight
- This branch's PR: the changelog heading names the tag, and this file.

## Next action
qdt moves its pin to `v0.3.0` when the owner decides (CowboyQuantum/qdt #140, #142, and #169 for the run
option). Its adapter must handle `State`'s rates as arrays, the 12-gate layout and device format 3
(`CHANGELOG.md`, Breaking).

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
- **Tagging:** this environment pushed `v0.3.0`; it refused `v0.2.0` (ADR 0025). Whether a session can
  tag depends on its credentials.
