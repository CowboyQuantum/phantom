# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** none in flight. Step 8's phantom side is released as `v0.2.0` (on `f72daac`); Steps 1–7 built
(their open items are on #16, #17, #25, #28, #32, #34, #37: sign-offs, catalog entries, P1 rows, the solver).
**Branch:** `issue-47/sensor-barriers-bath` (#47) from `main` at `84f9efa`; its pull request is open for review.
**Last session:** 2026-10-29 — the owner asked for the charge sensor dot to have its own barriers and a bath that
can be tuned: built as #47, [ADR 0026](adr/0026-the-sensor-is-a-dot-with-its-own-barriers-and-bath.md), 0.3.0
(untagged, breaking: `CHANGELOG.md`).

## Done this session
- `two_dot_layout` gains `SB_L`, `SB_R`, `SR_L`, `SR_R` after the first eight gates (`src/phantom/layout/geometry.py`).
- A compiled sensor with leads is a dot like the qubit dots: E_c and β from the three-dot capacitance matrix
  (`charge/derived.py`), rates through its barriers with gate slopes and an attempt-rate ceiling, baths with a
  turn-on (`readout/sensor.py`), a lifetime-broadened current (`readout/broadening.py`), compiled in
  `compile/device.py`; `Device` format 2; `State`, `Knowable`, `SensorTruth` extended; contract and API snapshot
  updated. The typed sensor is unchanged for layouts without sensor reservoirs.
- Reference programmed point adds `SB_L 0.4, SB_R 0.3, SR_L 0.7, SR_R 0.7` (set for seed 2).
- `tests/test_sensor_leads.py` (18 P0); 155 tests; four protocols (`gallery/protocols/step7c.py`), 70 rendering
  within budget, all unsigned; two `proposed` reference entries; `docs/MODELS.md`, `measured.json` re-recorded.

## In flight
- The pull request for #47: CI's ubuntu reference run and the Windows comparison (first cross-platform check of
  the sensor slopes' central differences and the wing bisection, ADR 0014), and the preview at `/pr-<n>/`.

## Open issues from this session
- #48 sensor rates' sensitivity and the bath gate's cross-talk (no screening; plan §8 #1).
- #49 no per-device sensor-barrier tune; 4 of 12 seeds form no device at the common point.
- #50 the gate 1/f noise placeholder is about ten times typical and now visible on the sensor.
- #51 `python -m phantom.contract` cannot print on a cp1252 Windows console (bug).
- Open blockers: none. Open decisions: plan §8 #1 (solver); release cadence against qdt's pins; #41 later.

## Next action
The owner reviews #47's pull request (the changed and the four new protocols in its preview, ADR 0026), merges it;
then decides whether to tag `v0.3.0` (ADR 0025) and when qdt (CowboyQuantum/qdt #140, #142) moves its pin.

## Known gaps
- Owner-dependent, never blocking: sign-offs for all 70 protocols; catalog entries #1–#14; confirming the ten
  `proposed` reference entries.
- The dots' lead and interdot rates are constant at runtime; no screening until the FEM backend; placeholders
  (`calibration: none`) as listed in the step files and ADRs 0019–0026.
- Cross-platform agreement (ADR 0014): avoid digit-cancelling quantities; use closed forms, grid points, bisection.
- Pages source is `gh-pages` (root); per-PR previews under `/pr-<n>/`.
