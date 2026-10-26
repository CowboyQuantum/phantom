# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** 1–4 complete except the owner's sign-offs, catalog entries and P1 rows (and Step 4's lock-in
dynamics on the master equation); 5 (#32) and 6 (#34) built, their P1 rows open; 7 (#37, compile and optimum) in
progress → `docs/steps/07-compile-optimum.md`. Milestones: "Steps 1-2", "Steps 3-4", "Steps 5-6: valley, spin",
"Steps 7-8: compile, cutover".
**Branch:** `step-7/compile-device` (first Step 7 PR, part of #37). PRs #19–#36 merged.
**Last session:** 2026-10-26 — Steps 4, 5 and 6 closed (#31, #33, #35, #36); Step 7 opened: `compile`, the new
`device` layer, `optimum`, ADR 0023, three protocols, ten P0 tests.

## Done
- `compile` / `device` / `optimum` (Step 7): `compile_device(layout, CompileParams, seed)` runs the scenario draw,
  the basis with the draw applied (closed form; `fem` hook behind `[fem]`), the wiring map, the dots and CI model with
  the sensor on its flank, the WKB-action barriers and rates, the interdot t and its gate scale, the valley and g draws,
  and writes a `phantom.device.Device` (parameters only, canonical JSON, format and model versions; basis recomputed
  on demand). `Device.runtime_*` build the CI model, sensor, barrier, Hubbard model, rate equation, sensed double dot
  with the noise draw, and spin draw. `optimum.optimum(device)`: sensor flank (coarse scan + golden section),
  virtual-gate matrix, barrier voltage for a target t, PSB read point and Elzerman peak time, exchange sweet spot in
  the operating cell with T₂* and SWAP fidelity. `electrostatics.barrier_action` (WKB action along the cut). The
  layer order gains `device` between `spin` and `optimum` (rule 2, checker, ADR 0004 note). ADR 0023. Protocols
  `compiled-device-roundtrip`, `virtual-gates-map`, `optimum-sensor-flank` (unsigned, within budget).
- `spin` (Step 6): exchange from the Hubbard sector, Lindblad dynamics, PSB, blockade transport, funnel, Elzerman,
  noise spectroscopy, gate fidelity; ADR 0021, 0022; fourteen protocols; three proposed entries.
- `valley` (Step 5): Rice draw per dot, three signatures; ADR 0020; four protocols; two proposed entries.
- Steps 1–4: `layout`, `electrostatics`, `disorder`, `charge`, `readout`; 43 protocols; tooling; CI with the ubuntu
  reference and the Windows comparison within 1e-9; per-PR previews; `docs/MODELS.md` map, graph and ledger.

## In flight
- First Step 7 PR (this branch): everything under "Done" for Step 7; `docs/MODELS.md` (device, compile and optimum
  subgraphs; graph and ledger refreshed); roster section I; step file; indexes; `measured.json`. 114 P0 tests.

## Next action
- Second Step 7 PR: move the Step 2–6 protocols that build a device by hand onto `compile_device` (one compiled device
  per protocol family, the Step 4 rectangular WKB replaced by the action), an `optimum` protocol per remaining layer
  (coupling target, read point, sweet spot), and the three-dot layout question. Then Step 8 (`qdt` cutover): the
  adapter in qdt pins a phantom release and maps its objects onto `Device` and `optimum`.
- Open decisions: plan §8 #1 (solver) stays open; the `fem` backend raises until it is made.

## Known gaps
- Owner-dependent, never blocking: sign-offs for all 60 protocols; catalog entries #1–#14; confirmation (or
  replacement) of the eight proposed reference entries; the papers for the P1 rows of Steps 2–6.
- The clean textbook layout's interdot barrier is high: at B_M = 0.2 V the clean t is 0.03 µeV while draws with
  interface charges give 4–370 µeV (seeds 1–3); the gallery device is seed 2 (t = 52 µeV, lead rates 1e9–1e11 Hz).
  The barrier gate's exponential scale comes out 7–10 mV per factor e (steeper than the 10–100 mV the `Barrier`
  docstring quotes), a property of the closed-form potential without screening.
- Two WKB estimates coexist: the rectangular one in Step 4's protocols and `readout.rates`, the action in compiled
  devices, until the Step 4 protocols move to compiled devices. No screening until the FEM backend exists.
- Spin, valley, noise and disorder placeholders (`calibration: none`) as listed in the step files; the quasi-static
  detuning noise takes the whole 1/f band; Step 4's lock-in dynamics on the master equation still owed; amplifier
  noise absent from the SNR picture; the two-dot CI model cannot represent a plunger–barrier short.
- Pages source is `gh-pages` (root); the docs root and per-PR previews under `/pr-<n>/` are served.
