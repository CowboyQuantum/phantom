# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** 1–4 complete except the owner's sign-offs, catalog entries and P1 rows (and Step 4's lock-in
dynamics on the master equation); 5 (#32) and 6 (#34) built, their P1 rows open; 7 (#37, compile and optimum) built
except the solver decision → `docs/steps/07-compile-optimum.md`. Milestones: "Steps 1-2", "Steps 3-4", "Steps 5-6:
valley, spin", "Steps 7-8: compile, cutover".
**Branch:** `step-7/optimum-protocols` (second Step 7 PR, part of #37). PRs #19–#38 merged.
**Last session:** 2026-10-27 — Step 7 first PR merged (#38: `compile`, `device`, `optimum`); second PR: the Step 4
dynamics protocols on the WKB action, optimum protocols for coupling, read point and sweet spot, and a population of
twelve compiled devices.

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
- Second Step 7 PR (this branch): `readout.wkb_rate_from_action_hz`; the gallery's Step 4 dynamics helper uses the
  action integral (operating points of `stability-vs-sweep-rate` and `stability-direction` moved to B_L = B_R = −0.02 V
  to keep a 16 kHz lead rate; `plunger-vs-barrier-map` sweeps B_L from −0.25 V); `compile.interdot_tunnel_coupling_ev`
  public; render option `markers`; protocols `optimum-coupling-target`, `optimum-readout-point`, `optimum-sweet-spot`,
  `compiled-device-draws`; two P0 tests; ADR 0017 status note and ADR 0023 consequence updated. 116 P0 tests.

## Next action
- Step 8 (`qdt` cutover), phantom side: a changelog and the consumer contract (plan §4) documented as the public API,
  a contract test over that surface that qdt's adapter can also run, and a release candidate for the owner to tag.
  The qdt side (the `qdt-sim` adapter over `phantom.device.Device`, its fixtures regenerated) lives in the qdt
  repository and is coordinated there. File the Step 8 issue under "Steps 7-8: compile, cutover".
- Open decisions: plan §8 #1 (solver) stays open; the `fem` backend raises until it is made. Step 8: the cadence of
  phantom releases against qdt pins.

## Known gaps
- Owner-dependent, never blocking: sign-offs for all 64 protocols; catalog entries #1–#14; confirmation (or
  replacement) of the eight proposed reference entries; the papers for the P1 rows of Steps 2–6.
- The clean textbook layout's interdot barrier is high: at B_M = 0.2 V the clean t is 0.03 µeV while draws with
  interface charges give 4–370 µeV (seeds 1–3); the gallery device is seed 2 (t = 52 µeV, lead rates 1e9–1e11 Hz).
  The barrier gate's exponential scale comes out 7–10 mV per factor e (steeper than the 10–100 mV the `Barrier`
  docstring quotes), a property of the closed-form potential without screening.
- The rectangular WKB estimate remains only as the bound in `readout.rates`; every rate the gallery and the compiled
  device use comes from the action integral. No screening until the FEM backend exists.
- Across seeds 0–11 at the gallery operating point the compiled interdot coupling spans 0.003–900 µeV and seed 6's
  interface charges stop the sensor dot from forming: the placeholder charge density near a high clean barrier
  dominates the spread until catalog entry #1 sets it.
- Spin, valley, noise and disorder placeholders (`calibration: none`) as listed in the step files; the quasi-static
  detuning noise takes the whole 1/f band; Step 4's lock-in dynamics on the master equation still owed; amplifier
  noise absent from the SNR picture; the two-dot CI model cannot represent a plunger–barrier short.
- Pages source is `gh-pages` (root); the docs root and per-PR previews under `/pr-<n>/` are served.
