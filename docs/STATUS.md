# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** 8 (qdt cutover, #40), phantom side merged (PRs #42, #44); `v0.2.0` goes on `f72daac`. Steps 1–4
complete except the owner's sign-offs, catalog entries and P1 rows (and Step 4's lock-in dynamics on the master
equation); 5 (#32) and 6 (#34) built, their P1 rows open; 7 (#37) built except the solver decision (plan §8 #1).
Milestones: "Steps 1-2", "Steps 3-4", "Steps 5-6: valley, spin", "Steps 7-8: compile, cutover".
**Branch:** `issue-45/release-tagging` (ADR 0025). `main` is at `f72daac`, the merge of PR #44, CI green.
**Last session:** 2026-10-27 — Step 8's phantom side: the consumer surface, the contract test, the API snapshot,
the changelog and version 0.2.0 (ADR 0024), with qdt's six requests from its Analysis 0007 built in; the
sensor-asymmetry fix (#43, PR #44); who tags a release (ADR 0025).

## Done
- Steps 1–7 as before: `layout`, `electrostatics`, `disorder`, `charge`, `readout`, `valley`, `spin`, `compile`,
  `device`, `optimum` (ADRs 0001–0023); tooling; CI with the ubuntu reference and the Windows comparison.
- Step 8, phantom side (ADR 0024, `docs/CONTRACT.md`, `CHANGELOG.md`):
  - `phantom.api`, the one consumer import, above `compile` with `phantom.contract`; rule 2 amended in `CLAUDE.md`
    and enforced by `scripts/check_imports.py` (no layer or root module imports either); `model_versions()`.
  - `device.session`: `Device.at(v, v_sd)` → `State`; `Device.session(noise, dynamics, charge, noise_seed)` →
    `Session` with a forward clock: `sense(v, dwell_s, samples, v_sd)` (sample k at t₀ + (k+1)·dwell), `hold(v,
    duration_s)`, `wait(duration_s)`, `spin_readout(...)`; `Device.measure`. One sample's charge: a Markov jump
    trajectory or an independent Boltzmann draw (`sampled`, default) or ⟨N⟩ (`mean`). `noise_seed` gives a new
    realisation with the device's noise parameters; `None` replays the device's own.
  - `device.truth`: `Device.truth()` (charge lattice through the wiring, `SensorTruth`: peak, period, charging
    energy, lever arm, width, wing, sensitivity at a bias, as qdt's tune-up nodes name them) and
    `Device.knowable()` (material, gate names, layout, temperature, programmed point, sensor/barrier/reservoir gates).
  - `python -m phantom.contract [device.json ...]`: fifteen checks, exit 0 or 1; `tests/test_contract.py` (three
    draws), `tests/test_device_session.py` (closed forms: Boltzmann draws, trajectory switching against the
    two-state propagator, white floor under a dwell change, sensor truth against the sensor current, lattice
    through the wiring, PSB occupations). The contract holds on every compilable seed 0–11 (6 does not compile).
  - `scripts/api_surface.py` → `docs/api/surface.json`, `--check` in `check_all.sh`.
  - `RateEquation.rate_matrix` vectorised, bit-identical (checked on 1,200 matrices and four sweeps): a session with
    dynamics costs about 0.6–0.9 ms per sample.
  - Protocols `session-charge-trace`, `session-hold-relaxation` (unsigned, within budget).
- #43 (PR #44): the compiled sensor carries the wiring draw's coupling asymmetry (`device.compiled_sensor`, used by
  compile's set point and `Device.runtime_sensor`); `device` and `compile` MODEL_VERSION 2.
- 66 protocols render within budget, all unsigned; 136 P0 tests; eight reference entries, all `proposed`.

## In flight
- `issue-45/release-tagging` (#45): ADR 0025 (a session tags a release when the owner says so, on a green `main`
  with the changelog entry), the `CLAUDE.md` line, ADR 0024's status line, the documents that repeated the rule.

## Next action
- The owner tags `v0.2.0` on `f72daac` (2026-10-27: this session's environment refused the tag push with HTTP 403;
  its branch pushes work). Merge #45's pull request.
- The qdt side lives in CowboyQuantum/qdt: #141 (the sequencer interpreter moves into qdt-lab, needs no phantom
  release), then #140 and #142 against the tag (the adapter over `phantom.api.Session`, committed device files with
  SHA-256, `python -m phantom.contract devices/*.json` in qdt's CI). Coordinate there when the owner says so.
- Later, separate: spin driven by gate voltages over time (#41).
- Open decisions: plan §8 #1 (solver); Step 8's cadence of phantom releases against qdt pins (held by the owner).

## Known gaps
- Owner-dependent, never blocking: sign-offs for all 66 protocols; catalog entries #1–#14; confirmation (or
  replacement) of the eight proposed reference entries; the papers for the P1 rows of Steps 2–6.
- The compiled lead and interdot rates are constant over voltage at runtime (`State` says so); a session's charge
  dynamics use them everywhere.
- Compile moves the sensor gate onto its flank in one step and re-derives the dots, so a draw whose charge shifts
  with the sensor gate lands off the 0.9-width target (seeds 0 and 3: 0.87 and 0.79 widths).
- `spin_readout` takes a prepared triplet probability and does not read the session's charge or voltages (#41).
- A point evaluated alone and inside a batch agree within 1e-12 relative, not bit for bit (numpy's kernels); the
  contract states it.
- The clean textbook layout's interdot barrier is high: at B_M = 0.2 V the clean t is 0.03 µeV. Across seeds 0–11 at
  the gallery operating point the compiled t spans 0.003–900 µeV and seed 6's interface charges stop the sensor dot
  from forming; the placeholder charge density dominates the spread until catalog entry #1 sets it. The gallery
  device is seed 2 (t = 52 µeV, lead rates 1e9–1e11 Hz). The barrier gate's exponential scale is 7–10 mV per factor e,
  a property of the closed-form potential without screening.
- No screening until the FEM backend exists. Trapping hysteresis (`disorder.trapping`) is not part of the compiled
  device.
- Spin, valley, noise and disorder placeholders (`calibration: none`) as listed in the step files; the quasi-static
  detuning noise takes the whole 1/f band; Step 4's lock-in dynamics on the master equation still owed; amplifier
  noise absent from the SNR picture; the two-dot CI model cannot represent a plunger–barrier short.
- Cross-platform agreement (ADR 0014, 1e-9 relative): quantities that cancel digits (eigenvalue differences, finite
  differences, argmax ties on a flat maximum, interpolation through near-equal values) have each broken the Windows
  gallery comparison once; use closed forms, Hellmann–Feynman derivatives, grid points or secular equations instead.
- Pages source is `gh-pages` (root); the docs root and per-PR previews under `/pr-<n>/` are served.
