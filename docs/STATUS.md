# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** 8 (qdt cutover) not started; its phantom-side design is written under "Next action". Steps 1–4
complete except the owner's sign-offs, catalog entries and P1 rows (and Step 4's lock-in dynamics on the master
equation); 5 (#32) and 6 (#34) built, their P1 rows open; 7 (#37) built except the solver decision (plan §8 #1).
Milestones: "Steps 1-2", "Steps 3-4", "Steps 5-6: valley, spin", "Steps 7-8: compile, cutover" (#4 on GitHub).
**Branch:** `step-8/consumer-contract`, which holds only this STATUS commit; `main` is at the merge of PR #39.
PRs #19–#39 merged; no pull request is open.
**Last session:** 2026-10-26 to 2026-10-27 — Steps 4 to 7 built and merged (PRs #31, #33, #35, #36, #38, #39). The
owner stopped the session as it was filing the Step 8 issue and writing `src/phantom/device/session.py`; neither was
done. Start Step 8 when the owner asks for it.

## Done
- Steps 1–4: `layout`, `electrostatics`, `disorder`, `charge`, `readout`; tooling (`check_references`, `gallery.py`
  with budgets, sign-off hashes and `--record`, `model_graph.py --check`); CI with the ubuntu reference and the Windows
  comparison within 1e-9 (ADR 0010, 0014); per-PR gallery previews on Pages (ADR 0012); `docs/MODELS.md`.
- Step 5 `valley`: Rice draw per dot with the alloy-disorder scalings; magnetospectroscopy kink, T₁ hot spot,
  pulsed-gate spectroscopy (ADR 0020).
- Step 6 `spin`: exchange from the Hubbard singlet sector (MODEL_VERSION 2: secular equation by bisection, full
  relative precision at tiny coupling), Lindblad dynamics with quasi-static noise by quadrature, PSB window from the
  valley draw, blockade transport, funnel, Elzerman shots, CPMG noise spectroscopy, SWAP fidelity (ADR 0021, 0022).
- Step 7 `compile` / `device` / `optimum` (ADR 0023): `compile.compile_device(layout, CompileParams, seed)` →
  `phantom.device.Device` (parameters only, canonical JSON with format and model versions, the basis recomputed on
  demand, `runtime_*` builders for every layer's objects); the `device` layer sits between `spin` and `optimum`
  (rule 2, `scripts/check_imports.py`); closed-form backend, `fem` hook that raises naming plan §8 #1;
  `optimum.optimum(device)`: sensor flank (coarse and fine grid scan, no interpolation, for cross-platform
  agreement), virtual-gate matrix (raw-scaled, cross-talk and sensor compensated), barrier voltage for a target t,
  PSB read point and Elzerman read time, exchange sweet spot with T₂* and SWAP fidelity. All tunnel rates come from
  the WKB action along the barrier (`electrostatics.barrier_action`, `readout.wkb_rate_from_action_hz`); the
  rectangular estimate remains only as the bound.
- 64 protocols render within budget, all unsigned; 117 P0 tests; eight reference entries, all `proposed`.

## In flight
- Nothing. No branch other than this one carries unmerged work.

## Next action
Step 8, phantom side (`docs/steps/08-qdt-cutover.md`). The design the last session had prepared:
1. **File the Step 8 issue** under milestone "Steps 7-8: compile, cutover" (title "Step 8 — qdt cutover (phantom
   side)"; goal, qdt side and exit as below) and continue on this branch.
2. **`src/phantom/device/session.py`** (plan §4's runtime surface):
   - `State` from `Device.at(v_programmed)`: carried voltages after the wiring faults, thermal and ground-state
     occupation, sensor potential and noise-free current, interdot t at the barrier gate, exchange energy, and the
     compiled rates (left lead, right lead, interdot). Vectorised over leading axes.
   - `Session` from `Device.session(noise=True, dynamics=True)`: owns fresh noise processes from the draw
     (`runtime_sensed(1 ms)`), the rate equation, a clock that only moves forward, and the charge distribution.
     `sense(v_programmed, dwell_s)` samples at `t0 + k·dwell`, adds the gate noise to the carried voltages, propagates
     the charge with `RateEquation.sweep` carrying `(p, states)` between calls (equilibrium occupation when
     `dynamics=False`), applies sensor and readout noise (white readout noise scaled by √(1 ms / dwell)), and returns
     `Measurement(times_s, current_a, n_present)`. `wait(duration_s)` relaxes the charge at the last point.
     `spin_readout(p_triplet, detuning_ev, shots, t_integration_s, t1_read_s)` samples `SpinOutcome(prepared_triplet,
     in_02)` from `spin.psb.psb_occupation` with the receiving dot's valley splitting, relaxation during integration,
     and an rng seeded `[device.seed, 8]`. A fresh session replays the same noise, so identical calls on two fresh
     sessions give identical results.
   - On `Device`: `at`, `session`, `measure(v, dwell_s)` (a fresh session's `sense`), and `truth` (the full record).
3. **`src/phantom/api.py`**: the one consumer import surface (`compile_device`, `CompileParams`, `Device`, `State`,
   `Session`, `Measurement`, `SpinOutcome`, `optimum`, `Optimum`, `Layout`, `Gate`, `Layer`, `Stack`,
   `two_dot_layout`, `ValleyParams`, `SpinParams`, `__version__`). Extend `scripts/check_imports.py` so no layer
   imports `phantom.api` or `phantom.contract`.
4. **`src/phantom/contract.py`** with `run()` and `python -m phantom.contract` (exit status 0 or 1), which qdt's
   adapter also runs: compile is deterministic in the seed; the JSON round trip is text-identical; `at` has the
   documented shapes and finite values; two fresh sessions give identical measurements; time cannot go backwards;
   with `noise=False` and a long dwell, `sense` matches `at`'s equilibrium current; `optimum` returns all five parts;
   `truth` carries the seed, scenario and dots; provenance names every model version. `tests/test_contract.py` runs it.
5. **`scripts/api_surface.py`** writing `docs/api/surface.json` (every public name, its signature, dataclass fields
   and public methods) with `--check` in `scripts/check_all.sh`, so an API change shows up as a diff that needs a
   changelog entry.
6. **`CHANGELOG.md`**, version 0.2.0 in `pyproject.toml` and `phantom.__version__`, ADR 0024 (the contract, semantic
   versioning for 0.x with breaking changes listed, the surface snapshot, release cadence left to the owner),
   `docs/CONTRACT.md` in the mkdocs nav, the step file, `docs/MODELS.md` and the indexes.
7. **The owner tags the release** after the PR merges. The qdt side (`qdt-sim` mapped onto `phantom.api.Device`,
   pinned to the tag, its physics deleted, its fixtures regenerated) lives in the public repository
   `CowboyQuantum/qdt`; coordinate it there when the owner says so.
- Open decisions: plan §8 #1 (solver); Step 8's cadence of phantom releases against qdt pins.

## Known gaps
- Owner-dependent, never blocking: sign-offs for all 64 protocols; catalog entries #1–#14; confirmation (or
  replacement) of the eight proposed reference entries; the papers for the P1 rows of Steps 2–6.
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
