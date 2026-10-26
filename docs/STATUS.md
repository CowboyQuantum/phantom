# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** 1 (#16), 2 (#17), 3 (#25) complete except the owner's sign-offs, catalog entries and P1 rows;
4 (#28, readout physics) in progress → `docs/steps/04-readout.md`. Milestones: "Steps 1-2", "Steps 3-4".
**Branch:** `step-4/noise-draws-and-dynamics` (third Step 4 PR, part of #28). PRs #19–#30 merged.
**Last session:** 2026-10-25 — Step 1 core and validation tooling built: CI and Hubbard charge models, SET sensor,
noise, lock-in channel, qdt-sim reproduction fixture; `check_references.py`; `gallery.py` with budgets and
sign-off hashes; eleven protocols rendering; CI gallery jobs (reference + Windows within tolerance), Pages job; ADR 0010.

## Done
- `charge`: `ConstantInteractionParams` (qdt-sim's model migrated, plus exact lever arms, addition energies,
  triple points), `Dot`, `Barrier`, `HubbardParams` (two-site spinful, t → 0 reproduces CI exactly; P0 tested).
- `readout`: `SensorSETParams`, noise processes with `trace()` bit-identical to qdt-sim's `sample()` sequences,
  `SensedDoubleDot` (sweeps with raster times and noise), `LockinParams` / `lockin_sweep_a_per_v`
  (`sensor_lockin` channel), `transport` (phenomenological turn-on and pinch-off, `calibration: none`).
- `tests/fixtures/qdt_sim_defaults.npz` (synthetic, from qdt commit ae6dd2c): maps, diamonds and noise sequences
  reproduced to 1e-12 / bit for bit. 38 P0 tests green.
- `scripts/check_references.py` in `check_all.sh`; `tests/conftest.py` skips `literature` tests whose entry is
  not `confirmed`; two `proposed` entries (van der Wiel 2003 RMP; Connors 2019 PRB) awaiting the owner's check.
- `gallery/` package + `scripts/gallery.py`: protocol YAML loader, counters and timing, canonical array hash,
  sign-off check, budget ledger, PNG + `index.html`. Fourteen Step 1 protocols in `docs/gallery/protocols/`, all
  within budget locally in `--strict` mode; all `unsigned`.
- Owner-driven revisions: `plunger-vs-barrier-map` in the lead-barrier configuration (dot 1 alone has a reservoir)
  with `readout.loading` (rate-limited loading, latching); four chopped-gate maps `stability-lockin-p1/-p2/-common/
  -detuning` through `lockin_sweep_dynamic_a_per_v` with a `LeadRateLandscape` (latching below the chop
  frequency, lifetime broadening far above it; three limit tests). `docs/MODELS.md` measurement-flow diagram and
  generated module graph and the runtime ledger (per-model micro-benchmarks and per-protocol costs recorded in `docs/gallery/measured.json` by `gallery.py --record`; `scripts/model_graph.py --check` in `check_all`), ADR 0011.
- CI: `gallery` (ubuntu, the reference, ADR 0010) and `gallery-windows` (compares its arrays to the reference within
  1e-12) upload `gallery-<os>` artifacts; `pages` job builds mkdocs + rendered gallery from `main`.

## In flight
- Third Step 4 PR: `disorder.noise_draw` (per-device 1/f amplitude, Poisson switchers, white floor; in the scenario
  draw), `SensedDoubleDot.from_noise_draw`, ADR 0019, protocols `switcher-histogram`, `drift-repeated-cuts`;
  `plunger-vs-barrier-map` re-rendered on the master equation with WKB rates per row; `LeadBarrier` retired;
  four P0 tests.
- Merged as PR #30, second Step 4 PR: `readout.rf` (lumped tank + sensor resistance: reflection, resonance, matching),
  `readout.single_shot` (consecutive integrated shots on the noise processes, SNR, fidelity; white noise rescaled by
  1/√dt), ADR 0018, protocols `rf-resonance`, `rf-sensor-flank`, `single-shot-histogram`, `snr-vs-integration-time`,
  five P0 tests.
- Merged as PR #29, first Step 4 PR: `readout.rates` (WKB), `readout.master_equation` (RateEquation with detailed balance, sweep
  propagation), `electrostatics.barrier_width_nm` at the lead Fermi level, reservoir gates in `two_dot_layout`,
  CI MODEL_VERSION 2 (duplicate clipped states no longer over-weighted; qdt-sim fixture compared away from the
  empty corner), ADR 0017, protocols `stability-vs-sweep-rate`, `stability-direction`, `lead-rate-time-trace`,
  six P0 tests (detailed balance, latching threshold, direction dependence, metastability without an interdot rate).
- Merged as PR #27, second Step 3 PR: `disorder.wiring` (open/shorted gates as a gate matrix, leak network, sensor-coupling factors),
  `disorder.trapping` (history-dependent trapped fraction, hysteresis loop), both in the scenario draw; ADR 0016;
  protocols `fault-signatures-stability`, `gate-leakage-pair`, `gate-hysteresis-trapping`; six P0 tests.
- Merged as PR #26, first Step 3 PR: `disorder` layer (interface charges as image-charge energy offsets, misalignment and width error as
  polygon transforms, seeded serialisable scenario draw, severity 0 bit-identical to clean), `GateBasis.offset_ev`,
  `phantom.constants` below every layer, ADR 0015, three protocols (`potential-map-disorder`,
  `stability-coarse-2dot-disorder`, `unintentional-dots`), six P0 tests; catalog draft entry gains `calibration: none`.
- Merged as PR #24, second Step 2 PR: `gate-leakage-iv`, `ohmic-iv`, `turn-on-global`,
  `turn-on-single-gate-layout`, `gate-hysteresis`, `crosstalk-lines` rendering; `readout.transport` gains
  `channel_transmission_current_a`, `gate_leakage_current_a`, `ohmic_current_a`; `Gate.area_nm2`; 3 tests.
- Merged as PR #23, first Step 2 PR: `layout` (Gate, Layer, Stack, MaterialClass, Layout JSON, `two_dot_layout`), `electrostatics`
  (solid-angle closed form for any polygon, GateBasis, dots and barriers from the landscape, offline Schrödinger),
  `charge.derived.derive_double_dot`, ADR 0013, reference entry `davies1995-patterned-2deg` (proposed), four Step 2
  protocols rendering (potential map, channel cut, layout-derived stability map, pinch-off crosstalk), 8 new P0 tests.
- Step 1 sign-offs still pending. PRs #21 and #22 merged.

## Known gaps
- All Step 2/3 pictures changed with the reservoir gates and the CI fix (re-rendered, unsigned). Operating points
  carry `R_L: 0.7, R_R: 0.7`; the lead rate spans 1e2–1e10 Hz as the outer barrier goes 0.0–0.2 V.
- Step 4 still owes the lock-in dynamics on the master equation (`readout.loading.LeadRateLandscape` and
  `latched_occupation` remain for it; a master-equation step per half-cycle costs minutes per map) and the two P1
  rows (a charge-noise spectrum, a latching signature) once their entries exist. Amplifier noise is absent from the
  SNR picture (device-side bound). Noise draws are placeholders until entries #1 and #12 are written.
- A plunger shorted to a neighbouring barrier opens that barrier and dissolves the dot into the reservoir (or merges the
  dots); the two-dot CI model cannot represent it, so the fault protocol shows a plunger–sensor-gate short instead.
- Step 3 generators run on the catalog draft's placeholder densities (`calibration: none`); the owner's catalog entries
  (#1–#14) set the real distributions. Strain and g-factor draws wait for `spin`; switcher and drift rates tie to the
  catalog in Step 4. No screening in the charge potentials (upper bound) until Step 7.
- ADR 0014: cross-platform tolerance raised from 1e-12 to 1e-9 after the layout-derived stability maps differed by
  2e-11 between platforms (derived chains amplify last-bit differences at sharp steps). Typed chains still agree to 3e-14.
- Step 2 still owes the P1 reproduction against a published layout (its reference entry does not exist yet; the
  owner names the paper) and the owner's sign-offs of its protocols. `pinch-off-barrier` from the layout is
  covered by `turn-on-single-gate-layout` (same sweep seen from the barrier side); say so if a separate picture is wanted.
- Pages source is `gh-pages` (root); the docs root and per-PR previews are served.
- Pages: ADR 0012 moves publishing to the `gh-pages` branch with per-PR previews at `/pr-<n>/`; the repository's
  Pages source must be set to that branch (root) once it exists, by the owner if the session is not permitted.
- No protocol is signed off. The owner reviews the artifact, we iterate in chat, then a session writes
  `docs/gallery/signoff/<id>.yaml` at the owner's instruction (ADR 0008).
- The two reference entries are `proposed`; no P1 test exists yet. The owner confirms or replaces them.
- `turn-on-single-gate` and `pinch-off-barrier` are phenomenological shapes until Step 2.
- Cross-platform bit identity failed on the first CI run of PR #19 (ten of eleven protocols; even pure `exp`/`tanh`
  curves differ in the last bit between the Ubuntu and Windows numpy builds). ADR 0010 now makes `ubuntu-latest`
  the reference platform for sign-off hashes and requires other platforms to agree within 1e-12 relative; the
  `gallery-windows` job compares to the reference artifact. Rule 3 in CLAUDE.md amended accordingly.
- `docs/catalog/entries/interface-fixed-charge.yaml` still has free-text sources and no `calibration` field
  (warned, not failed, because it is `draft`).
- Step 1 exit criteria still open: sign-offs; the benchmark is met locally (0.02 ms/point on
  `stability-coarse-2dot`, budget 0.1) but not yet confirmed on CI.

## Next action
Owner reviews the three Step 3 pictures on the PR preview and merges; writes or confirms catalog entries so the
generators get real densities; names the published layout for Step 2's P1 row and a MOS/SiGe unintentional-dot
comparison for Step 3's P1 row; signs off protocols from Steps 1–3. Next session: open Step 5 (valley: valley splitting draws per material and dot size, magnetospectroscopy and
pulsed-gate spectroscopy protocols); file its issue and milestone. The lock-in dynamics on the master equation
wait for a cheaper propagation (a fixed-rate analytic half-cycle update) and are a Step 4 follow-up.

Previous next action, still valid for Step 1:
When CI on PR #19 is green: owner reviews the gallery artifact and the two reference entries; the session iterates
on protocols from the review, then writes sign-offs as instructed. In parallel, start Step 2 on
`step-2/layout-electrostatics` (issue #17): `layout` dataclasses, Davies–Larkin–Stern closed-form potential,
lever arms derived from the layout, `potential-map-layout` and `potential-cut-along-channel` protocols first.
