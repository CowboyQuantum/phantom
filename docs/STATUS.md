# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** 1 (issue #16, awaiting sign-offs) and 2 (issue #17, in progress) together →
`docs/steps/01-charge-readout-core.md`, `docs/steps/02-layout-electrostatics.md`. Milestone: "Steps 1-2".
**Branch:** `step-2/electrical-health` (second Step 2 PR, part of #17). PRs #19–#23 merged.
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
- Second Step 2 PR (`step-2/electrical-health`): `gate-leakage-iv`, `ohmic-iv`, `turn-on-global`,
  `turn-on-single-gate-layout`, `gate-hysteresis`, `crosstalk-lines` rendering; `readout.transport` gains
  `channel_transmission_current_a`, `gate_leakage_current_a`, `ohmic_current_a`; `Gate.area_nm2`; 3 tests.
- Merged as PR #23, first Step 2 PR: `layout` (Gate, Layer, Stack, MaterialClass, Layout JSON, `two_dot_layout`), `electrostatics`
  (solid-angle closed form for any polygon, GateBasis, dots and barriers from the landscape, offline Schrödinger),
  `charge.derived.derive_double_dot`, ADR 0013, reference entry `davies1995-patterned-2deg` (proposed), four Step 2
  protocols rendering (potential map, channel cut, layout-derived stability map, pinch-off crosstalk), 8 new P0 tests.
- Step 1 sign-offs still pending. PRs #21 and #22 merged.

## Known gaps
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
Owner reviews the Step 2 pictures on the preview of the second Step 2 PR and merges; names the published layout for
the P1 row; signs off protocols from Steps 1 and 2 as reviewed. Then Step 3 (disorder and the fault models, issue to
be filed) opens, which needs the catalog entries (#1–#14) the owner is writing.

Previous next action, still valid for Step 1:
When CI on PR #19 is green: owner reviews the gallery artifact and the two reference entries; the session iterates
on protocols from the review, then writes sign-offs as instructed. In parallel, start Step 2 on
`step-2/layout-electrostatics` (issue #17): `layout` dataclasses, Davies–Larkin–Stern closed-form potential,
lever arms derived from the layout, `potential-map-layout` and `potential-cut-along-channel` protocols first.
