# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** 1 — Charge and readout core → `docs/steps/01-charge-readout-core.md`, issue #16. Step 2 (issue #17)
runs together with it. Milestone: "Steps 1-2: charge, readout, layout, electrostatics".
**Branch:** `step-1/charge-readout-core` (PR #19 open, part of #16).
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
  generated module graph (`scripts/model_graph.py --check` in `check_all`), ADR 0011.
- CI: `gallery` (ubuntu, the reference, ADR 0010) and `gallery-windows` (compares its arrays to the reference within
  1e-12) upload `gallery-<os>` artifacts; `pages` job builds mkdocs + rendered gallery from `main`.

## In flight
- PR #19 (part of #16, does not close it): awaiting the owner's review of the fourteen rendered protocols (gallery
  report published as a private artifact page from this machine; CI artifact `gallery-ubuntu-latest`).

## Known gaps
- **GitHub Pages is not enabled**; the session could not enable it (needs the owner: Settings → Pages → Source:
  GitHub Actions). The `pages` job has `continue-on-error: true` until then; remove it after enabling.
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
When CI on PR #19 is green: owner reviews the gallery artifact and the two reference entries; the session iterates
on protocols from the review, then writes sign-offs as instructed. In parallel, start Step 2 on
`step-2/layout-electrostatics` (issue #17): `layout` dataclasses, Davies–Larkin–Stern closed-form potential,
lever arms derived from the layout, `potential-map-layout` and `potential-cut-along-channel` protocols first.
