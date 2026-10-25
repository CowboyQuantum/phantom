# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** 1 — Charge and readout core → `docs/steps/01-charge-readout-core.md`, issue #16. Step 2 (issue #17)
runs together with it. Milestone: "Steps 1-2: charge, readout, layout, electrostatics".
**Branch:** `step-1/charge-readout-core` (PR #19 open, part of #16).
**Last session:** 2026-10-25 — Step 1 core and validation tooling built: CI and Hubbard charge models, SET sensor,
noise, lock-in channel, qdt-sim reproduction fixture; `check_references.py`; `gallery.py` with budgets and
sign-off hashes; eleven protocols rendering; CI gallery job, cross-platform hash comparison, Pages job.

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
  sign-off check, budget ledger, PNG + `index.html`. Eleven Step 1 protocols in `docs/gallery/protocols/`, all
  within budget locally in `--strict` mode; all `unsigned`.
- CI: `gallery` job on both platforms uploads `gallery-<os>` artifacts; `gallery-hashes-agree` fails if the array
  hashes differ between platforms; `pages` job builds mkdocs + rendered gallery from `main`.

## In flight
- PR #19 (part of #16, does not close it): awaiting CI (first cross-platform hash comparison) and the owner's
  review of the eleven rendered protocols in the `gallery-*` artifact.

## Known gaps
- **GitHub Pages is not enabled**; the session could not enable it (needs the owner: Settings → Pages → Source:
  GitHub Actions). The `pages` job has `continue-on-error: true` until then; remove it after enabling.
- No protocol is signed off. The owner reviews the artifact, we iterate in chat, then a session writes
  `docs/gallery/signoff/<id>.yaml` at the owner's instruction (ADR 0008).
- The two reference entries are `proposed`; no P1 test exists yet. The owner confirms or replaces them.
- `turn-on-single-gate` and `pinch-off-barrier` are phenomenological shapes until Step 2.
- Cross-platform bit-identity of the arrays is asserted by CI for the first time in PR #19; if `gallery-hashes-agree`
  fails, the BLAS-dependent paths (`@` in `induced_charge`, `einsum`, `eigh`) are the suspects.
- `docs/catalog/entries/interface-fixed-charge.yaml` still has free-text sources and no `calibration` field
  (warned, not failed, because it is `draft`).
- Step 1 exit criteria still open: sign-offs; the benchmark is met locally (0.02 ms/point on
  `stability-coarse-2dot`, budget 0.1) but not yet confirmed on CI.

## Next action
When CI on PR #19 is green: owner reviews the gallery artifact and the two reference entries; the session iterates
on protocols from the review, then writes sign-offs as instructed. In parallel, start Step 2 on
`step-2/layout-electrostatics` (issue #17): `layout` dataclasses, Davies–Larkin–Stern closed-form potential,
lever arms derived from the layout, `potential-map-layout` and `potential-cut-along-channel` protocols first.
