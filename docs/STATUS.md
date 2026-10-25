# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** 0 — Bootstrap, DONE → `docs/steps/00-bootstrap.md`. Next to open: Step 1 and Step 2 together.
**Branch:** `step-0/close-bootstrap` (merge to `main`, then branch for Step 1)
**Last session:** 2026-10-25 — closed Step 0: catalog issues filed, full LICENSE, CI confirmed green on both platforms.

## Done
- Step 0 in full. Remote `CowboyQuantum/phantom`, `uv.lock` tracked, CI green on ubuntu and windows (run 40153965906).
- Catalog issues #1–#14 filed by `scripts/issues_catalog.py --apply`; labels `catalog` and `layer:<name>` created.
- `LICENSE` is the canonical Apache-2.0 text (sha256 cfc7749b…23d30); appendix placeholders left as shipped.

## In flight
- This branch, awaiting merge to `main`.

## Known gaps
- No `roadmap`/`blocker` labels and no milestone exist yet; `/session-start` expects them. Create them when Steps 1 and 2 open.
- Steps 1 and 2 have no GitHub issue; CLAUDE.md's protocol wants `gh issue view <n>`. File one roadmap issue per step.
- CLAUDE.md names branches `issue-<n>/…`; `/session-start` names them `step-N/…`. Pick one and fix the other (small ADR or CLAUDE.md edit).
- `docs/catalog/entries/` has one `draft` example; the rest are issues #1–#14 for the human.

## Next action
Merge this branch. Then open Step 1 (`docs/steps/01-charge-readout-core.md`: migrate `qdt-sim`'s double-dot model and
noise into `charge` and `readout`) and Step 2 (`docs/steps/02-layout-electrostatics.md`: layout + closed-form
electrostatics), which can run together: file a roadmap issue for each, create the labels and milestone, then start
Step 1 on `step-1/charge-readout-core`.
