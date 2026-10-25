# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** 0 — Bootstrap → `docs/steps/00-bootstrap.md`
**Branch:** `main`
**Last session:** 2026-10-24 (chat) — repo scaffolded from plan v0.1: layers, ADR 0001–0006, catalog template, issue scripts, CI.

## Done
- Everything in Step 0 except: pushing to the remote, the full LICENSE text, filing the catalog issues.

## In flight
- Nothing.

## Known gaps
- `uv.lock` must be generated on first checkout; CI has never run.
- `docs/catalog/entries/` has one `draft` example; the rest are issues for the human.

## Next action
Push; `uv lock`; confirm CI green on both platforms; run `scripts/issues_catalog.py --apply`; then open
Step 1 (migrate `qdt-sim`'s model and noise into `charge` and `readout`) and Step 2 (layout + closed-form
electrostatics), which can run together.
