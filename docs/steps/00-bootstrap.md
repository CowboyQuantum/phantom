# Step 0 — Bootstrap

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
Repo exists; CI green on an empty P0 suite; layer skeletons with docstrings; catalog template; ADR 0001–0006; the issue set for the catalog filed.

## Deliverables
- uv project, `src/phantom/<layer>/__init__.py` for every layer with a docstring stating what it owns
- `scripts/check_imports.py` (layers downward only; no qdt, no instrument SDKs), `check_all.sh`, CI
- `docs/catalog/` template and schema; `scripts/issues_catalog.py` filed
- Full Apache-2.0 LICENSE text; remote on `CowboyQuantum/phantom`

## P0 tests (CI: analytic limits and invariants)
- workspace installs; checks and an empty pytest pass in CI

## P1 tests (CI, `literature`: published results, cited)
- none

## P2 (internal, never an exit criterion)
- comparisons with device statistics, when available

## Exit criteria
- CI green on `main`
- catalog issues exist on GitHub
- ADR 0001–0006 merged

## Open decisions
- plan §8 #4 catalog format

## Status
DONE
