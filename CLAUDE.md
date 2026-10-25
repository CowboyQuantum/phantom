# CLAUDE.md — phantom

Standing law for every agent session. Short on purpose: rules and where to look, not explanations.
Why: `docs/adr/`. What next: `docs/STATUS.md`. The protocol is the one qdt uses (and `jig` will carry).

## Session protocol (mandatory)
**Start** — read, in this order, nothing else until the task requires it: this file; `docs/STATUS.md`;
`gh issue view <n>` for the issue STATUS names; the `docs/steps/` file that issue links; the layer's
docstring in `src/phantom/<layer>/__init__.py`. Restate the task in one sentence; if it differs from the
issue, stop and say so.
**During** — one branch per issue: `step-N/<short-name>` for a roadmap-step issue, `issue-<n>/<short-name>` otherwise; conventional commits scoped by layer
(`feat(charge): …`, `test(disorder): …`, `docs(adr): …`); PR body says `Closes #<n>`.
**End** — `uv run pytest` and `bash scripts/check_all.sh` green; rewrite `docs/STATUS.md` (never append);
`docs/MODELS.md` updated if a model was added or rewired, and `scripts/gallery.py --record` re-run if a model or
protocol changed (ADR 0011); an ADR for any decision a future session
could ask "why?" about; comment on the issue; commit.

## Hard rules (CI enforces; see `scripts/`)
1. `phantom` imports nothing from qdt, from any measurement software, or from any instrument SDK. Ever.
   It is numpy/scipy in, numpy out. (ADR 0003)
2. Layers import only downward: `layout` ← `electrostatics` ← `disorder` ← `valley` ← `charge` ←
   `readout` ← `spin` ← `optimum`; `compile` may import everything; nothing imports `compile`. (ADR 0004)
3. **Reproducibility contract:** every model takes an explicit seed and carries a `model_version`; the same
   parameters, seed and version give identical output, bit for bit, on the same platform and numpy build, and
   agree within 1e-12 relative across platforms; `ubuntu-latest` in CI is the reference for sign-off hashes. (ADR 0003, 0010)
4. **Compile once, run cheap:** anything slower than ~1 ms per gate point belongs in `compile`, which
   produces parameters the runtime layers consume. Runtime layers never call a solver. (ADR 0001)
5. Every physical quantity carries its unit in the name or the docstring; parameters are dataclasses with
   documented ranges and a reference id where one exists.
6. Provenance: no real device, wafer, site, person or proprietary parameter values. Published values cite
   their reference entry; internal device statistics never enter this repo and are not a check. (ADR 0005, 0007)
7. Nothing in this repo knows what a user, score, challenge, agent or instrument is. Scenarios here are
   *physics draws*; gallery protocols are *sweeps*, not tuning steps. What is done with them is a consumer's business.
8. **Published results are the only standard.** Every claim of agreement with physics names an id in
   `docs/references/`; a `literature` test loads its curve from there; a session writes entries as
   `proposed`, only a human sets `confirmed`. (ADR 0007)
9. **The gallery is signed by a human.** Every protocol in `docs/gallery/protocols/` has a `budget`; CI
   renders it and checks its array hash against `docs/gallery/signoff/`; a session writes a sign-off only
   when the human says so. No rendered image is committed; every PR publishes a preview at `/pr-<n>/` on the
   Pages site and review happens there. (ADR 0008, 0009, 0012)

## Discipline
- **Two-rung ladder** (ADR 0007): P0 analytic limits and invariants in CI; P1 reproductions of
  `confirmed` reference items, in CI, marked `literature(id, item)`. No P2. A quantity with no published
  standard is `calibration: none` and is never called realistic.
- **Step exit** = P0 and P1 green **and** every gallery protocol the step's layers touch signed off.
- **Budgets** (ADR 0009): counts exact, wall time within 3× in CI; raising a budget is explained in the PR.
- **Statistical realism, not prediction** (ADR 0002): a phantom behaves like a device of its class; its
  disorder is a draw. Never claim to predict a particular device.
- **The defect catalog is data** (`docs/catalog/`, ADR 0006): a human fills it; code reads it.
- **qdt is one consumer among several.** Its adapter lives in qdt (`qdt-sim`), pins a phantom release,
  and maps qdt's objects onto phantom's. Changes here that break the adapter are breaking changes.

## Commands
```
uv sync --all-extras                      # dev install
uv run pytest                             # P0 + P1 (what CI runs)
bash scripts/check_all.sh                 # layer imports + no-qdt + provenance + references + model graph + lock + lint
uv run python scripts/gallery.py          # render protocols to build/gallery/ and verify sign-offs (Step 1+)
uv run mkdocs serve                       # browse docs and the rendered gallery
```

## Where to look
| Need | Read |
|---|---|
| A map of the layers and the two tiers | `docs/DEVELOPMENT_PLAN.md` §2–§4 |
| The current step's spec | `docs/steps/<NN-name>.md` (index: `docs/steps/README.md`) |
| Why it is built this way | `docs/adr/` |
| The defect and material catalog | `docs/catalog/` |
| How the models chain into a measurement | `docs/MODELS.md` (ADR 0011) |
| The papers that are the standard | `docs/references/` (schema in `SCHEMA.md`) |
| The protocols, sign-offs and budgets | `docs/gallery/` |
| How qdt consumes phantom | ADR 0003 |
