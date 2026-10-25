# CLAUDE.md — phantom

Standing law for every agent session. Short on purpose: rules and where to look, not explanations.
Why: `docs/adr/`. What next: `docs/STATUS.md`. The protocol is the one qdt uses (and `jig` will carry).

## Session protocol (mandatory)
**Start** — read, in this order, nothing else until the task requires it: this file; `docs/STATUS.md`;
`gh issue view <n>` for the issue STATUS names; the `docs/steps/` file that issue links; the layer's
docstring in `src/phantom/<layer>/__init__.py`. Restate the task in one sentence; if it differs from the
issue, stop and say so.
**During** — one branch per issue (`issue-<n>/<short-name>`); conventional commits scoped by layer
(`feat(charge): …`, `test(disorder): …`, `docs(adr): …`); PR body says `Closes #<n>`.
**End** — `uv run pytest` and `bash scripts/check_all.sh` green; rewrite `docs/STATUS.md` (never append);
an ADR for any decision a future session could ask "why?" about; comment on the issue; commit.

## Hard rules (CI enforces; see `scripts/`)
1. `phantom` imports nothing from qdt, from any measurement software, or from any instrument SDK. Ever.
   It is numpy/scipy in, numpy out. (ADR 0003)
2. Layers import only downward: `layout` ← `electrostatics` ← `disorder` ← `valley` ← `charge` ←
   `readout` ← `spin` ← `optimum`; `compile` may import everything; nothing imports `compile`. (ADR 0004)
3. **Reproducibility contract:** every model takes an explicit seed and carries a `model_version`; the same
   parameters, seed and version give identical output, bit for bit, on every platform CI runs. (ADR 0003)
4. **Compile once, run cheap:** anything slower than ~1 ms per gate point belongs in `compile`, which
   produces parameters the runtime layers consume. Runtime layers never call a solver. (ADR 0001)
5. Every physical quantity carries its unit in the name or the docstring; parameters are dataclasses with
   documented ranges and a literature source where one exists.
6. Provenance: no real device, wafer, site, person or proprietary parameter values. Published values cite
   their paper; internal statistics (P2) never enter this repo. (ADR 0005)
7. Nothing in this repo knows what a user, score, challenge, agent or instrument is. Scenarios here are
   *physics draws*; what is done with them is a consumer's business.

## Discipline
- **Test ladder P0/P1/P2** (ADR 0005): P0 analytic limits and invariants in CI; P1 reproductions of
  published results, in CI, marked `literature`; P2 comparisons with device statistics, internal, marked
  `internal`, never an exit criterion here.
- **Statistical realism, not prediction** (ADR 0002): a phantom behaves like a device of its class; its
  disorder is a draw. Never claim to predict a particular device.
- **The defect catalog is data** (`docs/catalog/`, ADR 0006): a human fills it; code reads it.
- **qdt is one consumer among several.** Its adapter lives in qdt (`qdt-sim`), pins a phantom release,
  and maps qdt's objects onto phantom's. Changes here that break the adapter are breaking changes.

## Commands
```
uv sync --all-extras                      # dev install
uv run pytest                             # P0 + P1 (what CI runs)
uv run pytest -m internal                 # P2, with the internal statistics mounted
bash scripts/check_all.sh                 # layer imports + no-qdt + provenance + lock + lint
```

## Where to look
| Need | Read |
|---|---|
| A map of the layers and the two tiers | `docs/DEVELOPMENT_PLAN.md` §2–§4 |
| The current step's spec | `docs/steps/<NN-name>.md` (index: `docs/steps/README.md`) |
| Why it is built this way | `docs/adr/` |
| The defect and material catalog | `docs/catalog/` |
| How qdt consumes phantom | ADR 0003 |
