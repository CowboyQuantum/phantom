# 0011 — A running map of how the models fit together lives in the docs and is kept current

**Status:** Accepted
**Date:** 2026-10-25
**Step:** 1 onward

## Context
By the end of Step 1 the readout layer already has seven modules and the charge layer three, and a
measurement is produced by chaining several of them (occupation, lead rates, latching, sensor, noise,
demodulation). As later steps add electrostatics, disorder, valley and spin, a human reading the repo
needs one place that says how these pieces combine to produce a measurement, and that place must not
drift from the code.

## Decision
- `docs/MODELS.md` holds two diagrams. The **measurement flow** is drawn by hand in Mermaid: boxes are
  models, arrows say which feeds which, placeholders that a later step replaces are dashed and name
  that step, and a table maps every gallery protocol to its path through the boxes. The **module
  dependency graph** is generated from the imports in `src/phantom` by `scripts/model_graph.py`.
- A pull request that adds, removes or rewires a model updates the hand-drawn diagram and the protocol
  table in the same PR. `check_all.sh` runs `scripts/model_graph.py --check`, so the generated graph
  can never be stale; the hand-drawn one is reviewed like any other doc.
- The page is in the mkdocs nav as "Models" and renders the Mermaid diagrams on the docs site.
- **Runtimes are recorded on the same page.** `scripts/gallery.py --record` writes
  `docs/gallery/measured.json`: the measured cost of every protocol in that run and a per-model
  micro-benchmark (microseconds per point for each runtime model), with commit, platform and date.
  `scripts/model_graph.py` renders it into the page beside the declared budgets, and `--check` keeps the
  table in sync with the file. A session re-records when it changes a model or a protocol. The numbers are
  from the recording machine; the reference platform's numbers stay in the CI artifact (ADR 0010).

## Consequences
- Easier: a newcomer, or the owner reviewing a gallery picture, can trace which models produced it.
- Harder: one more file to touch when a model changes; the hand-drawn diagram relies on review.
- Must now: CLAUDE.md names the file in the End protocol and in "Where to look".

## Alternatives considered
- Only the generated import graph — rejected: imports do not say what flows where, or what a placeholder is.
- Only the hand-drawn diagram — rejected: it would drift; the generated graph anchors it to the code.
