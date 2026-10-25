# 0005 — Test ladder P0 (analytic limits), P1 (published results), P2 (device statistics, internal)

**Status:** Superseded by 0007 (P2 removed; P0 and P1 stand)
**Date:** 2026-10-24
**Step:** all

## Context
A physics library needs its own notion of "proven," distinct from qdt's instrument ladder. The claims are
statistical (ADR 0002), so the highest level compares distributions with real device statistics that are
not public.

## Decision
- **P0** — analytic limits and invariants, in CI on every commit (Coulomb-diamond slopes, honeycomb
  geometry from a capacitance matrix, detailed balance, trace preservation, zero-severity draws reproduce
  the clean device bit for bit).
- **P1** — reproductions of published results, cited by DOI, in CI, marked `literature`, with the paper's
  stated tolerance.
- **P2** — comparisons with real device statistics (fault rates, sensor behaviour, tune-up time
  distributions), marked `internal`, run where that data lives, never committed here, and never an exit
  criterion for a step.
A step exits when its P0 and P1 rows pass.

## Consequences
- Easier: the repo stays public and provenance-clean; fidelity has a citable basis.
- Harder: some claims (fault rates) are only checkable at P2.
- Must now: every P1 test names its source; `-m 'not internal'` is the default pytest selection.

## Alternatives considered
- Only P0 — rejected: invariants do not establish realism.
