# 0009 — Every gallery protocol carries a time and compute budget, measured in CI

**Status:** Accepted
**Date:** 2026-10-25
**Step:** all, from Step 1

## Context
A main application of phantom is the backend of a simulated measurement suite. A consumer sweeping a
gate at interactive rates needs to know, per kind of measurement, how long a point costs and how much
computation it buys. ADR 0001 sets the rule (compile once, run cheap, ~1 ms per gate point) but nothing
tracks it per measurement, and regressions would be found by the consumer.

## Decision
- **Budget per protocol.** Each protocol file (ADR 0008) declares a `budget` block with two phases,
  matching ADR 0001's tiers:
  - `compile`: wall seconds allowed for the once-per-instance setup, and the number of solver calls it
    may make.
  - `run`: wall milliseconds per point, the number of points, the number of model evaluations per point,
    and `solver_calls: 0` always (runtime layers never call a solver, rule 4).
- **The gallery runner measures.** `scripts/gallery.py` times both phases and counts evaluations and
  solver calls through a counter in the runner (not in the library), and writes `budget.json` with
  declared and measured values, per protocol, per platform.
- **Counts are strict, time is loose.** Evaluation and solver counts are deterministic and CI compares
  them exactly. Wall time on shared runners is noisy, so CI fails only above a stated factor of the
  budget (3×, written in the runner), and prints the full table every run. Locally `--strict` applies
  the factor 1 for benchmarking.
- **Visible.** The budget table is part of the rendered gallery page and the PR artifact, next to the
  pictures, so cost is reviewed with realism.
- **Changing a budget is a decision.** Raising a protocol's budget goes in the PR description with the
  reason; lowering it is free.

## Consequences
- Easier: a consumer can read the cost of each measurement kind from one table; a change that makes a
  sweep slower is caught at the PR, with its cause attributed to a protocol.
- Harder: every protocol author must think about cost up front; the runner carries counting hooks.
- Must now: `budget` is a required block in the protocol schema; Step 1's benchmark deliverable (10⁴
  gate points/s/core) is expressed as the budget of the two-dot stability protocol.

## Alternatives considered
- A single global benchmark — rejected: it hides which measurement regressed and says nothing to a
  consumer choosing a sweep.
- Strict wall-time gates in CI — rejected: flaky on shared runners; counts carry the determinism.
