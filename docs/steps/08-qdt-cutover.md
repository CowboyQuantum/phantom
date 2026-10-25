# Step 8 — qdt cutover

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
qdt's `qdt-sim` becomes an adapter that pins a phantom release and contains no physics; qdt's golden files are regenerated once; the seam is enforced on both sides.

## Deliverables
- in phantom: a release tag and changelog; the consumer contract (plan §4) documented as the public API
- in qdt (coordinated issue): `qdt-sim` maps `StandConfig`/`ChannelMap` onto `phantom.Device`, implements the protocols over it, pins the release; the physics modules deleted; `make_sim_fixtures.py` regenerated; the vocabulary and import checks updated

## P0 tests (CI: analytic limits and invariants)
- phantom's API surface is covered by a contract test that qdt's adapter also runs

## P1 tests (CI, `literature`: published results, cited)
- none

## P2 (internal, never an exit criterion)
- comparisons with device statistics, when available

## Exit criteria
- qdt CI green against the pinned phantom with `qdt-sim` under 300 lines
- qdt's tune-up on the simulator reproduces its committed scenario outcomes after regeneration

## Open decisions
- cadence of phantom releases vs. qdt pins

## Status
NOT STARTED
