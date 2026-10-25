# 0001 — Two tiers: compile once per device instance, run cheap at every gate point

**Status:** Accepted
**Date:** 2026-10-24
**Step:** all

## Context
Faithful electrostatics (3D Poisson with screening) costs seconds to an hour per gate configuration;
atomistic valley splitting costs hours; a tune-up sweeps ten thousand gate points. No model that solves
physics at runtime can serve a tuning loop. Laplace is linear in the gate voltages, so a basis of potentials
can be solved once; Hubbard parameters are smooth functions of barrier height; valley splitting is a draw.

## Decision
Two tiers. `compile` runs the expensive physics once per device instance and emits a serialisable
`Device` of parameter fields and functions; the runtime layers evaluate it at a gate point in under ~1 ms
and never call a solver. Closed-form electrostatics (polygonal gates over a 2DEG) is the first compile
path; FEM through a wrapped external solver is the second. Anything slower than the runtime budget moves
to `compile`.

## Consequences
- Easier: ≥ 10⁴ gate points per second per core; scenarios are served pre-compiled; the runtime is pure
  numpy and portable.
- Harder: two code paths for electrostatics; a compiled device must be versioned and serialised; nonlinear
  screening is approximated at runtime (parameter functions) rather than solved.
- Must now: a benchmark in CI; `model_version` on every compiled device; the runtime budget in `CLAUDE.md`.

## Alternatives considered
- Solve physics at runtime with a fast solver — rejected: still ~10³× too slow for a loop.
- Precompute full stability maps — rejected: the loop needs dynamics (rates, noise, latching), not images.
