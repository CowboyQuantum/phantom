# Step 7 — Compile and optimum

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
The offline tier: FEM electrostatics for real stacks and overlapping gates through a wrapped external solver, Thomas–Fermi screening, and a `compile()` that produces a serialisable device; `optimum()` for every layer.

## Deliverables
- `compile`: layout + material + catalog draw + seed → `Device` (versioned, serialised); closed-form or FEM electrostatics selectable
- `optimum`: best operating point for sensor sensitivity, best lattice (virtual-gate matrix), target interdot coupling; with `spin`, achievable readout fidelity and coherence
- a solver wrapper behind the `[fem]` extra; never a solver written here

## P0 tests (CI: analytic limits and invariants)
- closed-form and FEM agree for a shallow single-layer layout within tolerance
- `optimum` is attained by a brute-force search on a small device

## P1 tests (CI, `literature`: reproduce an item of a `confirmed` reference entry; ADR 0007)
- none

## Gallery protocols (rendered by CI, signed off per protocol by a human; ADR 0008)
- to be listed before the step starts; each carries a `budget` (ADR 0009)

## Exit criteria
- a compiled device round-trips through serialisation and reproduces its diagrams
- `optimum` exists for every layer present

## Open decisions
- plan §8 #1 solver choice; #3 serialisation format

## Status
NOT STARTED
