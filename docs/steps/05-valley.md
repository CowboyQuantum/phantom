# Step 5 — Valley

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
Valley splitting per dot as a statistical draw for the material and dot size, with its consequences for Pauli spin blockade and relaxation.

## Deliverables
- `valley`: a distribution of valley splitting for a Ge profile and dot size (Ge-fluctuation-limited theory); per-dot draws; a hook for an interface-step model
- `charge` and later `spin` consume the draw: PSB window and lifting, relaxation hot spots

## P0 tests (CI: analytic limits and invariants)
- the distribution's mean and spread scale with Ge concentration and dot size as the theory states

## P1 tests (CI, `literature`: reproduce an item of a `confirmed` reference entry; ADR 0007)
- a published valley-splitting distribution for a stated heterostructure reproduces within tolerance

## Gallery protocols (rendered by CI, signed off per protocol by a human; ADR 0008)
- to be listed before the step starts; each carries a `budget` (ADR 0009)

## Exit criteria
- PSB readout quality varies between draws for physical reasons
- the P1 reproduction passes

## Open decisions
- whether to model interface steps explicitly

## Status
NOT STARTED
