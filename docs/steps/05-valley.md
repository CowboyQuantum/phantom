# Step 5 — Valley

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
Valley splitting per dot as a statistical draw for the material and dot size, with its consequences for Pauli spin blockade and relaxation.

## Deliverables
- `valley`: a distribution of valley splitting for a Ge profile and dot size (Ge-fluctuation-limited theory); per-dot draws; a hook for an interface-step model
- `charge` and later `spin` consume the draw: PSB window and lifting, relaxation hot spots

## P0 tests (CI: analytic limits and invariants)
- the distribution's mean and spread scale with Ge concentration and dot size as the theory states

## P1 tests (CI, `literature`: published results, cited)
- a published valley-splitting distribution for a stated heterostructure reproduces within tolerance

## P2 (internal, never an exit criterion)
- comparisons with device statistics, when available

## Exit criteria
- PSB readout quality varies between draws for physical reasons
- the P1 reproduction passes

## Open decisions
- whether to model interface steps explicitly

## Status
NOT STARTED
