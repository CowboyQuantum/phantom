# Step 2 — Layout and closed-form electrostatics

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
A gate layout and a material stack give V(r; gate vector) in microseconds; dots and barriers emerge from the potential; lever-arm fields and the capacitance matrix are derived, not typed.

## Deliverables
- `layout`: gate polygons, stack (2DEG depth, dielectric layers), material class (`SiGe`, `MOS`) with their default parameters
- `electrostatics`: Davies–Larkin–Stern closed-form potential for polygonal gates over a 2DEG; per-gate basis; superposition; a 2D effective-mass Schrödinger solve on the potential for dot positions and level spacings
- `charge` reads lever arms and barrier heights from `electrostatics` instead of constants

## P0 tests (CI: analytic limits and invariants)
- a single rectangular gate's potential matches the closed form at known points
- lever arms from the basis equal the capacitance-matrix values for a textbook two-dot layout
- a two-dot layout yields a honeycomb whose geometry matches Step 1's CI model within tolerance

## P1 tests (CI, `literature`: reproduce an item of a `confirmed` reference entry; ADR 0007)
- lever arms and charging energies for a published single-layer SiGe layout within the paper's stated range

## Gallery protocols (rendered by CI, signed off per protocol by a human; ADR 0008)
- to be listed before the step starts; each carries a `budget` (ADR 0009)

## Exit criteria
- a layout file is the only input needed to produce a stability diagram
- the P1 reproduction passes

## Open decisions
- gate overlap and multi-layer gates: closed-form limits vs. FEM (Step 7)

## Status
NOT STARTED
