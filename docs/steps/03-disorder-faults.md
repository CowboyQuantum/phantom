# Step 3 — Disorder and the fault models

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
Catalog entries become generators that perturb the compiled potential and parameters, so unintentional dots, leaky barriers, asymmetric lever arms, charge switchers and drift appear for physical reasons and at catalog-set rates.

## Deliverables
- `disorder`: static random interface/oxide charges (density per material class) added to the potential basis; gate misalignment and width error as basis transforms; strain and g-factor perturbations as parameter draws
- fault generators for every catalog entry marked `runtime`, each with a severity knob
- a scenario draw: material + catalog + seed → a perturbed device, serialisable

## P0 tests (CI: analytic limits and invariants)
- each generator's perturbation has the documented sign, scale and spatial extent
- a draw with zero severity reproduces the clean device bit for bit

## P1 tests (CI, `literature`: reproduce an item of a `confirmed` reference entry; ADR 0007)
- the density of unintentional dots per unit length vs. interface charge density matches a published MOS/SiGe comparison within an order of magnitude

## Gallery protocols (rendered by CI, signed off per protocol by a human; ADR 0008)
- to be listed before the step starts; each carries a `budget` (ADR 0009)

## Exit criteria
- every catalog entry marked `runtime` has a generator and a P0 test showing its signature in a stability diagram
- the fault signatures listed in the catalog are visible in rendered diagrams reviewed by the human

## Open decisions
- which catalog entries need the offline tier (Step 7) rather than a runtime generator

## Status
NOT STARTED
