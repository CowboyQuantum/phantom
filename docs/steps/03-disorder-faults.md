# Step 3 — Disorder and the fault models

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
Catalog entries become generators that perturb the compiled potential and parameters, so unintentional dots, leaky barriers, asymmetric lever arms, charge switchers and drift appear for physical reasons and at catalog-set rates.

## Deliverables
- `disorder`: static random interface/oxide charges (density per material class) added to the potential basis as image-charge energy offsets (done, ADR 0015); gate misalignment and width error as polygon transforms (done); strain and g-factor perturbations as parameter draws (waits for `spin`)
- fault generators for every catalog entry marked `runtime`, each with a severity knob: interface charges, misalignment and width error, open and shorted gates, gate-to-gate leakage, sensor-coupling asymmetry, trapping hysteresis (done on placeholders, ADR 0015/0016); valley, strain, dopants, drift wait for their layers or Step 4
- a scenario draw: material + catalog + seed → a perturbed device, serialisable (`disorder.scenario`, done; parameters are the catalog draft's placeholders until the entries are written)

## P0 tests (CI: analytic limits and invariants)
- each generator's perturbation has the documented sign, scale and spatial extent
- a draw with zero severity reproduces the clean device bit for bit

## P1 tests (CI, `literature`: reproduce an item of a `confirmed` reference entry; ADR 0007)
- the density of unintentional dots per unit length vs. interface charge density matches a published MOS/SiGe comparison within an order of magnitude

## Gallery protocols (rendered by CI, signed off per protocol by a human; ADR 0008)
- `potential-map-disorder`: the landscape of the clean layout beside three draws, gate outlines overlaid (rendering)
- `stability-coarse-2dot-disorder`: the layout-derived coarse map, clean beside three draws each re-derived from its own landscape (rendering)
- `unintentional-dots`: mean extra minima in the channel against charge density over many draws, MOS class (rendering; the P1 row's curve)
- `fault-signatures-stability`: the coarse map with clean wiring and three named faults: a weakly coupled dot, a shorted plunger–barrier pair, an open barrier (rendering)
- `gate-leakage-pair`: a gigaohm leak between two neighbouring gates against the healthy curve (rendering)
- `gate-hysteresis-trapping`: the turn-on swept up and down with trapping, the loop against the clean baseline (rendering)
- to come: switcher and drift rates tied to the catalog (Step 4); fault signatures in the lock-in and plunger-barrier maps as the owner asks

## Exit criteria
- every catalog entry marked `runtime` has a generator and a P0 test showing its signature in a stability diagram
- the fault signatures listed in the catalog are visible in rendered diagrams reviewed by the human

## Open decisions
- which catalog entries need the offline tier (Step 7) rather than a runtime generator

## Status
IN PROGRESS
