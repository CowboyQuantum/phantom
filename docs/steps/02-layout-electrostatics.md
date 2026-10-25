# Step 2 — Layout and closed-form electrostatics

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
A gate layout and a material stack give V(r; gate vector) in microseconds; dots and barriers emerge from the potential; lever-arm fields and the capacitance matrix are derived, not typed.

## Deliverables
- `layout`: gate polygons, stack (2DEG depth, dielectric layers), material class (`SiGe`, `MOS`) with their default parameters; JSON in and out; the textbook `two_dot_layout()` (done)
- `electrostatics`: Davies–Larkin–Stern closed-form potential for polygonal gates over a 2DEG as the solid angle of any polygon; per-gate basis; superposition; dots as minima with harmonic radii, barriers as saddles; a 2D effective-mass Schrödinger solve, offline (done; ADR 0013)
- `charge` reads lever arms and barrier heights from `electrostatics` instead of constants: `charge.derived.derive_double_dot` builds the CI parameters from a layout (done)
- reference entries (`proposed`) for the closed-form potential and for the published layout the P1 row reproduces (ADR 0007)

## P0 tests (CI: analytic limits and invariants)
- a single rectangular gate's potential matches the closed form at known points
- lever arms from the basis equal the capacitance-matrix values for a textbook two-dot layout
- a two-dot layout yields a honeycomb whose geometry matches Step 1's CI model within tolerance

## P1 tests (CI, `literature`: reproduce an item of a `confirmed` reference entry; ADR 0007)
- lever arms and charging energies for a published single-layer SiGe layout within the paper's stated range

## Gallery protocols (rendered by CI, signed off per protocol by a human; ADR 0008)
- `potential-map-layout`: V(x, y) in the 2DEG plane at a named gate setting, gate polygons overlaid; two wells and the barriers between them visible (rendering)
- `potential-cut-along-channel`: V along the dot-to-dot axis for several barrier voltages; the barrier height falls as the gate opens (rendering)
- `stability-coarse-2dot-layout`: the coarse map with nothing typed but the layout, the operating point and the sensor's per-electron coupling; the sensor gate is set on its flank first (rendering)
- `gate-leakage-iv`: leakage current vs. each gate; ~0 in the safe window, onset at the extremes
- `ohmic-iv`: current vs. source–drain bias with all gates accumulated; linear
- `turn-on-global`: current vs. all gates swept together; threshold and saturation, MOS vs. SiGe scale
- `pinch-off-crosstalk`: pinch-off voltage of one barrier vs. another gate's voltage; a straight line whose slope is the cross-capacitance ratio (rendering; step shape phenomenological)
- `crosstalk-lines`: sensor signal vs. every gate pair; the full cross-capacitance matrix from the line slopes
- `gate-hysteresis`: one gate up then down; the two traces overlap for a clean device (the loop appears in Step 3)
- still to do: `gate-leakage-iv`, `ohmic-iv`, `turn-on-global`, `crosstalk-lines`, `gate-hysteresis`; re-rendering `turn-on-single-gate` and `pinch-off-barrier` from the layout
Full roster: `docs/gallery/PROTOCOLS.md`. Each carries a `budget`; each is unsigned until the owner has reviewed it in the PR artifact.

## Exit criteria
- a layout file is the only input needed to produce a stability diagram
- the P1 reproduction passes against a `confirmed` entry
- every protocol above rendered and signed off

## Open decisions
- gate overlap and multi-layer gates: closed-form limits vs. FEM (Step 7)

## Status
IN PROGRESS
