# 0013 — Closed-form electrostatics: solid-angle potentials, effective depth, dots as harmonic minima

**Status:** Accepted
**Date:** 2026-10-25
**Step:** 2

## Context
Step 2 must turn a gate layout into the potential, the lever arms and the capacitance picture that the
charge layer consumes, in microseconds at runtime (ADR 0001), without a field solver. Several modelling
choices are needed to get there; each is an approximation a later step replaces, and each is the kind
of thing a future session would ask "why?" about.

## Decision
- **Potential.** Davies–Larkin–Stern half-space electrostatics (reference entry
  `davies1995-patterned-2deg`): the potential under a gate at voltage V is V times the solid angle the
  gate polygon subtends at the point, over 2π. Implemented for any planar polygon as a fan of triangles
  (Van Oosterom–Strackee), not only rectangles; the rectangle arctangent form is kept as a test.
  Screening by the electrons is not included; Step 7's solver adds it.
- **Stack.** The formula assumes one dielectric. A layered stack enters through an effective depth,
  Σ t_i · ε_well / ε_i: the thin-film approximation, exact for the normal field away from edges.
- **Dots.** Dots are the deepest local minima of the electron potential energy inside the layout's
  channel region; the sensor dot likewise in its region. The dot radius is the harmonic ground-state
  length from the curvature at the minimum; its own capacitance is the disk value from that radius;
  its total capacitance divides by (1 − Σ_j α_j) so that the gate capacitances are part of it; the
  interdot mutual capacitance is the two-sphere estimate 4πε r₁ r₂ / s. Leads are not counted.
- **Lever arms** are the basis values at the dot minimum (the dot as a point). The CI gate
  capacitances are then C_gate = C_dot · α, which makes the CI model's lever arms exactly the
  electrostatic ones; tested as the derivative of the dot energy with respect to each gate voltage.
- **Barriers** are the maximum of the energy along the straight line between two minima, measured
  from the higher minimum.
- **Compile tier.** Building the basis grid and finding dots is compile-tier (milliseconds); the
  runtime consumes the resulting `ConstantInteractionParams`. The finite-difference Schrödinger solve
  exists to check the harmonic estimate and for `compile`; runtime never calls it (rule 4).
- **Validation.** P0 tests cover the closed form against its rectangle special case and known limits,
  the lever-arm derivative identity, the honeycomb period against the derived capacitance, and the
  Schrödinger spectrum against the harmonic oscillator. The P1 row of Step 2 compares lever arms and
  charging energies with a published layout once its reference entry is confirmed.

## Consequences
- Easier: a layout file alone yields potentials, lever arms and a stability diagram; numbers are in the
  right range for silicon dots (lever arms 0.05–0.25, charging energies of order 10 meV).
- Harder: every quantity above the lever arms is an estimate with a stated approximation; fault models
  that depend on screening or on the exact dot shape wait for Step 7.
- Must now: `charge.derived` is the only place the CI parameters are built from a layout; the
  approximations are repeated in its docstring; changes to them bump `MODEL_VERSION`.

## Alternatives considered
- A Poisson solver at runtime — rejected: ADR 0001.
- Rectangles only — rejected: real layouts have rotated and L-shaped gates; the solid angle costs the same.
- Typing the capacitances per layout — rejected: the point of Step 2 is that they are derived.
