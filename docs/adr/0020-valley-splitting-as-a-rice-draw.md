# 0020 — Valley splitting is a Rice draw per dot from alloy disorder, with a deterministic part for MOS; three signatures

**Status:** Accepted
**Date:** 2026-10-26
**Step:** 5

## Context
Silicon dots have two low-lying valleys whose splitting varies from dot to dot and sets whether Pauli
spin blockade works and where spin relaxation has hot spots. Atomistic calculation is out of scope
(plan §3). The alloy-disorder theory of Si/SiGe wells gives a statistical description that is both
simple and the one the current literature uses to compare devices.

## Decision
- **Distribution** (`valley.splitting`): the valley coupling is a complex number Δ_det + σ(n₁ + i n₂)
  with n₁, n₂ standard normal; the splitting |Δ| follows a Rice distribution. σ scales as 1/ℓ with the
  dot radius, as √(x̄(1 − x̄)) with the germanium concentration the electron samples, and as the square
  root of the interface width (overlap with the alloy); Δ_det is small for SiGe and, for MOS, large and
  linear in the vertical field with a smaller random part from roughness. Reference values are
  placeholders (`calibration: none`) until `paqueletwuetz2022-valley-alloy`, `losert2023-valley-strategies`
  and `yang2013-mos-valley` are confirmed; the P1 row sets σ and Δ_det from the measured distribution.
- **Interface steps** are not modelled explicitly (open decision of the step file resolved for now): a
  `step_factor` hook scales σ, so a later step model can plug in without changing the draw's shape.
- **Signatures** (`valley.signatures`): the second electron's addition energy rises as +½ g μ_B B until
  g μ_B B = E_V and falls after (magnetospectroscopy kink at B* = E_V / g μ_B); the spin relaxation rate
  has a Lorentzian hot spot where E_Z = E_V on a background that the spin layer will make field-dependent;
  pulsed-gate spectroscopy splits a transition by the pulse amplitude and adds a weaker excited-state
  line E_V / α away (idealised, fixed weight).
- **Per-dot draws** take the dot's harmonic radius from `electrostatics.dots`; `compile` (Step 7) will
  attach them to the compiled device, and `spin` (Step 6) consumes them for PSB.

## Consequences
- Easier: a device's valley splittings are draws with the right dot-to-dot spread and the right size and
  material dependence; the fraction of dots unusable for PSB falls out of the distribution.
- Harder: the deterministic part and the absolute scale are placeholders; interface steps are a hook.
- Must now: P0 tests for the scalings and the Rayleigh/Rice limits; protocols `valley-splitting-draws`,
  `magnetospectroscopy-transition`, `t1-hot-spot`, `pulsed-gate-spectroscopy`.

## Alternatives considered
- A fixed splitting per material — rejected: ADR 0002, and PSB statistics are the point.
- Tight-binding valley calculations — out of scope by plan §3.
