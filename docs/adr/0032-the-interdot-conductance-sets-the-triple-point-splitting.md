# 0032 — The interdot conductance sets the triple-point splitting: Kemble's transmission of the interdot barrier, Golden & Halperin's f(g), an effective interdot capacitance in the CI model

**Status:** Accepted
**Date:** 2026-10-31
**Step:** 8 follow-up (issue #66; after ADR 0029–0031)

## Context
After #58–#60 the owner still did not see the behaviour they were tuning for: "if I increase the barrier
gate between plunger gates, the interdot width (in anti-detuning) should get bigger until I form a
superdot underneath the barrier gate. The opposite should also be true, a lower barrier gate voltage should
shrink the interdot width until all I see is a grid" (2026-10-31). That width is the triple points'
separation. On the reference device, the fractional splitting f = 2C_m/(C_t + C_m) moved only from 0.14
to 0.16 over ±0.2 V of B_M:

- **Geometry.** C_m came from the sphere formula of ADR 0013, which B_M changed by about 0.7 aF.
- **Hybridization.** ADR 0029's coupling rounds the interdot line but shifts the triple points by at most
  2t/E_add ≈ 0.15.
- **No superdot.** A single tunnelling level cannot produce the merged dot at all.

In double-dot transport the splitting follows the interdot conductance through quantum charge
fluctuations (Waugh et al. 1995; review: van der Wiel et al. 2003, the honeycomb from a grid to one dot).
Golden & Halperin give f(g) for N_ch channels of conductance g each, with a weak and a strong limit that
do not interpolate reliably for N_ch = 2.

The owner chose (2026-10-31):
- a conductance-driven splitting entering the CI model as an effective interdot capacitance;
- a smooth interpolation between the published limits;
- folding it into 0.4.0.

Four forces shaped the design:

- **Continuity to the merge.** The conductance has to reach 1 continuously as the dots merge. WKB's
  e^(−2S) stops at S = 0, and the rates' logistic halves it.
- **Bit-exact limit.** Without conductance nothing may change.
- **Gate charge.** The tunnelling relaxes charge quantization. It must not move the gates' induced charge.
- **References.** Every published curve enters as a `proposed` reference entry (ADR 0007).

## Decision
The interdot conductance per channel is Kemble's transmission of the interdot barrier, and Golden &
Halperin's f(g) sets an effective interdot capacitance in the CI model. Specifically:

- **The barrier.**
  - Compile records the top ΔU of the interdot barrier above its fixed tunnelling level (ADR 0027),
    negative below it: the highest U on the straight cut between the dots, sampled bilinearly.
  - It also records ΔU's slopes against every gate, with the dots held, and its curvature ħω_b (a
    parabola fitted ±4 nm around the top). The shape grid gains a `barrier_top_ev` column, the top with
    the dots moved.
  - Device fields `interdot_barrier_top_ev`, `interdot_barrier_top_slopes_ev_per_v` and
    `interdot_barrier_hbar_omega_ev`; `DEVICE_FORMAT_VERSION` 5.
- **The conductance.** g = 1/(1 + exp(2πΔU/ħω_b)) (Kemble), with N_ch = 2 (spin). It is ½ when the
  level reaches the top, and it goes continuously towards 1 as the barrier sinks to the merge.
- **The splitting.**
  - f_GH(g) joins (2 ln 2/π²)·N_ch·g (g ≤ 0.15) and 1 + (16e^γ/π³)(1 − g)·ln(1 − g) (g ≥ 0.85) with a
    monotone cubic (`calibration: none`).
  - The total is f = f_C + (F_max − f_C)·f_GH(g), with f_C = 2C_m/(C_t + C_m) of the geometric
    capacitances and F_max = 0.98.
- **The CI model.**
  - `ConstantInteractionParams.splitting` adds ΔC_m·[[1, −1], [−1, 1]] to the capacitances in the
    energies, with ΔC_m = (F_max − f_C)·f_GH·(C_t + C_m)/(2(1 − f)). That is exactly zero at g = 0.
  - The gate-induced charge keeps the geometric capacitances.
  - The CI model's measured splitting along the symmetric axis is f (P0).
- **The consumer surface.** `State.interdot_conductance` and `State.interdot_splitting`;
  `Device.runtime_splitting()`.
- **References.** `golden1996-peak-splitting`, `waugh1995-tunable-coupling` and `kemble1935-wkb-barrier`
  as `proposed` entries, and the honeycomb evolution as an item of `vanderwiel2003-double-dot-transport`.
- **Release.** 0.4.0. Model versions: `charge.splitting` 1, `charge.shape` 2, `compile.device` 6,
  `device` 7, `device.session` 8.

## Consequences
- Easier:
  - Opening B_M takes the honeycomb from small hexagons, through elongated hexagons, to the diagonal lines
    of one superdot, and closing it shrinks the segments back (`honeycomb-vs-barrier`).
  - On the reference device g rises from 0.002 at the operating point to 0.99 at B_M +50 mV, and f from
    0.16 to 0.94. Where the dots merge, f holds at 0.98.
  - The exchange (ADR 0031) and every other CI quantity see the same capacitances.
- Harder:
  - **No exact grid.** Closing B_M shrinks f only to the geometric f_C ≈ 0.14. The sphere-formula C_m has no
    screening by the barrier gate (#48), so the "grid" is small hexagons, not a grid.
  - **Steep rise.** f rises within about 40 mV of B_M. Kemble's transmission changes by e for every
    ħω_b/2π ≈ 0.7 meV of barrier height, and the unscreened barrier moves fast with the gates (#48).
  - **The f(g) curve.** Between g = 0.15 and 0.85 it is an interpolation, which the authors themselves
    call unreliable for N_ch = 2. A digitised measurement (Waugh et al.) would replace it as a P1 row.
  - **Mixed models.** The splitting is the metallic-dot theory (many levels) applied to semiconductor dots
    of about 23 electrons. It sits on top of ADR 0029's single-level hybridization, and the two overlap by
    at most 2t/E_add.
  - **Asymmetry.** For unequal dots f uses the mean total capacitance (the symmetric-equivalent
    definition).
- Must now:
  - Device files from earlier 0.4.0 commits do not load (format 5).
  - A test that checks the shape's own capacitances turns the splitting off.

## Alternatives considered
- **A better electrostatic C_m from the potential (overlapping charge densities).** Rejected by the owner.
  It is closer to the geometry, but it cannot reach the superdot without a solver, and it moves f little.
- **Fitting measured f(g) (Waugh, Livermore) now.** Deferred. The curve has to be digitised and confirmed by
  a human before it can be the standard (ADR 0007).
- **The WKB transmission e^(−2S), or the rates' logistic, as the conductance.** Rejected. Either stops at a
  half (or at one, abruptly) when the barrier falls below the level, so f would jump at the merge.
- **Raising t's ceiling, or several levels per dot.** Rejected. The single-level splitting is 2t, and
  the merged dot needs the many-level charge-fluctuation physics that f(g) summarises.
- **Moving the induced charge with the effective capacitance.** Rejected. Tunnelling relaxes charge
  quantization and does not change how much charge the gates induce.
