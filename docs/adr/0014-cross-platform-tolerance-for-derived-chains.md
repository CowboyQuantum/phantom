# 0014 — The cross-platform tolerance is 1e-9 relative, because derived parameter chains amplify last-bit differences

**Status:** Accepted (amends 0010)
**Date:** 2026-10-26
**Step:** 2

## Context
ADR 0010 set the cross-platform agreement of gallery arrays at 1e-12 relative, measured per array as the
largest difference over the largest magnitude. Every Step 1 protocol met it with margin (worst 3e-14).
The first Step 2 protocols that derive their parameters from a layout did not: `stability-coarse-2dot-layout`
and `crosstalk-lines` differ between ubuntu-latest and windows-latest by 1.5e-11 and 2.0e-11. The chain
solid angle (arctangents) → lever arms → capacitance matrix inverse → gate-voltage of the sensor flank →
Boltzmann occupation turns a last-bit difference in an elementary function into a shift of a transition
position of order 1e-14 V, and at a pixel that sits on a sharp step that shift is a relative change of
order 1e-11 in the current. Nothing is wrong with either platform's result; the tolerance was tighter
than the arithmetic allows for derived chains with sharp features.

## Decision
- The cross-platform tolerance (`CROSS_PLATFORM_RTOL` in `gallery/runner.py`) is 1e-9 relative, per
  array, with the same definition as ADR 0010. A genuine change to physics or parameters moves arrays by
  1e-6 or more, so the check still catches every real change; what it stops flagging is arithmetic noise
  amplified through derived chains.
- Everything else in ADR 0010 stands: bit-for-bit per platform, ubuntu-latest the reference for sign-off
  hashes, the Windows job compared to the reference.
- Rounding derived parameters to suppress the noise was considered and rejected again for the reason in
  ADR 0010: a value on a rounding boundary flips.

## Consequences
- Easier: layout-derived protocols pass on both platforms without special cases.
- Harder: a change that moves an array by less than 1e-9 relative is invisible to the Windows check;
  the reference platform's exact hash (sign-off) still sees it.
- Must now: CLAUDE.md rule 3 and the gallery README say 1e-9; the measured worst cases are recorded here
  (3e-14 for typed chains, 2e-11 for derived chains, both on 2026-10-26).

## Alternatives considered
- Per-protocol tolerances — rejected: one number is easier to reason about and 1e-9 is loose enough.
- Comparing occupations instead of currents — rejected: the sensor current is what the protocol measures.
