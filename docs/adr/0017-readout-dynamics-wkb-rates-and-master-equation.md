# 0017 — Readout dynamics: WKB lead rates from the landscape, a master equation with detailed balance, latching as lag

**Status:** Accepted; the rectangular rate estimate is superseded by the WKB action integral (ADR 0023)
**Date:** 2026-10-26
**Step:** 4

## Context
Step 1 modelled rate-limited loading with a placeholder (`readout.loading`): an exponential rate in a
barrier gate and a coin flip per dwell. Step 4 owes the physics: rates that follow from the layout's
barriers, a proper equation of motion for the charge states, and the latching and sweep-direction
effects a real device shows, all cheap enough for the runtime tier (ADR 0001).

## Decision
- **Rates** (`readout.rates`): WKB through a rectangular barrier of the height and width read off the
  landscape along the line from the dot to the reservoir (or to the other dot), at the dot's level
  (minimum plus zero-point energy); the attempt frequency is the dot's confinement frequency over 2π,
  clipped to 1 GHz–10 THz. `electrostatics.barrier_width_nm` supplies height and width. A smooth
  barrier transmits more than a rectangular one, so these are lower bounds within an order of magnitude.
- **Master equation** (`readout.master_equation.RateEquation`): states are the CI occupation window;
  lead transitions at Γ·f(ΔU/kT) and Γ·(1 − f), interdot at Γ_m·f(±ΔU/kT), so every pair of rates
  satisfies detailed balance and the stationary distribution is the Boltzmann one the CI model already
  uses (P0 test). Propagation along a sweep by the matrix exponential per sample; probability is
  carried between occupation windows by state label.
- **Latching is the lag**: when Γ · dwell ≪ 1 the distribution cannot follow, lines fade and shift with
  sweep direction; no separate latching rule. The Step 1 placeholder stays for the protocols that use
  it until they are re-rendered, and is marked superseded.
- **Not yet**: the RF tank response and catalog-tied fluctuator statistics are the next Step 4 PRs; the
  interdot rate is phenomenological (inelastic) until the Hubbard coherence is wired in.

## Consequences
- Easier: sweep-rate and direction dependence come out of the same model as the ground state; the
  reference for "how fast can this device be swept" is a derived rate, not a typed one.
- Harder: ~25×25 matrix exponential per sample (tens of microseconds); maps of 10⁴–10⁵ points take
  seconds, inside budget but the slowest protocols so far.
- Must now: P0 tests for detailed balance and for the latching threshold; protocols
  `stability-vs-sweep-rate`, `stability-direction`, `lead-rate-time-trace`.

## Alternatives considered
- Keep the coin-flip placeholder — rejected: no detailed balance, no direction dependence, no derived rate.
- Full Lindblad with coherences — deferred to the spin step; charge dynamics here are incoherent.
