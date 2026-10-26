# Step 4 — Readout physics

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
Dynamics, not just ground states: lead tunnel rates, rate equations, latching and sweep-rate dependence, the charge sensor and RF reflectometry response, and noise from fluctuator ensembles.

## Deliverables
- `readout`: WKB tunnel rates to leads from barrier heights (done, `readout.rates`, with `electrostatics.barrier_width_nm` at the lead Fermi level); a master equation over the few-electron states with detailed balance (done, `readout.master_equation`); latching when sweeps outrun rates (emerges from the lag, done); sensor SET with Coulomb-blockade peaks and capacitive shifts (Step 1); RF tank-circuit response (to do)
- fluctuator ensemble with a distribution of switching rates producing 1/f and telegraph noise with the right statistics (Step 1 has the ensemble; tying its rates to the catalog is to do)
- the textbook layout gained reservoir accumulation gates R_L, R_R so lead barriers are the outer barrier gates (done)

## P0 tests (CI: analytic limits and invariants)
- detailed balance in the master equation
- latching appears above a sweep-rate threshold predicted from the rates
- sensor peak position shifts by the capacitive coupling times the charge change

## P1 tests (CI, `literature`: reproduce an item of a `confirmed` reference entry; ADR 0007)
- a published charge-noise spectrum (amplitude at 1 Hz, slope) for the material class
- a published latching signature in a stability diagram

## Gallery protocols (rendered by CI, signed off per protocol by a human; ADR 0008)
- `stability-vs-sweep-rate`: the layout-derived coarse map at 1 ms, 10 µs and 100 ns per point with WKB lead rates; equilibrium, lag, frozen (rendering)
- `stability-direction`: the same map scanned forward and backward at a lagging dwell, with the charge difference (rendering)
- `lead-rate-time-trace`: charge relaxation after a plunger step for three lead barriers; the slope is the derived lead rate (rendering)
- to come: `rf-resonance`, `rf-sensor-flank`, `single-shot-histogram`, `snr-vs-integration-time`; `switcher-histogram` and `drift-repeated-cuts` with catalog-tied rates

## Exit criteria
- sweep-rate and direction change the measured diagram the way they do on a real device
- both P1 reproductions pass

## Open decisions
- how much of RF reflectometry to model (tank + sensor impedance) vs. a transfer function

## Status
IN PROGRESS
