# Step 4 — Readout physics

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
Dynamics, not just ground states: lead tunnel rates, rate equations, latching and sweep-rate dependence, the charge sensor and RF reflectometry response, and noise from fluctuator ensembles.

## Deliverables
- `readout`: WKB tunnel rates to leads from barrier heights; a master equation over the few-electron states; latching when sweeps outrun rates; sensor SET with Coulomb-blockade peaks and capacitive shifts; RF tank-circuit response
- fluctuator ensemble with a distribution of switching rates producing 1/f and telegraph noise with the right statistics

## P0 tests (CI: analytic limits and invariants)
- detailed balance in the master equation
- latching appears above a sweep-rate threshold predicted from the rates
- sensor peak position shifts by the capacitive coupling times the charge change

## P1 tests (CI, `literature`: reproduce an item of a `confirmed` reference entry; ADR 0007)
- a published charge-noise spectrum (amplitude at 1 Hz, slope) for the material class
- a published latching signature in a stability diagram

## Gallery protocols (rendered by CI, signed off per protocol by a human; ADR 0008)
- to be listed before the step starts; each carries a `budget` (ADR 0009)

## Exit criteria
- sweep-rate and direction change the measured diagram the way they do on a real device
- both P1 reproductions pass

## Open decisions
- how much of RF reflectometry to model (tank + sensor impedance) vs. a transfer function

## Status
NOT STARTED
