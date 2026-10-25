# Step 6 — Spin

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
Few-qubit spin dynamics under pulses: gate-dependent g, exchange J(ε, t_c), spin–orbit/EDSR or micromagnet driving, hyperfine, T₁, T₂*, and spin-blockade readout with fidelity.

## Deliverables
- `spin`: Hamiltonian from `charge` (J from Hubbard), Zeeman with g(V), driving terms, Lindblad with dephasing from charge noise through J and g sensitivities; Monte Carlo over noise for shot statistics
- readout fidelity and visibility from PSB plus the valley draw

## P0 tests (CI: analytic limits and invariants)
- trace preservation and positivity; Rabi and Ramsey in the noiseless limit are analytic

## P1 tests (CI, `literature`: reproduce an item of a `confirmed` reference entry; ADR 0007)
- published J(ε) curve; a published T₂* vs. charge-noise relation

## Gallery protocols (rendered by CI, signed off per protocol by a human; ADR 0008)
- to be listed before the step starts; each carries a `budget` (ADR 0009)

## Exit criteria
- a two-qubit exchange gate's fidelity depends on the compiled device's noise and valley draws
- both P1 reproductions pass

## Open decisions
- 3-qubit cost; which semi-analytic gate models beyond that

## Status
NOT STARTED
