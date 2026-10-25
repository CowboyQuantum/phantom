# Step 6 — Spin

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
Few-qubit spin dynamics under pulses: gate-dependent g, exchange J(ε, t_c), spin–orbit/EDSR or micromagnet driving, hyperfine, T₁, T₂*, and spin-blockade readout with fidelity.

## Deliverables
- `spin`: Hamiltonian from `charge` (J from Hubbard), Zeeman with g(V), driving terms, Lindblad with dephasing from charge noise through J and g sensitivities; Monte Carlo over noise for shot statistics
- readout fidelity and visibility from PSB plus the valley draw

## P0 tests (CI: analytic limits and invariants)
- trace preservation and positivity; Rabi and Ramsey in the noiseless limit are analytic

## P1 tests (CI, `literature`: published results, cited)
- published J(ε) curve; a published T₂* vs. charge-noise relation

## P2 (internal, never an exit criterion)
- comparisons with device statistics, when available

## Exit criteria
- a two-qubit exchange gate's fidelity depends on the compiled device's noise and valley draws
- both P1 reproductions pass

## Open decisions
- 3-qubit cost; which semi-analytic gate models beyond that

## Status
NOT STARTED
