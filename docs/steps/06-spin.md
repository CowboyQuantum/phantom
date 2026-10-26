# Step 6 — Spin

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
Few-qubit spin dynamics under pulses: gate-dependent g, exchange J(ε, t_c), spin–orbit/EDSR or micromagnet driving, hyperfine, T₁, T₂*, and spin-blockade readout with fidelity.

## Deliverables
- `spin`: Hamiltonian from `charge` (J from Hubbard), Zeeman with g(V), driving terms, Lindblad with dephasing from charge noise through J and g sensitivities; Monte Carlo over noise for shot statistics
- readout fidelity and visibility from PSB plus the valley draw

## P0 tests (CI: analytic limits and invariants)
- trace preservation and positivity; Rabi and Ramsey in the noiseless limit are analytic — `tests/test_spin.py` (also: T₁ and T_φ rates, the quasi-static Gaussian envelope by quadrature, J against second-order perturbation theory at small t, J monotonic toward the anticrossing and exponential in the barrier, the sweet spot dJ/dε = 0, the detuning axis and noise projection, the PSB window at E_V and its draw dependence, the seeded g draw)

## P1 tests (CI, `literature`: reproduce an item of a `confirmed` reference entry; ADR 0007)
- published J(ε) curve — waits for `petta2005-singlet-triplet` or `reed2016-symmetric-exchange` (both proposed) to be confirmed and digitised
- a published T₂* vs. charge-noise relation — waits for `yoneda2018-coherence-charge-noise` (proposed)

## Gallery protocols (rendered by CI, signed off per protocol by a human; ADR 0008)
- `exchange-vs-detuning`, `exchange-vs-barrier` — J from the Hubbard sector with the second-order and exponential laws overlaid, dJ/dε
- `exchange-fingerprint` — singlet return against (ε, τ) with the draw's quasi-static detuning noise
- `exchange-house-cells` — singlet return against (P1, P2) at fixed τ with the barrier open: cells bounded by J contours and the honeycomb lines
- `esr-line`, `rabi-chevron`, `ramsey-fringes` — the g draw's two lines, the Lindblad chevron, Ramsey with the draw's T₂*
- `psb-readout-window` — singlet and triplet charge steps for two valley draws
- second PR: `psb-bias-triangles`, `spin-funnel`, `elzerman-trace`, `charge-noise-spectroscopy-cpmg`, `charge-noise-spectroscopy-j-vs-ramsey`
- each carries a `budget` (ADR 0009); all unsigned

## Exit criteria
- a two-qubit exchange gate's fidelity depends on the compiled device's noise and valley draws (gate fidelity waits for the second PR; the fingerprint visibility and the PSB window already do)
- both P1 reproductions pass (need the confirmed entries)

## Open decisions
- 3-qubit cost; which semi-analytic gate models beyond that — the Lindblad step is 4ⁿ × 4ⁿ, so three spins cost 64 × 64 per expm, still cheap; left open until a three-dot layout exists

## Status
IN PROGRESS — `spin.exchange`, `spin.dynamics`, `spin.psb` built (ADR 0021), ten P0 tests, eight protocols rendering and within budget; second PR, P1 rows and gate fidelity open.
