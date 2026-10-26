# Step 6 — Spin

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
Few-qubit spin dynamics under pulses: gate-dependent g, exchange J(ε, t_c), spin–orbit/EDSR or micromagnet driving, hyperfine, T₁, T₂*, and spin-blockade readout with fidelity.

## Deliverables
- `spin`: Hamiltonian from `charge` (J from Hubbard), Zeeman with g(V), driving terms, Lindblad with dephasing from charge noise through J and g sensitivities; Monte Carlo over noise for shot statistics
- readout fidelity and visibility from PSB plus the valley draw

## P0 tests (CI: analytic limits and invariants)
- blockade: the forward base is suppressed and lifted beyond E_V or at zero field, reverse bias unblocked; funnel at J = E_Z; Elzerman shots average to the analytic curve; filter function Ramsey, echo and CPMG limits and white-noise invariance; the SWAP gate is perfect without noise and degrades with σ_ε, with the sweet spot at the symmetric point — `tests/test_spin_blockade_noise.py`
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
- `psb-bias-triangles` — sequential-tunnelling triangles with the blocked base up to the valley draw, lifted at zero field, none in reverse
- `spin-funnel` — S–T₋ funnel along J(ε) = E_Z with the g draw's S–T₀ loss
- `elzerman-trace` — spin-dependent tunnel-out shots and the analytic read window
- `charge-noise-spectroscopy-cpmg`, `charge-noise-spectroscopy-j-vs-ramsey` — the qubit as a spectrometer of the draw's 1/f noise; T₂*(ε) ∝ 1/(dJ/dε)
- `exchange-gate-fidelity` — SWAP fidelity against operating detuning for the noiseless device and two noise draws
- each carries a `budget` (ADR 0009); all unsigned

## Exit criteria
- a two-qubit exchange gate's fidelity depends on the compiled device's noise and valley draws — met for the noise draw (`exchange-gate-fidelity`, P0 test) and the valley draw through the readout window (`psb-readout-window`); the compiled device itself is Step 7
- both P1 reproductions pass (need the confirmed entries)

## Open decisions
- 3-qubit cost; which semi-analytic gate models beyond that — the Lindblad step is 4ⁿ × 4ⁿ, so three spins cost 64 × 64 per expm, still cheap; left open until a three-dot layout exists

## Status
IN PROGRESS — `spin` complete for a pair (ADR 0021, 0022): exchange, dynamics, PSB, blockade transport, funnel, Elzerman, noise spectroscopy, gate fidelity; fifteen P0 tests; fourteen protocols rendering and within budget. Open: the two P1 rows (confirmed entries), the three-qubit cost decision, sign-offs.
