# 0022 — Blockade transport as sequential cycles with a spin-branch average; the funnel as two-level mixing; Elzerman shots as exponential draws; noise spectroscopy by filter functions on the draw's spectrum; gate fidelity as the entanglement fidelity of the Lindblad channel

**Status:** Accepted
**Date:** 2026-10-26
**Step:** 6

## Context
The second half of Step 6 needs transport with Pauli blockade, the spin funnel, single-shot readout
traces, charge-noise spectroscopy with the qubit, and a gate fidelity that reacts to the device's draws
(the step's exit criterion), all within rule 4's cost and rule 3's reproducibility.

## Decision
- **Bias triangles** (`spin.blockade`): sequential tunnelling through the electron and hole cycles at
  the (1,1)/(0,2) vertex, each cycle's current e / (¼ τ_S + ¾ τ_T) with τ the sum of inverse rates of its
  steps; lead steps carry Fermi factors from the constant-interaction energies, the interdot step is
  downhill-allowed with a thermal tail; the triplet branch passes only beyond E_V (the valley-excited
  (0,2) triplet) or through a mixing rate quenched by field as 1/(1 + (B/B_mix)²); reverse bias starts
  from the (0,2) singlet and has no blockade. Not a master equation over all states: ADR 0017's rate
  equation covers zero-bias dynamics; a finite-bias master equation with spin is a later cost decision.
- **Funnel** (`spin.mixing`): S–T₋ and S–T₀ as two-level anticrossings with time-averaged transfer
  ½ Δ²/(Δ² + δ²); the mixing Δ is a placeholder; ΔE_Z from the g draw; J from `spin.exchange`.
- **Elzerman** (`spin.elzerman`): shots as exponential draws in a fixed order (deterministic per seed);
  the analytic empty-probability curve is the ensemble limit and a P0 test.
- **Noise spectroscopy** (`spin.noise_spectroscopy`): W_N = exp(−∫ S_ω F_N/(2πf)² df) with the CPMG
  filter function evaluated numerically on a log grid; the spectrum is the draw's 1/f law with the
  amplitude defined over the sensor band and extended above it, so the qubit and the sensor share one
  amplitude; projection through dJ/dε and the lever arms.
- **Gate fidelity** (`spin.gates`): entanglement fidelity of the propagated Lindblad superoperator
  against the target unitary (which includes the nominal J and the g draw's Zeeman difference), averaged
  over the quasi-static detuning offset by quadrature; F_avg = (dF_e + 1)/(d + 1).
- **EDSR** (`spin.dynamics`): the drive strength ½ g μ_B (dB/dx) x_ac with a dipole displacement
  estimate; a placeholder until a micromagnet layout exists.

## Consequences
- Easier: the exit criterion holds — the SWAP fidelity varies with the noise draw and the readout
  visibility with the valley draw; the CPMG spectrum and the sensor PSD are one draw.
- Harder: the blockade picture has no co-tunnelling, no excited orbital lines and constant inelastic
  rates; the funnel is time-averaged; three-qubit cost stays open.
- Must now: P0 tests for each (blocked base, lifting by E_V and by mixing, no reverse blockade; funnel
  position; shot average; Ramsey, echo and CPMG limits and the white-noise invariance; perfect gate
  without noise, degrading with σ_ε, sweet spot); protocols `psb-bias-triangles`, `spin-funnel`,
  `elzerman-trace`, `charge-noise-spectroscopy-cpmg`, `charge-noise-spectroscopy-j-vs-ramsey`,
  `exchange-gate-fidelity`.

## Alternatives considered
- A finite-bias master equation with spin states for the triangles — correct but heavier; deferred.
- Monte Carlo over noise traces for CPMG — the filter-function integral is exact for Gaussian noise
  and deterministic.
