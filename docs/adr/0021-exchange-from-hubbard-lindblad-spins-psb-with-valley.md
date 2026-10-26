# 0021 — Exchange is read from the Hubbard two-electron sector; spins evolve under a Lindblad propagator with quasi-static noise averaged by quadrature; Pauli spin blockade's window is the valley draw

**Status:** Accepted
**Date:** 2026-10-26
**Step:** 6

## Context
The spin layer must tie the qubit's behaviour to the charge model, the noise draw and the valley draw
rather than to free parameters: the exchange J comes from the same tunnel coupling and charging
energies that shape the stability diagram, the dephasing from the same charge noise that the sensor
sees, and the readout window from the same valley splitting that `valley` draws per dot. The models
must stay cheap (rule 4) and bit-reproducible (rule 3).

## Decision
- **Exchange** (`spin.exchange`): J(V) = E(|↑↑⟩) − E(lowest two-electron singlet) from the six-state
  two-electron block of the Hubbard Hamiltonian at each gate point; no perturbative formula at runtime
  (the second-order 2t²/Δ_R + 2t²/Δ_L is a P0 limit and an overlay). The detuning axis is the on-site
  energy difference realised through the two plungers' lever-arm block at fixed sum; the symmetric point
  is ε_i = −U_i/2 − V₁₂. dJ/dε by a centred difference; detuning noise σ_ε from the gate-referred rms
  through (α₁ⱼ − α₂ⱼ).
- **Dynamics** (`spin.dynamics`): rotating-frame Hamiltonians ½δσ_z + ½Ωσ_x per spin plus J(S₁·S₂ − ¼);
  a Lindblad superoperator with √(1/T₁) σ₋ and √(1/2T_φ) σ_z per spin, exponentiated once per time step
  with `expm` and applied repeatedly (deterministic across platforms, unlike a non-Hermitian
  eigendecomposition). Slow noise (the 1/f band of the draw, the hyperfine field) is quasi-static: a
  Gaussian offset averaged by Gauss–Hermite quadrature, which gives the exp(−(τ/T₂*)²) envelope with
  T₂* = √2 ħ/σ_E exactly; no Monte Carlo at runtime, so the output is a deterministic function of the
  draw. The g-factors are a seeded draw (`SpinDraw`); dg/dV, the hyperfine rms, T₁ and the Markovian T_φ
  are placeholders (`calibration: none`) until the Step 6 entries are confirmed.
- **Pauli spin blockade** (`spin.psb`): the singlet's (0,2) occupation steps at the anticrossing
  (√2 t matrix element, thermal average with the upper branch), the triplet's at ε = E_V of the
  receiving dot; the contrast is the window; the visibility multiplies in the triplet's survival
  exp(−t_int/T₁). The sensor's own charge-readout error (Step 4) is a separate factor for the consumer.
- **Not yet**: three-qubit cost, EDSR/micromagnet driving terms, the spin-funnel, Elzerman readout,
  CPMG noise spectroscopy and the J-versus-Ramsey check (second Step 6 PR), J-gate fidelity.

## Consequences
- Easier: the fingerprint, house cells, J(ε), J(V_B), ESR lines, Rabi, Ramsey and the PSB window all
  follow from one device's draws; the fingerprint visibility and the Ramsey T₂* change with the noise
  draw, the PSB window with the valley draw, the ESR splitting with the g draw.
- Harder: quasi-static averaging ignores the spectral shape within the band (CPMG will need the
  spectrum itself); the Lindblad step cost grows as 4ⁿ for n spins.
- Must now: P0 tests for the perturbative limit, trace and positivity, Rabi, Ramsey, T₁, the PSB window
  and the draw dependence; protocols `exchange-vs-detuning`, `exchange-vs-barrier`, `exchange-fingerprint`,
  `exchange-house-cells`, `esr-line`, `rabi-chevron`, `ramsey-fringes`, `psb-readout-window`.

## Alternatives considered
- Heitler–London or Hund–Mulliken J from orbitals — a compile-tier solver; the Hubbard sector is exact
  for the two-site model phantom already uses and costs one 6×6 diagonalisation.
- Monte Carlo over noise traces at runtime — rejected for cost and reproducibility; quadrature over a
  quasi-static Gaussian is exact for the envelope and deterministic. Trace-based noise returns for CPMG.
