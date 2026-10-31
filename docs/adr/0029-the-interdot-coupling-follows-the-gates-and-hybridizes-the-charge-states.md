# 0029 — The interdot coupling follows the gates and hybridizes the charge states: t from half the slopes of ln Γ_m, capped at ħω₀/π; eigenstates of the CI window in the statics and a Pauli master equation over them in the dynamics

**Status:** Accepted
**Date:** 2026-10-31
**Step:** 8 follow-up (issue #58, the first of #58–#60; after ADR 0027)

## Context
In a simulated tune-up the owner found that "the barrier gates do not actually modulate the interdot
coupling AT ALL" (2026-10-31). On 0.3.0:

- **t on paper only.** `State.tunnel_coupling_ev` was t₀·exp((V_BM − V₀)/v_scale): it followed the barrier
  gate alone, had no ceiling (15 meV at B_M + 50 mV on the reference device) and entered nothing a sensor
  reads. It reached only the exchange and the transport placeholders.
- **A classical charge.** The sensor read the classical CI model and its rate equation. Opening B_M did
  not change an interdot line, apart from the interdot *rate* of ADR 0027 and a fixed cross-capacitance.
- **The owner's decisions (2026-10-31).** Hybridize the charge states in the statics and the dynamics.
  Reshape the dots with a compiled grid in #59. Follow with the exchange in #60. Release all three as
  0.4.0.

Four forces shaped the design:

- **Cost.** `State` is evaluated on grids of 10⁴–10⁵ points, and a session pays one matrix exponential
  per sample. A per-point diagonalisation must stay at microseconds.
- **Consistency.** ADR 0027 already compiles the slopes of ln Γ_m against every gate at a fixed tunnelling
  level. Γ_m ∝ T and t ∝ √T come from the same cut.
- **Limits.** At t → 0 the result must be the CI model, which the P0 tests and every earlier protocol
  rely on.
- **Reproducibility (ADR 0014).** Eigenvalues near 0.05 eV with kT = 13 µeV lose digits in a Boltzmann
  weight unless the matrix is conditioned.

## Decision
The interdot coupling follows every gate, and the CI window's charge states of one total charge
hybridize through it. The statics use the window's eigenstates, and the dynamics are a Pauli master
equation over them. Specifically:

- **t at the gates.** t(V) = t_max·σ(ln(t_ref/(t_max − t_ref)) + ½ s_m·ΔV):
  - t_ref is the compiled t at the operating point.
  - s_m are the slopes of ln Γ_m of ADR 0027.
  - t_max = ħω₀/π is compile's t = ħω₀√T/π at a transmission of one.
  - The form is `saturating_rates_hz`'s. `charge.hybrid.InterdotCoupling` carries it, and
    `Device.runtime_coupling()` builds it from fields the device file already has, so the device format
    is unchanged.
- **Hybridization.**
  - H_N = Σ U(N)|N⟩⟨N| − t Σ(|N₁+1, N₂−1⟩⟨N₁, N₂| + h.c.) on the window's states of each total charge,
    spinless with one level per dot.
  - The window is the CI window's (2w+1)² grid with its centre held at or above zero, and a negative
    occupation left out.
  - Each block is a chain of at most 2w+1 states, diagonalised in batches with its lowest level taken out
    first, so the eigenvalues keep their digits.
  - Charges are summed relative to the window's centre.
- **Statics.** With a `coupling`, `ConstantInteractionParams` gives the occupation as the diagonal of the
  thermal density matrix of the eigenstates. ⟨N⟩, the sensor's potential and current, and the lock-in
  follow from it. The ground state is the dominant charge state of the lowest eigenstate. Without a
  coupling the model is unchanged, bit for bit.
- **Dynamics.** `RateEquation.eigen_rate_matrix` gives the rates between eigenstates:
  - Lead i adds an electron from k to l at Γ_i·|⟨l|c†_i|k⟩|²·f((E_l − E_k)/kT) and removes it at the
    complement.
  - Within a total charge, k relaxes to l at Γ_m·g_kl·f(·), with g_kl = Σ_b w_bk·w_bl the product of the
    two eigenstates' weights on each bond b between charge states one moved electron apart. That is the
    sum of |⟨l|X_b|k⟩|² + |⟨l|Z_b|k⟩|², the bond's transfer and its dipole.
  - g is exactly Γ_m's classical rule at t = 0, and Γ_m for an isolated pair at every detuning.
  - Detailed balance holds in the eigenenergies.
  - When the gates move, the distribution is projected onto the new eigenstates, q_l = Σ_k p_k|⟨l|k⟩|²,
    with coherences dropped. `sweep` and `Session` share this.
- **What a sample reads.** A `sampled` session draws an eigenstate, and the sensor reads its ⟨N⟩. Its total
  charge is an integer and its split need not be. `Measurement.n_present` is that charge, and the
  contract checks the total. `State.tunnel_coupling_ev` is the gated t.
- **Release.** 0.4.0 (breaking). Model versions: `charge.hybrid` 1, `readout.master_equation` 2,
  `device` 5, `device.session` 5, `optimum.targets` 2.

## Consequences
- Easier:
  - Opening B_M widens an interdot line from the thermal step to ½(1 − ε/Ω·tanh(Ω/2kT)). On the reference
    device the 10–90 % width goes from 0.22 mV to 0.65 mV as t goes from 0.05 to 30 µeV
    (`interdot-line-vs-barrier`).
  - Closing B_M is what isolates a dot. With B_L closed alone, dot 2's lead loads dot 1 at
    Γ_R·t²(1/Δ + 1/Δ')². The electron-like and hole-like paths add, as in elastic cotunnelling, at about
    10⁶ Hz against Γ_L ≈ 10⁴ Hz. With both closed, a sweep latches and a parked point telegraphs
    (`session-closing-barrier`).
  - `optimum(device).coupling` inverts the same t a session uses.
- Harder:
  - **Validity of the projection.** The projection assumes a dwell long against ħ/Ω and dephasing in the
    eigenbasis. At |ε| ≲ t it shares the electron between the dots even when t is so small that the
    coherent transfer would take longer than the dwell (t ≲ ħ/dwell, sub-neV at microsecond dwells). On the
    reference device t stays above 0.05 µeV, and a sweep lands that close to a degeneracy rarely. There is
    no charge-dephasing rate and no Landau–Zener step.
  - **Plunger sensitivity.** t inherits the unscreened sensitivity of #48. The plungers move it as much as
    B_M, and two cells 18 mV apart on B_M differ by a factor of 30 in t at their interdot lines.
  - **Fixed shape and spin.** The dots' capacitances, lever arms and positions stay those of the
    operating point until #59. The model is spinless, so there is no √2 on singlet tunnelling (that is
    #60's sector) and no valley structure.
  - **Cost.** A grid point costs about 7 µs instead of 1.3 µs in the statics, and a sample a few small
    diagonalisations more.
  - **Exchange.** J still comes from the Hubbard model with the barrier-only exponential, so
    `State.exchange_ev` and `State.tunnel_coupling_ev` disagree until #60.
- Must now:
  - A consumer that labels charge states from `n_present` rounds it. An interdot line is no longer a step.
  - A test that isolates a dot closes its own barrier and B_M.
  - The contract's `sample_semantics` and `hold_moves_the_charge` checks compare the total charge.

## Alternatives considered
- **Hybridizing statics only.** Rejected by the owner. A session's charge would disagree with `State`
  near every interdot line.
- **A separate slope compile for t (central differences of `interdot_tunnel_coupling_ev` per gate).**
  Rejected. It costs 24 more barrier actions per compile, its level moves with the gates (the rule ADR
  0027 rejected), and ADR 0027 already shows the two agree on B_M.
- **Full density-matrix (Bloch–Redfield) dynamics.** Rejected for cost: 25² states per sample. A
  charge-dephasing rate and a Landau–Zener step would lift the projection's limit, and are left for when
  a protocol needs them.
- **Incoherent resonant tunnelling between localized states (a Lorentzian Γ(ε)).** Rejected. It cannot
  give the polarization width √(ε² + 4t²) at t ≫ kT, which is what a tune-up reads.
- **One 25×25 diagonalisation per point.** Rejected. Different total charges are degenerate across
  blocks, so the eigenvectors could mix them, and it costs more.
- **An integer `n_present` (the eigenstate's dominant charge state) beside a non-integer sensor input.**
  Rejected: `n_present` would no longer be what entered the sensor.
