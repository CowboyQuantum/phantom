# 0031 — The exchange of the dots' valence electrons follows the gates: the singlet secular equation on the CI window's local sector at the gated coupling and the shaped dots; the Pauli readout at the gated coupling

**Status:** Accepted
**Date:** 2026-10-31
**Step:** 8 follow-up (issue #60, the last of #58–#60; after ADR 0029 and 0030)

## Context
After ADR 0029 and 0030 the barrier gates set the interdot coupling and the dots' shape, but not the
exchange. `State.exchange_ev` came from the two-site Hubbard model (ADR 0021). That model holds 0–2
electrons per dot and took t from the barrier gate alone, with no ceiling. A compiled device's dots hold
about 23 electrons each, so away from the Hubbard model's own symmetric point J described an unrelated
charge configuration. At the reference device's programmed point it read about 13 meV, whether t was
0.14 or 45 µeV. The Pauli spin blockade (PSB) readout used the compiled t, not the t at the gates.
The bisection behind J also divided zero by zero where t underflowed: a RuntimeWarning, with a finite
J.

The owner asked for the barrier gate to set the exchange frequency (2026-10-31).

Three forces shaped the design:

- **Precision.** J spans twelve decades. The secular equation solved relative to E(1,1) keeps its
  relative digits (ADR 0021), and the charge gaps have to keep theirs.
- **Consistency.** The CI window already carries the gated t (ADR 0029) and the shaped capacitances
  (ADR 0030).
- **Scope.** The spin layer's protocols and the coherence optimum use typed or few-electron Hubbard
  models at their symmetric point, which sits far from the operating point. There, a t extrapolated
  from the operating point's slopes would mean nothing.

## Decision
J at a gate point is the singlet–triplet splitting of the dots' valence electrons, from the
two-electron secular equation on the CI window's local sector. Specifically:

- **The sector.** The reference is the ground charge state (N₁, N₂) of the hybridized window. The
  partners are (N₁+1, N₂−1), Δ_L above, and (N₁−1, N₂+1), Δ_R above, at the CI energies at the gates.
  - Each gap is formed directly from N − n_g, as E(δᵀC⁻¹(N − n_g) + ½δᵀC⁻¹δ) with C = C(V), so it keeps
    its digits.
  - A partner with a negative occupation is left out (Δ = +inf).
  - t is the gated coupling.
  - `spin.exchange.local_exchange_ev`; the secular solver is shared (`singlet_exchange_ev`). Where
    t = 0 it returns J = −min(0, Δ_R, Δ_L) without dividing zero by zero.
- **Meaning.** This is the exchange of a qubit pair where both dots hold an odd number of electrons:
  the effective (1,1). Elsewhere the valence spins are paired, J describes no qubit, and it stays
  finite. Like the Hubbard model, the sector has spinless orbitals and one valley.
- **The consumer surface.**
  - `State.exchange_ev` is this J, and `State.exchange_frequency_hz` = J/h.
  - `Session.spin_readout` reads at the interdot coupling of the last point the session was given.
  - The few-electron Hubbard models, `runtime_hubbard` and the coherence optimum keep the barrier-gate
    exponential at their symmetric point.
- **Release.** 0.4.0. Model versions: `spin.exchange` 3, `device.session` 7.

## Consequences
- Easier:
  - The barrier gate sets the exchange frequency. On the reference device's (23, 28) cell, J/h goes
    from about 1 kHz with B_M 50 mV below the programmed point to 230 MHz at it and to tens of GHz as B_M
    opens. It follows 2t²(1/Δ_R + 1/Δ_L) while t ≪ Δ and falls below it as t nears the gaps.
  - Across an odd–odd cell, J rises toward both anticrossings (`exchange-frequency-vs-barrier`).
  - The readout's contrast closes with the barrier.
- Harder:
  - **Jumps.** J jumps where the ground charge state changes, because its reference sector changes with
    it. Inside a cell it is smooth.
  - **No spin structure.** There are no valleys or orbital excitations in the sector, so for multi-electron
    valence shells J is the single-orbital estimate.
  - **Two exchanges.** The Hubbard models in the spin layer still give the few-electron J of ADR 0021 at
    their symmetric point, so a device's `State.exchange_ev` and `optimum(device).coherence` describe
    different charge configurations.
  - **Implausibly large J.** t inherits the unscreened sensitivity of #48, so J reaches tens of GHz
    within tens of mV.
- Must now:
  - A consumer reads J and J/h per point.
  - A consumer reads the effective (1,1) from the parity of `n_ground`.

## Alternatives considered
- **Hubbard models placed at the device's valence electrons per point.** Rejected: a 16×16
  diagonalisation per point, the same physics as the secular equation, and it loses J's digits
  (ADR 0021).
- **J as NaN outside odd–odd cells.** Rejected: the contract requires `State` to be finite. The parity is
  available from `n_ground`.
- **The spin layer's sweet spot and gate fidelity on the device's local sector.** Deferred. It would move
  the coherence optimum to the operating point's cell, a change to a signed-off-class protocol family,
  and should be its own issue.
- **The compiled t in the readout.** Rejected: the owner asked for the readout to follow the barrier gates
  too.
