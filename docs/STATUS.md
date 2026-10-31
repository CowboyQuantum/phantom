# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** none in flight. Released: `v0.2.0` (on `f72daac`) and `v0.3.0` (on `6920f19`). 0.4.0 is
untagged and complete once this branch merges: #58 (merged as #61), #59 (merged as #62) and #60 (this
branch). Steps 1–7 are built; their open items (sign-offs, catalog entries, P1 rows, the solver) are on
#16, #17, #25, #28, #32, #34 and #37.
**Branch:** `issue-60/exchange-follows-barriers` (PR against `main`). `main` is at `204de96`, the merge of
#62.
**Last session:** 2026-10-31. In a simulated tune-up the owner found that the barrier gates did not change
the interdot coupling at all. They asked for the barrier gates to set the interdot coupling, the dots'
positions and shape, and the exchange frequency. Decisions taken that day:
- **t in the signal:** hybridize the charge states, in the statics and the dynamics.
- **Dot shape:** a grid compiled along the barrier gates and interpolated at runtime.
- **Slicing:** #58 → #59 → #60, released together as 0.4.0.

## Done this session
- **#58 (ADR 0029, #61):** t follows every gate and hybridizes the charge states. An interdot line widens
  with B_M, and closing B_L alone no longer isolates dot 1.
- **#59 (ADR 0030, #62):** the barrier gates move and reshape the dots through a grid compiled along each
  of them: positions, sizes, lever arms, capacitances, t and the rates; the merge is flagged.
  - It also fixed two cross-platform weaknesses that its Windows comparisons exposed.
  - The hybridized window is ±4 in the statics and follows t in the dynamics, which removes the seam at
    the middle of interdot segments.
  - The session propagates with `readout.master_equation.propagator`, which keeps every entry's
    relative digits.
- **#60 ([ADR 0031](adr/0031-the-exchange-of-the-valence-electrons-follows-the-gates.md), this
  branch):** the exchange of the dots' valence electrons follows the gates.
  - `spin.exchange.local_exchange_ev` solves the singlet secular equation on the CI window's local
    sector: (N₁, N₂) with (N₁±1, N₂∓1), at the gated t and the shaped C(V). Its gaps are formed from
    N − n_g directly.
  - `State.exchange_ev` uses it, and `State.exchange_frequency_hz` = J/h. `Session.spin_readout`
    reads at the last point's t.
  - The secular solver (`singlet_exchange_ev`) no longer divides zero by zero where t = 0.
  - On the reference device's (23, 28) cell, J/h runs from about 1 kHz to tens of GHz along B_M.
  - New protocol `exchange-frequency-vs-barrier`. Tests: `tests/test_local_exchange.py`, 5 P0.
- **Release `v0.3.0`:** tagged and released as for `v0.2.0` (PR #57).
- **Totals:** 199 tests; 75 protocols within budget, all unsigned; contract OK; ledger re-recorded.

## In flight
- The PR for #60 waits for the owner's review of its preview at `/pr-<n>/`. Start with
  `exchange-frequency-vs-barrier`.

## Next action
The owner reviews and merges #60. 0.4.0 is then complete: its tag (ADR 0025), and when qdt moves its
pin, are the owner's decisions. qdt's adapter must handle:
- device format 4;
- the eigenstate charge in `n_present` (round it to label charge states);
- the new `State` fields;
- J per point.

## Open issues
#48 the unscreened sensitivity of the rates, t and now J; #49 per-device sensor-barrier tune; #50 gate
1/f placeholder; #51 cp1252 console; #41 spin driven by voltages. Candidates for new issues:
- the spin layer's sweet spot and gate fidelity on the device's local sector (ADR 0031, deferred);
- the plungers and the sensor's β on the shape grid (ADR 0030);
- a charge-dephasing rate and a Landau–Zener step for the hybridized dynamics (ADR 0029).

Open blockers: none. Open decisions: plan §8 #1 (solver); release cadence against qdt's pins.

## Known gaps
- **Owner-dependent, never blocking:** sign-offs for all 75 protocols; catalog entries #1–#14;
  confirming the ten `proposed` reference entries.
- **Exchange (ADR 0031):**
  - J jumps where the ground charge state changes.
  - The sector has spinless orbitals and one valley.
  - The spin layer's Hubbard models and the coherence optimum keep the few-electron J at their
    symmetric point.
- **Shape (ADR 0030):**
  - The gates' changes add, with no cross terms.
  - C_m comes from the sphere formula, so a merged double dot is held and flagged, not modelled.
  - The compile tier's curvature has about 1 % pixel noise.
  - The sensor's β, the attempt rates and the plungers are not on the grid.
- **Hybridization (ADR 0029):** the projection assumes a dwell long against ħ/Ω, with no charge
  dephasing; the model is spinless.
- **Rates and t:** linear in the non-barrier gates; no screening until the FEM backend exists.
- **Sensor process:** first-order split; sequential noise; no joint dot–sensor energy (ADR 0028).
- **Placeholders** (`calibration: none`): as listed in the step files and ADRs 0019–0031.
- **Cross-platform (ADR 0014):**
  - Avoid digit-cancelling quantities. Condition diagonalisations and sum charges relative to a centre.
  - Keep window truncation far below 1e-9 at rint ties.
  - Propagate generators with non-negative steps. Form gaps from N − n_g.
- **Pages:** the source is `gh-pages` (root); per-PR previews are under `/pr-<n>/`.
