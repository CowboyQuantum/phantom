# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** none in flight. Released: `v0.2.0` (on `f72daac`) and `v0.3.0` (on `6920f19`, #47, #53,
#54). 0.4.0 is untagged: it carries #58, and #59 and #60 are planned for it. Steps 1–7 are built; their open
items (sign-offs, catalog entries, P1 rows, the solver) are on #16, #17, #25, #28, #32, #34 and #37.
**Branch:** `issue-58/interdot-hybridization` (PR against `main`). `main` is at `746ab04`, the merge of #57.
**Last session:** 2026-10-31. The owner found in a simulated tune-up that the barrier gates did not
change the interdot coupling at all, and asked for them to set the interdot coupling, the dots' positions
and shape, and later the exchange. Decisions taken that day (via questions):
- **t in the signal:** hybridize the charge states, in the statics and the dynamics.
- **Dot shape:** compile a grid along the barrier gates and interpolate at runtime (#59).
- **Slicing:** three issues in order, #58 → #59 → #60, released together as 0.4.0.

## Done this session
- **#58 ([ADR 0029](adr/0029-the-interdot-coupling-follows-the-gates-and-hybridizes-the-charge-states.md)):
  the interdot coupling follows the gates and hybridizes the charge states.**
  - `charge.hybrid`: `InterdotCoupling` (the compiled t moved by ½ the slopes of ln Γ_m, capped at
    ħω₀/π) and the CI window's eigenstates, one chain per total charge, diagonalised in batches with each
    block's lowest level taken out.
  - `ConstantInteractionParams.coupling`; `Device.runtime_coupling()`; `RateEquation.eigen_rate_matrix`,
    `project` and `carried_states`: a Pauli master equation over the eigenstates, shared by `sweep` and
    `Session`. Without a coupling, every output is 0.3.0's, bit for bit.
  - `Measurement.n_present` in `sampled` mode is an eigenstate's ⟨N⟩. The contract checks its integer
    total, and the hold check compares the total charge.
  - Closing B_L alone no longer isolates dot 1: dot 2's lead loads it at Γ_R t²(1/Δ + 1/Δ')². The tests
    and `session-closing-barrier` now isolate it with B_L and B_M closed.
  - New protocol `interdot-line-vs-barrier`. `optimum.coupling` inverts the gated t.
  - Tests: `tests/test_interdot_hybridization.py`, 8 P0, plus the cotunnelling bypass in
    `test_dot_rates.py`.
- **Release `v0.3.0`:** tagged on `6920f19` after CI passed there; GitHub Release made as for `v0.2.0`;
  noted on #53 and #54. PR #57 recorded it.
- **Totals:** 183 tests; 73 protocols within budget, all unsigned; contract OK; ledger re-recorded.

## In flight
- The PR for #58 waits for the owner's review of its preview at `/pr-<n>/`. Start with
  `interdot-line-vs-barrier` and `session-closing-barrier`.

## Next action
After #58 merges: #59. Compile re-solves the dots on a grid along B_L, B_M and B_R (positions, C_dot,
C_gate, the actions behind t and the rates), the runtime interpolates them, and the dots merge as B_M
opens. Device format 4. Then #60: J from the singlet secular equation on the local (1,1)-like sector of
the CI window at the gated t, U and V₁₂, and PSB at the gated t. Whether to tag 0.4.0, and when qdt
moves its pin, are the owner's decisions.

## Open issues
#58, #59, #60 (this line of work); #48 rates' and now t's unscreened sensitivity (the plungers move t as
much as B_M); #49 per-device sensor-barrier tune; #50 gate 1/f placeholder; #51 cp1252 console; #41 spin
driven by voltages. Open blockers: none. Open decisions: plan §8 #1 (solver); release cadence against
qdt's pins.

## Known gaps
- **Owner-dependent, never blocking:** sign-offs for all 73 protocols; catalog entries #1–#14;
  confirming the ten `proposed` reference entries.
- **Hybridization (ADR 0029):**
  - The projection between samples assumes a dwell long against ħ/Ω. At |ε| ≲ t it shares the electron
    even when t is below ħ/dwell. There is no charge dephasing rate and no Landau–Zener step.
  - The model is spinless, with one level per dot.
  - `State.exchange_ev` still uses the barrier-only exponential, until #60.
- **Rates and t:** linearised exponents, within about 20 % of the compile tier over ±10 mV. No screening
  until the FEM backend exists.
- **Sensor process:** first-order split within a dwell; sequential noise; no joint dot–sensor energy
  (ADR 0028).
- **Placeholders** (`calibration: none`): as listed in the step files and ADRs 0019–0029.
- **Cross-platform (ADR 0014):** avoid digit-cancelling quantities. Condition diagonalisations (take
  each block's lowest level out first) and sum charges relative to a centre.
- **Pages:** the source is `gh-pages` (root); per-PR previews are under `/pr-<n>/`.
