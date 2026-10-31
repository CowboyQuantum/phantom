# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** none in flight. Released: `v0.2.0` (on `f72daac`) and `v0.3.0` (on `6920f19`). 0.4.0 is
untagged: it carries #58 (merged as #61) and #59 (this branch), and #60 is planned for it. Steps 1–7 are
built; their open items (sign-offs, catalog entries, P1 rows, the solver) are on #16, #17, #25, #28, #32,
#34 and #37.
**Branch:** `issue-59/barrier-shaped-dots` (PR against `main`). `main` is at `663a7ab`, the merge of #61.
**Last session:** 2026-10-31. The owner found in a simulated tune-up that the barrier gates did not
change the interdot coupling at all. They asked for the barrier gates to set the interdot coupling, the
dots' positions and shape, and later the exchange. Decisions taken that day:
- **t in the signal:** hybridize the charge states, in the statics and the dynamics.
- **Dot shape:** a grid compiled along the barrier gates and interpolated at runtime.
- **Slicing:** #58 → #59 → #60, released together as 0.4.0.

## Done this session
- **#58 (ADR 0029, merged as #61):** the interdot coupling follows every gate and hybridizes the charge
  states. See `CHANGELOG.md`.
- **#59 ([ADR 0030](adr/0030-the-barrier-gates-reshape-the-dots-through-a-compiled-grid.md)): the barrier
  gates reshape the dots.**
  - **Compile.** The double dot is re-derived on 17 nodes along each of B_L, B_M and B_R (±0.2 V). Minima
    are refined between grid points (`refined_dot`), and the capacitance helper is
    `capacitance_matrix_af`. Each grid ends where the potential stops forming two dots, and that end is
    marked merged. `charge.shape.ShapeGrid`; `Device.shape`; device format 4.
  - **Runtime.** PCHIP per gate, with the changes added and held beyond each end.
    `ConstantInteractionParams.shape` gives C(V) and n_g = C(V)(α(V)·V − Δo(V))/E. `InterdotCoupling` and
    `GateRates` take their exponents along the barrier gates from the grid. Without a shape, everything
    is 0.3.0's.
  - **State** gains `dot_positions_nm`, `dot_radii_nm`, `addition_energies_ev` and `dots_merged`.
    `optimum.coupling` bisects t along B_M.
  - **On the reference device,** B_M +0.15 V draws the dots 18 nm together, C_m grows by 40 % and dot 1's
    addition energy falls by 1 meV. The dots merge at about +0.13 V from the tuned point.
  - **Tests** now follow one transition by its electron count, because barrier gates move the plunger
    lines (absolute charge). The isolated-dot configuration is B_L −0.11 and B_M −0.05.
  - New protocol `dot-shape-vs-barrier`. Tests: `tests/test_dot_shape.py`, 9 P0.
- **Release `v0.3.0`:** tagged and released as for `v0.2.0` (PR #57).
- **Totals:** 192 tests; 74 protocols within budget, all unsigned; contract OK; ledger re-recorded.

## In flight
- The PR for #59 waits for the owner's review of its preview at `/pr-<n>/`. Start with
  `dot-shape-vs-barrier`, then `interdot-line-vs-barrier` and `session-closing-barrier`.

## Next action
After #59 merges: #60. The exchange J comes from the singlet secular equation on the local (1,1)-like
sector of the CI window, at the gated t and at the U and V₁₂ of the shaped C(V). PSB reads at the gated t.
The exchange's bisection should also stop dividing by zero along B_M (a RuntimeWarning on `main` too; J
stays finite). Whether to tag 0.4.0, and when qdt moves its pin, are the owner's decisions.

## Open issues
#59 and #60 (this line of work); #48 the unscreened sensitivity of the rates and t (the plungers move t
as much as B_M); #49 per-device sensor-barrier tune; #50 gate 1/f placeholder; #51 cp1252 console; #41
spin driven by voltages. Open blockers: none. Open decisions: plan §8 #1 (solver); release cadence
against qdt's pins.

## Known gaps
- **Owner-dependent, never blocking:** sign-offs for all 74 protocols; catalog entries #1–#14;
  confirming the ten `proposed` reference entries.
- **Shape (ADR 0030):**
  - The gates' changes add, with no cross terms.
  - C_m comes from the sphere formula, so a merged double dot never becomes one dot; points past the
    merge are held and flagged.
  - The compile tier's curvature has pixel noise of about 1 %.
  - The sensor's β and the attempt rates stay at the operating point. The plungers are not on the grid.
  - A rate along a barrier gate levels off at about half the attempt rate.
- **Hybridization (ADR 0029):**
  - The projection between samples assumes a dwell long against ħ/Ω. There is no charge dephasing and
    no Landau–Zener step.
  - The model is spinless.
  - `State.exchange_ev` still uses the barrier-only exponential, until #60.
- **Rates and t:** linear in the gates other than the barrier gates, within about 20 % over ±10 mV. No
  screening until the FEM backend exists.
- **Sensor process:** first-order split within a dwell; sequential noise; no joint dot–sensor energy
  (ADR 0028).
- **Placeholders** (`calibration: none`): as listed in the step files and ADRs 0019–0030.
- **Cross-platform (ADR 0014):** avoid digit-cancelling quantities. Condition diagonalisations and sum
  charges relative to a centre.
- **Pages:** the source is `gh-pages` (root); per-PR previews are under `/pr-<n>/`.
