# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** none in flight. Released: `v0.2.0` (on `f72daac`) and `v0.3.0` (on `6920f19`). 0.4.0 is
untagged and complete once this branch (#66) merges. It carries #58 (merged as #61), #59 (merged as #62)
and #60 (merged as #63). Steps 1–7 are built; their open items (sign-offs, catalog entries, P1 rows, the solver)
are on #16, #17, #25, #28, #32, #34 and #37.
**Branch:** `issue-66/interdot-splitting` (PR #67 against `main`). `main` is at `2f9f2d6`, the merge of #63.
**Last session:** 2026-10-31. The owner asked for the barrier gates to set the interdot coupling, the dots'
shape and the exchange (#58–#60). They then found the tune-up behaviour they wanted still missing: the
interdot width (anti-detuning) should grow with B_M until the dots form a superdot, and shrink to a grid as
it closes. Decisions taken that day:
- **t in the signal:** hybridize the charge states, in the statics and the dynamics.
- **Dot shape:** a grid compiled along the barrier gates and interpolated at runtime.
- **Splitting:** conductance-driven, as an effective interdot capacitance, with a smooth interpolation
  between Golden & Halperin's limits.
- **Release:** all of it folded into 0.4.0.

## Done this session
- **#58 (ADR 0029, #61):** t follows every gate and hybridizes the charge states.
- **#59 (ADR 0030, #62):** the barrier gates move and reshape the dots through a compiled grid. Also
  fixed: the hybridized window's seam, and a propagator that keeps every entry's relative digits.
- **#60 (ADR 0031, #63):** J of the valence electrons at the gated t and the shaped dots;
  `State.exchange_frequency_hz`; the Pauli readout at the gated t.
- **#66 ([ADR 0032](adr/0032-the-interdot-conductance-sets-the-triple-point-splitting.md), this
  branch):** the interdot conductance sets the triple-point splitting.
  - g is Kemble's transmission of the interdot barrier. Its top ΔU above the fixed tunnelling level
    follows every gate through compiled slopes and the shape grid's new `barrier_top_ev` column.
  - f = f_C + (F_max − f_C)·f_GH(g), where f_GH joins Golden & Halperin's weak and strong limits
    monotonically.
  - f enters the CI energies as an effective interdot capacitance ΔC_m·[[1, −1], [−1, 1]]; the induced
    charge keeps the geometric capacitances.
  - `charge.splitting`; `State.interdot_conductance`, `interdot_splitting`; device format 5.
  - On the reference device g goes from 0.002 at the operating point to 0.99 at B_M +50 mV, and f from
    0.16 to 0.94, reaching 0.98 at the merge.
  - New protocol `honeycomb-vs-barrier`. Tests: `tests/test_interdot_splitting.py`, 9 P0.
  - Reference entries `golden1996-peak-splitting`, `waugh1995-tunable-coupling` and
    `kemble1935-wkb-barrier` (`proposed`), plus the item `honeycomb-vs-interdot-coupling` of
    `vanderwiel2003-double-dot-transport`.
- **Cross-platform (ADR 0033, this branch).** `session-closing-barrier` failed the Windows comparison at
  1.1–2.2e-9. The Windows arrays (CI now keeps them on failure) showed compile's own numbers differing by
  about 1e-12 between platforms, amplified about a thousand-fold by a latched transition read on the
  sensor's flank. The runtime protocols now run on a device compiled once and stored under
  `docs/gallery/devices/` (`--record-devices`); a P0 test checks it against a fresh compile within 1e-9.
- **Release `v0.3.0`:** tagged and released as for `v0.2.0` (PR #57).
- **Totals:** 210 tests; 76 protocols within budget, all unsigned; contract OK; ledger re-recorded.

## In flight
- PR #67 (#66): waits for the owner's review of `honeycomb-vs-barrier`.

## Next action
The owner reviews #67. After it merges, 0.4.0 is complete. Its tag and qdt's pin are the
owner's decisions. The owner should also confirm or correct the four `proposed` reference entries and
items (the DOIs are from memory).

Open choices for later:
- digitise Waugh et al.'s f(g) for a P1 row;
- screening of C_m by the barrier gate, so a closed B_M gives a true grid (#48);
- the sensor current's cross-platform margin on some Windows runners (≈1.1e-9 in
  `session-closing-barrier`; macOS also misses 1e-9 in untouched sensor protocols).

## Open issues
#48 the unscreened sensitivity of the rates, t and J, and now the geometric splitting floor; #49
per-device sensor-barrier tune; #50 gate 1/f placeholder; #51 cp1252 console; #41 spin driven by voltages.
Candidates for new issues:
- the spin layer's sweet spot on the device's local sector (ADR 0031);
- the plungers and the sensor's β on the shape grid (ADR 0030);
- charge dephasing and Landau–Zener steps (ADR 0029);
- the sensor current's cross-platform margin.

Open blockers: none. Open decisions: plan §8 #1 (solver); release cadence against qdt's pins.

## Known gaps
- **Owner-dependent, never blocking:** sign-offs for all 76 protocols; catalog entries #1–#14;
  confirming the thirteen `proposed` reference entries.
- **Splitting (ADR 0032):**
  - A closed B_M stops at the geometric f_C ≈ 0.14 (small hexagons, not a grid).
  - f rises within about 40 mV of B_M.
  - f(g) between g = 0.15 and 0.85 is an interpolation.
  - The metallic-dot theory is applied on top of the single-level hybridization.
- **Exchange (ADR 0031):** J jumps where the ground charge state changes; spinless, one valley; the spin
  layer's Hubbard models keep their few-electron J.
- **Shape (ADR 0030):** the gates' changes add; C_m comes from the sphere formula; about 1 % pixel noise;
  β, the attempt rates and the plungers are not on the grid.
- **Hybridization (ADR 0029):** the projection assumes a dwell long against ħ/Ω, with no charge
  dephasing; spinless.
- **Rates and t:** linear in the non-barrier gates; no screening until the FEM backend exists.
- **Sensor process:** first-order split; sequential noise; no joint dot–sensor energy (ADR 0028).
- **Placeholders** (`calibration: none`): as listed in the step files and ADRs 0019–0032.
- **Cross-platform (ADR 0014):** condition diagonalisations; sum charges relative to a centre; keep window
  truncation far below 1e-9 at rint ties; propagate generators with non-negative steps; form gaps from
  N − n_g.
- **Pages:** the source is `gh-pages` (root); per-PR previews are under `/pr-<n>/`.
