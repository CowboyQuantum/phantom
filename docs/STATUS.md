# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** none in flight. Released: `v0.2.0` (on `f72daac`), `v0.3.0` (on `6920f19`) and `v0.4.0`
(on `ae8d5db`). 0.5.0 (#70, ADR 0034) is merged on `main` at `cf9e1b2`, untagged. Steps 1–7 are built; their
open items (sign-offs, catalog entries, P1 rows, the solver) are on #16, #17, #25, #28, #32, #34 and #37.
**Branch:** `issue-72/puddles` (#72) from `main` at `cf9e1b2`; its pull request #73 is open for review. It carries
0.6.0 (untagged, breaking: `CHANGELOG.md`).
**Last session:** 2026-11-01, two pieces of work for the owner.
- **#70 (merged as #71, ADR 0034, 0.5.0):** the interdot barrier has one transmission, Kemble's of its top,
  shared by t, Γ_m and g; the gates' action on the top is screened by two placeholders (0.4 for B_M, times
  0.3 for every other gate); the top continues past the shape grid's merge node. The honeycomb-to-superdot
  transition spans 64 mV of B_M (was 30) and 126 mV of the plunger diagonal (was 15); a far-open B_M
  leaves one dot at every occupancy (#69 closed with it).
- **#72 (this branch, [ADR 0035](adr/0035-the-dots-grow-with-their-occupancy.md), 0.6.0):** the owner then
  asked for the dots' size against their occupancy, so that the triple-point width grows first at the
  highest cells and moves toward (1,1) as B_M opens. Each dot is now a classical puddle with R ∝ N^(1/3)
  anchored at the compiled dot (`charge.puddle`): its capacitances scale as √(n/N_ref) with the
  gate-induced charge (E_add falls as N^(−1/3), f_C grows, the classical N ∝ φ^(3/2)), and the interdot
  barrier is seen from the puddles' edges on the compiled cut (zero shift at the anchor, tens of meV
  higher at the first electrons). `State.puddle_radii_nm`; `CompileParams.puddles`; device format 7.

## Done this session
- ADR 0034 and ADR 0035; CHANGELOG 0.5.0 and 0.6.0; CONTRACT, MODELS (nodes IDB and PUD, protocol rows).
- `charge.barrier.InterdotBarrier` (#70) and `charge.puddle.PuddleGrowth` (#72); `ShapeGrid.quantity_along`;
  the scale threaded through `InterdotSplitting`, `InterdotCoupling`, `GateRates`, `RateEquation`.
- Compile stores the screening per gate, the cut between the dots, the turning points and the anchor
  occupancy; the runtime finds exactly the compiled numbers at the operating point.
- On the reference device at the programmed gates, from (1,1) to (28,35): puddle radii 5/12 → 17/40 nm,
  E_add 32/28 → 6/6 meV, f_C 0.06 → 0.16, the barrier top 42 meV above the level → below it. The superdot
  front: (26,29) at B_M 0.3 V, (19,19) at 0.5 V, (10,7) at 0.8 V, the first electrons at 1.0 V.
- Tests: `tests/test_interdot_barrier.py` (6) and `tests/test_puddles.py` (6); 222 tests. Protocols:
  `honeycomb-vs-barrier` and `interdot-line-vs-barrier` at the new B_M offsets, new `honeycomb-vs-occupancy`;
  77 protocols within budget, all unsigned; the stored gallery device and `measured.json` re-recorded.
- Reference entry `shikin1991-parabolic-dot-classical-density` (`proposed`; the citation was checked against
  the journal's listing, the density and the N ∝ R³ relation must be confirmed from the paper).

## In flight
- The pull request #73 for #72: CI's ubuntu run and the Windows gallery comparison (new closed forms: a square
  root of the induced charge, scaled capacitances, a linear interpolation on the compiled cut), and the
  preview at `/pr-<n>/`.

## Next action
The owner reviews #72's pull request (ADR 0035, `honeycomb-vs-occupancy` in its preview, the anchor at the
compiled dot), merges it, decides whether to tag `v0.6.0` (or `v0.5.0` first, on `cf9e1b2`; ADR 0025) and
when qdt (CowboyQuantum/qdt #172) moves its pin. qdt's adapter must handle device format 7,
`State.interdot_barrier_top_ev` and `puddle_radii_nm`, the merged region at about B_M +0.11 V from the
tuned point, and charge maps whose period shrinks toward high occupancy.

Open choices for later:
- anchor the puddle at one electron (R(1) = ℓ) once the solver gives the radius at every N; today's
  anchor reads the one-electron harmonic radius as the puddle of N_ref ≈ 23 electrons, which makes the
  (1,1) addition energy about 32 meV;
- extend the two screening factors to the lead and sensor barriers (#48);
- digitise a measured f(V_B), t(V_B) or E_add(N) curve for a P1 row;
- screening of C_m by the barrier gate (#48); compile's last digits across platforms (ADR 0033).

## Open issues
#48 unscreened lead, sensor and geometric sensitivities; #49 per-device sensor-barrier tune; #50 gate 1/f
placeholder; #51 cp1252 console; #41 spin driven by voltages; #64 a dot that must lose two electrons loses
one; #65 recompiling at the sensor set point moves it.

Open blockers: none. Open decisions: plan §8 #1 (solver); release cadence against qdt's pins.

## Known gaps
- **Owner-dependent, never blocking:** sign-offs for all 77 protocols; catalog entries #1–#14; confirming the
  fourteen `proposed` reference entries.
- **Puddles (ADR 0035):** classical N^(1/3), no shells; the compiled cut held where the shape grid moves
  the dots; the edge shift on the bare profile, not a self-consistent one; the anchor's one-electron
  radius read as the N_ref puddle; the induced charge is s·n up to the mutual capacitance's share (2 %).
- **Interdot barrier (ADR 0034):** the screening factors are set by hand; only the interdot barrier is
  screened; t reaches ħω₀/π where the barrier is gone; the continued top past the merge is linear.
- **Splitting (ADR 0032):** f(g) between g = 0.15 and 0.85 is an interpolation; the metallic-dot theory
  sits on the single-level hybridization.
- **Exchange (ADR 0031):** J jumps where the ground charge state changes; spinless, one valley.
- **Shape (ADR 0030):** the gates' changes add; C_m from the sphere formula; about 1 % pixel noise; β, the
  attempt rates and the plungers are not on the grid.
- **Hybridization (ADR 0029):** the projection assumes a dwell long against ħ/Ω; spinless.
- **Rates:** the lead rates are linear in the non-barrier gates with the bare potential's slopes (#48).
- **Sensor process:** first-order split; sequential noise; no joint dot–sensor energy (ADR 0028).
- **Placeholders** (`calibration: none`): as listed in the step files and ADRs 0019–0035.
- **Cross-platform (ADR 0014):** condition diagonalisations; sum charges relative to a centre; keep window
  truncation far below 1e-9 at rint ties; propagate generators with non-negative steps; form gaps from N − n_g.
- **Pages:** the source is `gh-pages` (root); per-PR previews are under `/pr-<n>/`.
