# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** none in flight. Released: `v0.2.0` (on `f72daac`), `v0.3.0` (on `6920f19`) and `v0.4.0`
(on `ae8d5db`). Steps 1–7 are built; their open items (sign-offs, catalog entries, P1 rows, the solver) are on
#16, #17, #25, #28, #32, #34 and #37.
**Branch:** `issue-70/barrier-transition-width` (#70, with #69) from `main` at `ad0d3fa`; its pull request #71 is
open for review. It carries 0.5.0 (untagged, breaking: `CHANGELOG.md`).
**Last session:** 2026-11-01. The owner, running qdt on `v0.4.0`, saw the honeycomb turn into a superdot
"basically instantly" along the plungers at a place that did not move with B_M, and asked for the transition
to happen "over a larger space" (#70). The measurements behind it: the barrier top moved at −152 meV/V with
B_M and −128/−110 meV/V with the plungers (nothing screens them, #48), so f went from 0.2 to 0.9 within 15 mV
of the plungers; and the shape grid held the top past the merge node while the plungers' part kept growing,
so above B_M ≈ 0.35 V every map was the same and the double dot re-formed at low occupancy (also #69).

## Done this session
- **#70 ([ADR 0034](adr/0034-the-interdot-barrier-has-one-transmission.md)):** the interdot barrier has one
  transmission, Kemble's of its top ΔU(V), and t = (ħω₀/π)·√T, Γ_m = f_a·T and g = T all follow it.
  - `charge.barrier.InterdotBarrier`: ΔU(V) = ΔU_ref + Σ σ_g·δ_g, screened compiled slopes for every gate and,
    along the barrier gates, the shape grid's top continued linearly past each end
    (`ShapeGrid.quantity_along`). The shape itself is still held and flagged past the merge.
  - Placeholders (`calibration: none`): `CompileParams.interdot_barrier_screening` 0.4 (the barrier gate) and
    `interdot_cross_screening` 0.3 (every other gate, times the first); 1 and 1 is the bare potential.
    `Device.interdot_barrier_screening` records σ per gate; device format 6.
  - The compiled record is the same transmission: `tunnel_coupling_ev` = (ħω₀/π)√T_ref (57 µeV on the
    reference device, was 46 from the WKB action), the interdot `CompiledBarrier`'s rate and screened Kemble
    slopes, `barrier_v_scale_v` = 2/slope (23 mV). `State.interdot_barrier_top_ev`.
  - On the reference device: f 0.2 → 0.9 over 64 mV of B_M (was ~30) and 126 mV of the plunger diagonal (was
    15); each plunger a quarter of B_M's action; t by e every 23 mV of B_M (was 9); nothing merges within
    ±0.15 V of the programmed plungers at the programmed B_M; with B_M at 5 V one dot at every occupancy; below
    the compiled plungers B_M opens t and f again (#69).
  - Tests: `tests/test_interdot_barrier.py` (6 P0); the isolated-dot tests close B_M by 0.12 V now (t falls
    more slowly), the session's two-state test uses the pair's total rate (the path through dot 2 counts
    where P1 has slowed the left lead, #48). 216 tests.
  - Protocols `honeycomb-vs-barrier` and `interdot-line-vs-barrier` map at the new B_M offsets; the stored
    gallery device and `measured.json` re-recorded; 76 protocols within budget, all unsigned; contract OK.
- Model versions: `charge.barrier` 1, `charge.shape` 3, `charge.splitting` 2, `charge.hybrid` 3,
  `readout.rates` 2, `compile.device` 7, `device` 8, `device.session` 9; `charge` 5, `readout` 5.

## In flight
- The pull request #71 for #70: CI's ubuntu run and the Windows gallery comparison (new closed forms: a weighted
  PCHIP sum, a linear continuation, Kemble's logistic), and the preview at `/pr-<n>/`.
- The measured ledger was recorded on a slow machine (compile 0.63 s against 0.32 before, the old code
  itself 1.7× slower here): the session protocols read 5 ms/pt, within 3× of budget but not by much.

## Next action
The owner reviews #70's pull request (ADR 0034, the two changed protocols in its preview, the screening
defaults), merges it, decides whether to tag `v0.5.0` (ADR 0025) and when qdt (CowboyQuantum/qdt #172) moves
its pin. qdt's adapter must handle device format 6, `State.interdot_barrier_top_ev`, and the merged region
sitting at about B_M +0.11 V instead of +30 mV.

Open choices for later:
- extend the two screening factors to the lead and sensor barriers (#48): a plunger still moves a lead rate
  by e every 5 mV while it moves t by e every 90 mV;
- digitise a measured f(V_B) or t(V_B) curve to replace the hand-set screening by a P1 row;
- screening of C_m by the barrier gate, so a closed B_M gives a true grid (#48);
- compile's own last digits across platforms (ADR 0033).

## Open issues
#48 unscreened lead, sensor and geometric sensitivities; #49 per-device sensor-barrier tune; #50 gate 1/f
placeholder; #51 cp1252 console; #41 spin driven by voltages; #64 a dot that must lose two electrons loses
one; #65 recompiling at the sensor set point moves it. #69 is addressed by this branch (closed with #70).

Open blockers: none. Open decisions: plan §8 #1 (solver); release cadence against qdt's pins.

## Known gaps
- **Owner-dependent, never blocking:** sign-offs for all 76 protocols; catalog entries #1–#14; confirming the
  thirteen `proposed` reference entries.
- **Interdot barrier (ADR 0034):** the screening factors are set by hand to the owner's requirement; only the
  interdot barrier is screened; t reaches ħω₀/π (1.3 meV) where the barrier is gone; the continued top past
  the merge is a linear extrapolation, harmless because T is 1 there.
- **Splitting (ADR 0032):** a closed B_M stops at the geometric f_C ≈ 0.14; f(g) between g = 0.15 and 0.85 is
  an interpolation; the metallic-dot theory sits on the single-level hybridization.
- **Exchange (ADR 0031):** J jumps where the ground charge state changes; spinless, one valley.
- **Shape (ADR 0030):** the gates' changes add; C_m from the sphere formula; about 1 % pixel noise; β, the
  attempt rates and the plungers are not on the grid.
- **Hybridization (ADR 0029):** the projection assumes a dwell long against ħ/Ω; spinless.
- **Rates:** the lead rates are linear in the non-barrier gates with the bare potential's slopes (#48).
- **Sensor process:** first-order split; sequential noise; no joint dot–sensor energy (ADR 0028).
- **Placeholders** (`calibration: none`): as listed in the step files and ADRs 0019–0034.
- **Cross-platform (ADR 0014):** condition diagonalisations; sum charges relative to a centre; keep window
  truncation far below 1e-9 at rint ties; propagate generators with non-negative steps; form gaps from N − n_g.
- **Pages:** the source is `gh-pages` (root); per-PR previews are under `/pr-<n>/`.
