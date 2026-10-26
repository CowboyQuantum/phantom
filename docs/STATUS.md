# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** 1 (#16), 2 (#17), 3 (#25), 4 (#28) complete except the owner's sign-offs, catalog entries and P1
rows (and Step 4's lock-in dynamics on the master equation); 5 (#32, valley) in progress → `docs/steps/05-valley.md`.
Milestones: "Steps 1-2", "Steps 3-4", "Steps 5-6: valley, spin".
**Branch:** `step-5/valley-draws` (first Step 5 PR, part of #32). PRs #19–#31 merged.
**Last session:** 2026-10-26 — Step 4 closed with its third PR (#31: noise draws, switcher and drift protocols,
plunger-vs-barrier on the master equation); Step 5 opened: `valley` layer, ADR 0020, four protocols, six P0 tests.

## Done
- `valley` (Step 5): `splitting.ValleyParams` with the alloy-disorder scalings (σ ∝ 1/radius, √(x̄(1−x̄)),
  √(interface width), a `step_factor` hook), a deterministic part that for MOS grows with the vertical field,
  `draw_valley_splitting_uev` (Rice draw per dot, seeded), `rice_mean`; `signatures` (second-electron
  magnetospectroscopy shift with the kink at g μ_B B = E_V, Lorentzian T₁ hot spot, pulsed-gate spectroscopy
  occupation with the excited-state line at E_V/α, its weight vanishing smoothly at zero amplitude). ADR 0020.
  Protocols `valley-splitting-draws`, `magnetospectroscopy-transition`, `t1-hot-spot`, `pulsed-gate-spectroscopy`
  (all unsigned, within budget). Proposed reference entries `paqueletwuetz2022-valley-alloy`, `yang2013-mos-valley`
  (written from memory; the owner confirms DOI, authors and figures). The reference values of both material
  classes are placeholders (`calibration: none`) until those entries are confirmed.
- Steps 1–4: `layout`, `electrostatics`, `disorder` (charges, geometry faults, wiring faults, trapping, noise draws,
  scenario draw), `charge` (CI, Hubbard, derived from a layout), `readout` (SET sensor, noise processes, lock-in
  channel, WKB rates, master equation, RF tank, single shot, transport placeholders); 39 protocols rendering; 89 P0
  tests; `check_references.py`, `gallery.py` with budgets, sign-off hashes, `--record`; `model_graph.py --check`;
  CI with the ubuntu reference and the Windows comparison within 1e-9 (ADR 0010, 0014); per-PR gallery previews
  on Pages (ADR 0012); `docs/MODELS.md` measurement-flow diagram, module graph and runtime ledger (ADR 0011).

## In flight
- First Step 5 PR (this branch): everything under "Done" for `valley`, plus `docs/MODELS.md` (valley subgraph in the
  hand-drawn map; generated graph and ledger refreshed), the step file, the ADR index, `measured.json`.

## Next action
- After this PR merges: Step 6 (`spin`): Zeeman and g-factor draws, Pauli spin blockade consuming the valley draw
  (the Step 5 exit criterion "PSB readout quality varies between draws" is tested there), exchange J(ε, V_B),
  the house-cell and fingerprint protocols the owner asked for (axes P1, P2 at the evolve coordinate with J held
  high; detuning vs J), relaxation background for the hot spot. File the Step 6 issue under "Steps 5-6: valley, spin".

## Known gaps
- Owner-dependent, never blocking: sign-offs for all 43 protocols; catalog entries #1–#14; confirmation (or
  replacement) of the five proposed reference entries (van der Wiel 2003, Connors 2019, Davies 1995, Paquelet Wuetz
  2022, Yang 2013); the papers for the P1 rows of Steps 2, 3, 4 and 5 (Step 5's P1 needs the digitised distribution
  figure of the confirmed entry; then σ_ref and Δ_det are fitted and the placeholders retired).
- Step 4 still owes the lock-in dynamics on the master equation (`readout.loading.LeadRateLandscape` and
  `latched_occupation` remain for it; a master-equation step per half-cycle costs minutes per map). Amplifier noise
  is absent from the SNR picture. Noise draws and Step 3 generators run on the catalog draft's placeholder values.
- Valley: interface steps are a hook (`step_factor`), not a model; the hot-spot background rate is a constant until
  `spin` supplies its field dependence; the pulsed-spectroscopy weights are fixed (no loading-rate model behind them);
  `compile` (Step 7) will attach the per-dot draws to a compiled device using the harmonic radii from `electrostatics`.
- A plunger shorted to a neighbouring barrier dissolves the dot; the fault protocol shows a plunger–sensor-gate
  short instead. No screening in the charge potentials until Step 7. Operating points carry `R_L: 0.7, R_R: 0.7`.
- Pages source is `gh-pages` (root); the docs root and per-PR previews under `/pr-<n>/` are served.
