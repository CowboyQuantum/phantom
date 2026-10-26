# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** 1 (#16), 2 (#17), 3 (#25), 4 (#28) complete except the owner's sign-offs, catalog entries and P1
rows (and Step 4's lock-in dynamics on the master equation); 5 (#32, valley) built, its P1 row open; 6 (#34, spin)
in progress → `docs/steps/06-spin.md`. Milestones: "Steps 1-2", "Steps 3-4", "Steps 5-6: valley, spin".
**Branch:** `step-6/exchange-and-psb` (first Step 6 PR, part of #34). PRs #19–#33 merged.
**Last session:** 2026-10-26 — Step 4 closed (#31), Step 5 opened and merged (#33: `valley` layer, ADR 0020, four
protocols), Step 6 opened: `spin` layer (exchange, Lindblad dynamics, PSB), ADR 0021, eight protocols, ten P0 tests.

## Done
- `spin` (Step 6, first PR): `exchange` (J from the Hubbard two-electron block, the detuning axis through the
  plunger lever arms, the symmetric point, dJ/dε, σ_ε from gate noise, Gauss–Hermite quasi-static averaging),
  `dynamics` (rotating-frame Hamiltonians, `Lindblad` with T₁ and T_φ collapse operators and `expm` stepping,
  `SpinParams`/`SpinDraw`/`draw_spin`, Rabi, Ramsey and singlet-return formulas, T₂* = √2ħ/σ_E), `psb` (singlet
  and triplet charge steps, the window at the receiving dot's valley splitting, visibility). ADR 0021. Protocols
  `exchange-vs-detuning`, `exchange-vs-barrier`, `exchange-fingerprint`, `exchange-house-cells`, `esr-line`,
  `rabi-chevron`, `ramsey-fringes`, `psb-readout-window` (unsigned, within budget). Proposed reference entries
  `petta2005-singlet-triplet`, `reed2016-symmetric-exchange`, `yoneda2018-coherence-charge-noise` (from memory).
- `valley` (Step 5): Rice draw per dot with the alloy-disorder scalings, MOS deterministic part, interface-step
  hook; magnetospectroscopy kink, T₁ hot spot, pulsed-gate spectroscopy; ADR 0020; four protocols; two proposed
  entries (`paqueletwuetz2022-valley-alloy`, `yang2013-mos-valley`).
- Steps 1–4: `layout`, `electrostatics`, `disorder`, `charge`, `readout`; 43 protocols; tooling (`check_references`,
  `gallery.py` with budgets, sign-offs, `--record`; `model_graph.py --check`); CI with the ubuntu reference and the
  Windows comparison within 1e-9; per-PR gallery previews on Pages; `docs/MODELS.md` map, graph and ledger.

## In flight
- First Step 6 PR (this branch): everything under "Done" for `spin`, `docs/MODELS.md` (spin subgraph; graph and
  ledger refreshed), the step file, the indexes, `measured.json`. 99 P0 tests.

## Next action
- Second Step 6 PR: `psb-bias-triangles` (transport with blockade on the master equation or the Hubbard
  thermal state), `spin-funnel` (S–T₋ mixing against (ε, B) with the g draw), `elzerman-trace` (spin-dependent
  tunnel-out on the lead rate), `charge-noise-spectroscopy-cpmg` (filter functions on the draw's 1/f spectrum,
  must join the sensor `noise-psd`), `charge-noise-spectroscopy-j-vs-ramsey` (T₂* ∝ 1/(dJ/dε)); a two-qubit
  exchange-gate fidelity from the Lindblad propagator with the draw's noise (exit criterion); EDSR/micromagnet
  driving term. Then Step 7 (`compile`).

## Known gaps
- Owner-dependent, never blocking: sign-offs for all 51 protocols; catalog entries #1–#14; confirmation (or
  replacement) of the eight proposed reference entries; the papers for the P1 rows of Steps 2–6 (Step 5 needs
  the digitised distribution; Step 6 needs J(ε) and T₂*-vs-noise curves).
- Spin placeholders (`calibration: none`): g spread and dg/dV, hyperfine rms, T₁ and the Markovian T_φ, the PSB
  triplet T₁ at the read point. The quasi-static detuning noise takes the whole 1/f band of the noise draw
  (18 µeV rms on the default device), which is on the large side; the band split (quasi-static vs. Markovian)
  waits for the CPMG protocol. Ramsey traces of the two dots coincide because the hyperfine placeholder dominates
  and is the same per dot.
- Step 4 still owes the lock-in dynamics on the master equation. Amplifier noise is absent from the SNR picture.
  Noise draws and Step 3 generators run on the catalog draft's placeholder values. Valley: interface steps are a
  hook; the hot-spot background is a constant; pulsed-spectroscopy weights are fixed.
- A plunger shorted to a neighbouring barrier dissolves the dot; the fault protocol shows a plunger–sensor-gate
  short instead. No screening in the charge potentials until Step 7. Operating points carry `R_L: 0.7, R_R: 0.7`.
- Pages source is `gh-pages` (root); the docs root and per-PR previews under `/pr-<n>/` are served.
