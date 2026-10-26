# STATUS

> Rewritten at the end of every session. Never appended. A fresh session continues from this file alone.

**Current step:** 1 (#16), 2 (#17), 3 (#25), 4 (#28) complete except the owner's sign-offs, catalog entries and P1
rows (and Step 4's lock-in dynamics on the master equation); 5 (#32, valley) built, its P1 row open; 6 (#34, spin)
built for a pair, its P1 rows open → `docs/steps/06-spin.md`. Milestones: "Steps 1-2", "Steps 3-4", "Steps 5-6:
valley, spin".
**Branch:** `step-6/blockade-transport-and-noise-spectroscopy` (second Step 6 PR, part of #34). PRs #19–#35 merged.
**Last session:** 2026-10-26 — Step 4 closed (#31), Step 5 built (#33), Step 6 first PR merged (#35: exchange,
Lindblad dynamics, PSB), second PR opened: blockade transport, funnel, Elzerman, noise spectroscopy, gate fidelity.

## Done
- `spin` (Step 6): `exchange` (J from the Hubbard two-electron block, detuning axis, symmetric point, dJ/dε, σ_ε,
  quadrature averaging), `dynamics` (rotating-frame Hamiltonians, `Lindblad`, spin draw, Rabi/Ramsey/singlet-return
  formulas, EDSR drive strength), `psb` (window at the valley splitting), `blockade` (sequential cycles at the
  (1,1)/(0,2) vertex with the triplet branch blocked below E_V and mixing quenched by field), `mixing` (funnel),
  `elzerman` (shots and the analytic read window), `noise_spectroscopy` (CPMG filter functions on the draw's 1/f
  spectrum), `gates` (SWAP fidelity as the Lindblad channel's entanglement fidelity, quadrature over σ_ε). ADR 0021,
  0022. Fourteen protocols (`exchange-vs-detuning`, `exchange-vs-barrier`, `exchange-fingerprint`,
  `exchange-house-cells`, `esr-line`, `rabi-chevron`, `ramsey-fringes`, `psb-readout-window`, `psb-bias-triangles`,
  `spin-funnel`, `elzerman-trace`, `charge-noise-spectroscopy-cpmg`, `charge-noise-spectroscopy-j-vs-ramsey`,
  `exchange-gate-fidelity`), all unsigned and within budget. Fifteen P0 tests. Proposed entries
  `petta2005-singlet-triplet`, `reed2016-symmetric-exchange`, `yoneda2018-coherence-charge-noise`.
- `valley` (Step 5): Rice draw per dot with the alloy-disorder scalings, MOS deterministic part, interface-step hook;
  magnetospectroscopy kink, T₁ hot spot, pulsed-gate spectroscopy; ADR 0020; four protocols; two proposed entries.
- Steps 1–4: `layout`, `electrostatics`, `disorder`, `charge`, `readout`; 43 protocols; tooling (`check_references`,
  `gallery.py` with budgets, sign-offs, `--record`; `model_graph.py --check`); CI with the ubuntu reference and the
  Windows comparison within 1e-9; per-PR gallery previews on Pages; `docs/MODELS.md` map, graph and ledger.

## In flight
- Second Step 6 PR (this branch): `spin.blockade`, `spin.mixing`, `spin.elzerman`, `spin.noise_spectroscopy`,
  `spin.gates`, EDSR helpers; six protocols; five P0 tests; ADR 0022; `docs/MODELS.md` (spin subgraph extended;
  graph and ledger refreshed); roster row `exchange-gate-fidelity`; step file; `measured.json`. 104 P0 tests.

## Next action
- Step 7 (`compile`): a compiled device from a layout, a scenario draw, a noise draw, valley and spin draws; the
  per-dot valley splittings attached from the harmonic radii; the lead rates at the read point feeding
  `spin.elzerman`; screening in the charge potentials; the runtime parameters every layer consumes. File the Step 7
  issue and a "Steps 7-8" milestone. Three-qubit cost (Step 6 open decision) waits for a three-dot layout.

## Known gaps
- Owner-dependent, never blocking: sign-offs for all 57 protocols; catalog entries #1–#14; confirmation (or
  replacement) of the eight proposed reference entries; the papers for the P1 rows of Steps 2–6.
- Spin placeholders (`calibration: none`): g spread and dg/dV, hyperfine rms, T₁, Markovian T_φ, the PSB triplet
  T₁, the blockade rates (lead, inelastic interdot, mixing and its field scale), the funnel mixing Δ, the Elzerman
  rates (parameters until the compiled device supplies the lead rates), the EDSR gradient and displacement. The
  quasi-static detuning noise takes the whole 1/f band of the draw (18 µeV rms on the default device); the CPMG
  spectrum extends the same amplitude to 1 GHz. The blockade picture has no co-tunnelling or excited orbital lines.
  Ramsey traces of the two dots coincide because the hyperfine placeholder dominates and is the same per dot.
- Step 4 still owes the lock-in dynamics on the master equation. Amplifier noise is absent from the SNR picture.
  Noise draws and Step 3 generators run on the catalog draft's placeholder values. Valley: interface steps are a
  hook; the hot-spot background is a constant; pulsed-spectroscopy weights are fixed.
- A plunger shorted to a neighbouring barrier dissolves the dot; the fault protocol shows a plunger–sensor-gate
  short instead. No screening in the charge potentials until Step 7. Operating points carry `R_L: 0.7, R_R: 0.7`.
- Pages source is `gh-pages` (root); the docs root and per-PR previews under `/pr-<n>/` are served.
