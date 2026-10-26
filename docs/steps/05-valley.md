# Step 5 — Valley

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
Valley splitting per dot as a statistical draw for the material and dot size, with its consequences for Pauli spin blockade and relaxation.

## Deliverables
- `valley`: a distribution of valley splitting for a Ge profile and dot size (Ge-fluctuation-limited theory); per-dot draws; a hook for an interface-step model
- `charge` and later `spin` consume the draw: PSB window and lifting, relaxation hot spots

## P0 tests (CI: analytic limits and invariants)
- the distribution's mean and spread scale with Ge concentration and dot size as the theory states — `tests/test_valley.py` (σ ∝ 1/radius, √(x̄(1−x̄)), √(interface width); Rayleigh and large-ν limits of the Rice mean; draws reproduce the scalings and are seed-reproducible)
- the signatures sit where the splitting says: magnetospectroscopy kink at E_V/(g μ_B) with slopes ±μ_B, hot-spot peak at the same field, pulsed-spectroscopy excited line at E_V/α

## P1 tests (CI, `literature`: reproduce an item of a `confirmed` reference entry; ADR 0007)
- a published valley-splitting distribution for a stated heterostructure reproduces within tolerance — waits for `paqueletwuetz2022-valley-alloy` (proposed) to be confirmed and its figure digitised; then σ_ref and Δ_det are set from it
- the MOS splitting's linear field dependence — waits for `yang2013-mos-valley` (proposed)

## Gallery protocols (rendered by CI, signed off per protocol by a human; ADR 0008)
- `valley-splitting-draws` — histograms over 4000 dots for Si/SiGe and MOS, fraction below 50 µeV, mean against radius with the 1/radius law
- `magnetospectroscopy-transition` — the second-electron transition against field for three draws, kink at each draw's splitting
- `t1-hot-spot` — T₁ against field for three draws, Lorentzian dip where the Zeeman energy equals the splitting
- `pulsed-gate-spectroscopy` — occupation against plunger and pulse amplitude, the excited-state line E_V/α from the pulsed line
- each carries a `budget` (ADR 0009); all four are unsigned

## Exit criteria
- PSB readout quality varies between draws for physical reasons (needs Step 6's spin layer to consume the draw)
- the P1 reproduction passes (needs the confirmed entry)

## Open decisions
- whether to model interface steps explicitly — resolved for now by ADR 0020: not explicitly; a `step_factor` hook scales σ so a step model can plug in later

## Status
IN PROGRESS — `valley.splitting` and `valley.signatures` built (ADR 0020), six P0 tests, four protocols rendering and within budget; the P1 rows and the spin consumer are open.
