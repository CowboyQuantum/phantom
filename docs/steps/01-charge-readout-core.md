# Step 1 — Charge and readout core

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
Migrate `qdt-sim`'s `DoubleDotModel` and noise models into `charge` and `readout` as the seed, strip every qdt dependency, make parameters explicit functions of gate voltage, prove the reproducibility contract, and stand up the validation tooling every later step uses: the reference library check, the gallery renderer with budgets, and the CI job that publishes it.

## Deliverables
- `charge`: `Dot`, `Barrier`, a two-dot constant-interaction model and a Hubbard model with t(V_barrier) exponential, U(dot size), V_ij; ground-state charge per gate point
- `readout`: `SensorSET` conductance vs. gate with capacitive coupling to dots; noise processes (white, telegraph, 1/f, fluctuator ensemble) with explicit time stepping
- Reproducibility: seed + `model_version` → bit-identical arrays; a cross-platform CI job
- `scripts/check_references.py` in `check_all.sh`; the first reference entries (`proposed`) for the quantities the P0/P1 rows below depend on (ADR 0007)
- `scripts/gallery.py`: runs every protocol, writes arrays, PNGs, `budget.json` and a mkdocs page to `build/gallery/`; verifies sign-off hashes; `gallery` dependency group with matplotlib (ADR 0008, 0009)
- CI: `gallery` job uploading the artifact on every PR; Pages deployment of the docs site from `main`
- Benchmark: ≥ 10⁴ gate points/s/core on a 2-dot Hubbard model, expressed as the `run` budget of `stability-coarse-2dot`

## P0 tests (CI: analytic limits and invariants)
- single-dot Coulomb diamonds: slopes equal the lever arms; addition energy equals U
- two-dot honeycomb: triple-point spacing from the CI capacitance matrix
- Hubbard limit t → 0 reproduces the CI model
- noise: PSD of the 1/f process matches its parameter within tolerance; telegraph rates match
- gallery arrays hash identically on both platforms; declared counts match measured counts

## P1 tests (CI, `literature`: reproduce an item of a `confirmed` reference entry; ADR 0007)
- none yet; candidates: honeycomb geometry and lever arms from a review of few-electron double dots; a 1/f charge-noise amplitude at 1 Hz for Si MOS or SiGe. Entries are proposed in this step and confirmed by a human before a test is written.

## Gallery protocols (rendered by CI, signed off per protocol by a human; ADR 0008)
- `turn-on-single-gate`: current vs. one plunger, others fixed; shows a threshold and saturation
- `pinch-off-barrier`: current vs. one barrier gate; shows exponential pinch-off
- `diamonds-single-dot`: current vs. (plunger, bias); Coulomb diamonds with slopes from the lever arms
- `stability-coarse-2dot`: charge state or sensor signal vs. (P1, P2), wide range; honeycomb
- `stability-fine-2dot`: same, zoomed on one triple-point pair; interdot line with tunnel-coupling broadening
- `sensor-coulomb-peak`: sensor conductance vs. sensor plunger; a peak and its slope
- `noise-time-trace`: sensor signal vs. time at a fixed point; telegraph and drift visible
- `noise-psd`: power spectral density of the trace; the 1/f slope and white floor
Each carries a `budget`; each is unsigned until the owner has reviewed it in the PR artifact.

## Exit criteria
- qdt-sim's committed scenario outputs are reproduced by the migrated core (one-off check, before Step 8 regenerates them)
- reproducibility job green on two platforms
- benchmark met, as the budget of `stability-coarse-2dot` within factor 1 locally
- every protocol above rendered and signed off

## Open decisions
- plan §8 #2 Hubbard vs. CI default

## Status
NOT STARTED
