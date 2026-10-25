# Step 1 — Charge and readout core

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
Migrate `qdt-sim`'s `DoubleDotModel` and noise models into `charge` and `readout` as the seed, strip every qdt dependency, make parameters explicit functions of gate voltage, and prove the reproducibility contract.

## Deliverables
- `charge`: `Dot`, `Barrier`, a two-dot constant-interaction model and a Hubbard model with t(V_barrier) exponential, U(dot size), V_ij; ground-state charge per gate point
- `readout`: `SensorSET` conductance vs. gate with capacitive coupling to dots; noise processes (white, telegraph, 1/f, fluctuator ensemble) with explicit time stepping
- Reproducibility: seed + `model_version` → bit-identical arrays; a cross-platform CI job
- Benchmark: ≥ 10⁴ gate points/s/core on a 2-dot Hubbard model

## P0 tests (CI: analytic limits and invariants)
- single-dot Coulomb diamonds: slopes equal the lever arms; addition energy equals U
- two-dot honeycomb: triple-point spacing from the CI capacitance matrix
- Hubbard limit t → 0 reproduces the CI model
- noise: PSD of the 1/f process matches its parameter within tolerance; telegraph rates match

## P1 tests (CI, `literature`: published results, cited)
- none yet

## P2 (internal, never an exit criterion)
- comparisons with device statistics, when available

## Exit criteria
- qdt-sim's committed scenario outputs are reproduced by the migrated core (one-off check, before Step 8 regenerates them)
- reproducibility job green on two platforms
- benchmark met

## Open decisions
- plan §8 #2 Hubbard vs. CI default

## Status
NOT STARTED
