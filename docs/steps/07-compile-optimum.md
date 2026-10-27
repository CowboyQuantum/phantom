# Step 7 — Compile and optimum

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
The offline tier: FEM electrostatics for real stacks and overlapping gates through a wrapped external solver, Thomas–Fermi screening, and a `compile()` that produces a serialisable device; `optimum()` for every layer.

## Deliverables
- `compile`: layout + material + catalog draw + seed → `Device` (versioned, serialised); closed-form or FEM electrostatics selectable — `compile.compile_device` → `phantom.device.Device` (its own layer between `spin` and `optimum`, ADR 0023); closed form built, FEM a hook behind `[fem]` until plan §8 #1 is decided
- `optimum`: best operating point for sensor sensitivity, best lattice (virtual-gate matrix), target interdot coupling; with `spin`, achievable readout fidelity and coherence — `optimum.optimum(device)`: sensor flank, virtual-gate matrix, barrier voltage for a target t, PSB read point and Elzerman peak time, exchange sweet spot with T₂* and SWAP fidelity
- a solver wrapper behind the `[fem]` extra; never a solver written here

## P0 tests (CI: analytic limits and invariants)
- closed-form and FEM agree for a shallow single-layer layout within tolerance — waits for the solver (plan §8 #1); the FEM hook's refusal is tested
- a compiled device round-trips text-identically and its reloaded runtime objects reproduce the stability map bit for bit; compile is deterministic in the seed; the clean device equals the Step 2 derivation; wiring faults map the programmed point — `tests/test_compile_optimum.py`
- `optimum` is attained by a brute-force search on a small device — sensor flank against a dense scan, virtual gates decouple the dots exactly, the coupling target inverts the compiled barrier, the read point against a dense detuning scan, the sweet spot against a detuning scan of the SWAP fidelity (`tests/test_compile_optimum.py`)

## P1 tests (CI, `literature`: reproduce an item of a `confirmed` reference entry; ADR 0007)
- none

## Gallery protocols (rendered by CI, signed off per protocol by a human; ADR 0008)
- `compiled-device-roundtrip` — the coarse map from the compiled and from the reloaded device, and their (zero) difference
- `virtual-gates-map` — the coarse map in raw and in virtual plunger coordinates
- `optimum-sensor-flank` — the sensor peak and the flank the optimum chooses
- `optimum-coupling-target` — t against the barrier gate at the compile tier and from the runtime exponential; the optimum's choice for three targets beside the compile-tier crossing
- `optimum-readout-point` — the blockade contrast with the chosen read point; the Elzerman empty probability with the chosen read time
- `optimum-sweet-spot` — SWAP infidelity through the chosen sweet spot, and J
- `compiled-device-draws` — twelve devices from one layout: coupling, lead rates, valley draws, and what the optimum reaches; devices that do not form are recorded
- each carries a `budget` (ADR 0009); all unsigned

## Exit criteria
- a compiled device round-trips through serialisation and reproduces its diagrams — met (P0 and `compiled-device-roundtrip`)
- `optimum` exists for every layer present — met for sensor, lattice, coupling, readout and coherence; the layers' protocols are not yet all re-rendered from a compiled device

## Open decisions
- plan §8 #1 solver choice — still open; the `fem` backend is a hook
- #3 serialisation format — decided (ADR 0023): canonical JSON, basis recomputed, format and model versions recorded

## Status
IN PROGRESS — `compile`, `device` and `optimum` built (ADR 0023), twelve P0 tests, seven protocols rendering; the Step 4 dynamics protocols use the compiled device's WKB action. Open: the solver choice and its agreement test (plan §8 #1), the owner's sign-offs.
