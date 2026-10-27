# Step 8 — qdt cutover

> Fixed headings. Agents read this file for the current step; the narrative is in `DEVELOPMENT_PLAN.md`.

## Goal
qdt's `qdt-sim` becomes an adapter that pins a phantom release and contains no physics; qdt's golden files are regenerated once; the seam is enforced on both sides.

## Deliverables
- in phantom: a release tag and changelog; the consumer contract (plan §4) documented as the public API — `phantom.api` (ADR 0024), `docs/CONTRACT.md`, `CHANGELOG.md`, version 0.2.0, the surface snapshot `docs/api/surface.json` checked in `check_all.sh`; the runtime surface `Device.at` → `State`, `Device.session` → `Session` (`sense` with bias and samples, `hold`, `wait`, `spin_readout`; a noise seed apart from the device seed; one sample's charge stated), `Device.truth` (charge lattice, sensor peak and wing) and `Device.knowable`; the release is `v0.2.0` at `f72daac` (who tags: ADR 0025)
- in qdt (coordinated issue): `qdt-sim` maps `StandConfig`/`ChannelMap` onto `phantom.Device`, implements the protocols over it, pins the release; the physics modules deleted; `make_sim_fixtures.py` regenerated; the vocabulary and import checks updated

## P0 tests (CI: analytic limits and invariants)
- phantom's API surface is covered by a contract test that qdt's adapter also runs — `phantom.contract` (`python -m phantom.contract [device.json ...]`), fifteen checks, run by `tests/test_contract.py` on three draws
- the session against closed forms: Boltzmann draws without dynamics, trajectory switching as the two-state propagator at three dwells, the white floor's density under a dwell change, a noise seed keeping the noise parameters, the sensor truth against the sensor current, the lattice through the wiring, the PSB occupations and relaxation — `tests/test_device_session.py`

## P1 tests (CI, `literature`: reproduce an item of a `confirmed` reference entry; ADR 0007)
- none

## Gallery protocols (rendered by CI, signed off per protocol by a human; ADR 0008)
- `session-charge-trace` — the charge a sample sees at a dot-lead degeneracy at three dwells, and the noisy read at 1 ms
- `session-hold-relaxation` — the mean charge after a hold relaxing on exp(−Γ_L t), and one trajectory
- each carries a `budget` (ADR 0009); both unsigned

## Exit criteria
- qdt CI green against the pinned phantom with `qdt-sim` under 300 lines — qdt side (CowboyQuantum/qdt#140, #141, #142)
- qdt's tune-up on the simulator reproduces its committed scenario outcomes after regeneration

## Open decisions
- cadence of phantom releases vs. qdt pins

## Status
IN PROGRESS — phantom side built (ADR 0024): the consumer surface, the contract test, the snapshot, the changelog and 0.2.0; two protocols unsigned. Open: the tag `v0.2.0` on `f72daac` (the owner's: this session's environment refused the tag push), the release cadence, the qdt side.
