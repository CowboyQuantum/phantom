# Architecture Decision Records

Immutable once merged (except a `Superseded by` status line). New decision = new ADR via `/adr <title>`.

| # | Title | Status |
|---|---|---|
| 0001 | Two tiers: compile once per device instance, run cheap at every gate point | Accepted |
| 0002 | Statistical realism, never prediction of a particular device | Accepted |
| 0003 | Dependency-free library with a reproducibility contract; qdt's adapter lives in qdt | Accepted |
| 0004 | Nine layers that import only downward; `compile` imports all | Accepted (amended by 0023, 0024) |
| 0005 | Test ladder P0 / P1 / P2 | Superseded by 0007 |
| 0006 | The defect and material catalog is data a human fills in | Proposed |
| 0007 | Published results only; every claim names an entry in the reference library | Accepted |
| 0008 | Gallery of measurement protocols, rendered by CI, signed off per protocol by a human | Accepted |
| 0009 | Every gallery protocol carries a time and compute budget, measured in CI | Accepted |
| 0010 | Bit-for-bit reproducibility per platform; across platforms within 1e-12; ubuntu-latest is the reference | Accepted |
| 0011 | A running map of how the models fit together lives in the docs and is kept current | Accepted |
| 0012 | Every pull request publishes a preview of the docs and the rendered gallery; main publishes the root | Accepted |
| 0013 | Closed-form electrostatics: solid-angle potentials, effective depth, dots as harmonic minima | Accepted |
| 0014 | The cross-platform tolerance is 1e-9 relative, because derived parameter chains amplify last-bit differences | Accepted |
| 0015 | Disorder generators: interface charges as image-charge offsets, lithography faults as polygon transforms, one seeded scenario draw | Accepted |
| 0016 | Wiring faults as a gate matrix and a leak network; sensor asymmetry as coupling draws; trapping as a history state | Accepted |
| 0017 | Readout dynamics: WKB lead rates from the landscape, a master equation with detailed balance, latching as lag | Accepted |
| 0018 | RF reflectometry as a lumped tank loaded by the sensor; single-shot readout as integrated noisy shots | Accepted |
| 0019 | Noise is drawn per device: a 1/f amplitude, a Poisson number of dominant switchers, a white floor | Accepted |
| 0020 | Valley splitting is a Rice draw per dot from alloy disorder, with a deterministic part for MOS; three signatures | Accepted |
| 0021 | Exchange from the Hubbard two-electron sector; Lindblad spins with quasi-static noise by quadrature; PSB window is the valley draw | Accepted |
| 0022 | Blockade transport as sequential cycles; funnel as two-level mixing; Elzerman shots; noise spectroscopy by filter functions; gate fidelity as the Lindblad channel's entanglement fidelity | Accepted |
| 0023 | A compiled device is a parameter record in its own layer; the basis is recomputed, not stored; JSON serialisation; the FEM backend is a hook; optimum reads the record | Accepted |
| 0024 | The consumer contract: `phantom.api` above `compile`, a session with a forward clock, knowable fields apart from truth, a contract test consumers run, a checked surface snapshot, 0.x releases | Accepted (amended by 0025) |
| 0025 | A session tags a release when the owner says so, on a green `main` with the changelog entry; the owner can always tag | Accepted |
| 0026 | The sensor is a dot with its own barriers and bath: compiled leads, a lifetime-broadened current, its charging energy and coupling from the layout | Accepted |
| 0027 | The qubit dots' rates follow the gates at runtime: compiled slopes of the WKB action at a fixed tunnelling level, capped at the attempt rate | Accepted |
| 0028 | The sensor dot's charge is a process in time: a birth–death chain over its levels on the session clock, counted electrons per dwell, back-action on the dots against its mean | Accepted |
| 0029 | The interdot coupling follows the gates and hybridizes the charge states: t from half the slopes of ln Γ_m, capped at ħω₀/π; eigenstates of the CI window in the statics and a Pauli master equation over them in the dynamics | Accepted |
| 0030 | The barrier gates reshape the dots through a compiled grid: the double dot re-derived along each barrier gate, its changes interpolated and added at runtime, held where the dots merge | Accepted |
| 0031 | The exchange of the dots' valence electrons follows the gates: the singlet secular equation on the CI window's local sector at the gated coupling and the shaped dots; the Pauli readout at the gated coupling | Accepted |
