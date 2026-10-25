# phantom — Development Plan & Roadmap

**Status:** v0.1, for review · **Owner:** Captain · **Date:** 2026-10-24
**Scope:** the physics of a gate-defined spin-qubit device as a library: from gate layout and material
stack to charge, readout and spin, with disorder and fabrication defects drawn from a catalog. Consumers:
qdt (as a simulated backend), agent evaluation, hosted training scenarios, device-physics studies.

## 1. Purpose
Real tune-up is hard because of what is not known about a device: its disorder, its defects, its valley
splittings. No model predicts a particular device's tune-up at runtime, and this one does not try. phantom
produces **device instances** that behave like real devices of their class — the right kinds of faults at
the right rates, the right sweep-rate dependence, the right readout behaviour — and runs them cheaply enough
to be swept ten thousand times by a tuning loop. Its claims are statistical (ADR 0002).

## 2. Two tiers (ADR 0001)
- **Compile once.** Expensive physics runs once per device instance: electrostatics for a layout and stack,
  screening, disorder placement, valley-splitting draws. Output: a *compiled device* — parameter fields and
  functions the runtime consumes (lever-arm fields, Hubbard parameters as functions of gate voltage, tunnel
  rates, sensor coupling, noise spectra).
- **Run cheap.** The runtime layers evaluate a compiled device at a gate point in well under a millisecond:
  charge ground state, rates, sensor signal, noise samples, spin dynamics for a few qubits.
A tune-up never waits on a solver; a scenario is compiled before it is served.

## 3. Layers (ADR 0004) — each imports only downward
| Layer | Owns | Runtime-capable form | Offline form |
|---|---|---|---|
| `layout` | Gate polygons, stack (2DEG depth, dielectrics), material class (SiGe, MOS) | — | — |
| `electrostatics` | V(r; gate vector) and its basis per gate | Closed-form polygonal-gate potentials over a 2DEG (Davies–Larkin–Stern) | FEM Poisson, Thomas–Fermi screening, via a wrapped external solver |
| `disorder` | Static potential perturbations and their generators, from the catalog | Precomputed perturbation fields | Placement draws |
| `valley` | Valley splitting per dot (Si) | A draw from a distribution for the material and dot size | — (atomistic is out of scope) |
| `charge` | Dots, barriers, few-electron Hubbard/CI, parameter functions t(V), U, V_ij | Diagonalize small matrices per gate point | — |
| `readout` | Lead tunnel rates, rate equations, latching, charge sensor (SET) response, RF reflectometry, noise processes (fluctuators, 1/f, telegraph, drift) | All cheap | — |
| `spin` | Zeeman with gate-dependent g, exchange J(ε, t_c), spin–orbit/EDSR, hyperfine, T₁, T₂*, Pauli spin blockade | Lindblad for ≤3 qubits; Monte Carlo over noise | — |
| `optimum` | For a compiled device: the best achievable operating point, lattice, sensor sensitivity, and (with spin) readout fidelity and coherence | Cheap | — |
| `compile` | Turns layout + material + a catalog draw + seed into a compiled device | — | Minutes per instance |

`qdt-sim`'s `DoubleDotModel` and noise models are the seed of `charge` and `readout` (Step 1 migrates them).

## 4. What a consumer gets (ADR 0003)
- `compile(layout, material, draw, seed) -> Device` — a compiled, serialisable device instance with a
  `model_version`.
- `Device.at(gate_vector) -> State` (charge configuration, rates, sensor signal); `Device.measure(...)`
  for readout and spin outcomes; noise advanced by explicit time steps.
- `Device.truth` — the hidden parameters, and `optimum(Device)` — the best achievable results. Consumers
  that must hide truth from a user (a tuning loop, a product) simply do not hand it over.
- Reproducibility: parameters + seed + version → identical output.
- No instruments, no measurement vocabulary, no users: the adapter that turns a `Device` into a qdt stand
  lives in qdt and pins a phantom release.

## 5. The defect and material catalog (ADR 0006)
`docs/catalog/` is data a human fills in: each entry names a defect or material trait, the layer it
perturbs, the parameter and its distribution, a severity knob, the signature it leaves in a tune-up, and a
source. `disorder` and `compile` read it. The catalog is the part only an experienced tuner can write; its
entries are filed as issues for that reason.

## 6. Validation: test ladder, reference library, gallery (ADR 0007, 0008, 0009; 0005 superseded)
| Level | What | Where |
|---|---|---|
| P0 | Analytic limits and invariants: a single dot's Coulomb diamonds, two-dot honeycomb geometry, lever arms from a known capacitance matrix, detailed balance in rate equations, Lindblad trace preservation | CI |
| P1 | Reproductions of items in `docs/references/` (confirmed by a human): valley-splitting distributions, charge-noise spectra, PSB fidelity vs. valley splitting, exchange vs. detuning | CI, marked `literature(id, item)` |
| Gallery | Measurement protocols (turn-on, pinch-off, diamonds, stability maps, sensor peaks, noise traces and spectra, later spin readout) rendered with a fixed seed; a human signs each off on its data hash; each carries a time and compute budget | CI renders and publishes; `docs/gallery/`; sign-off per protocol gates step exit |

No comparison with internal device statistics for now (ADR 0007). Quantities without a published standard are marked `calibration: none`.

## 7. Roadmap — operational detail in `docs/steps/`
| Step | Name | Exit in one line | Detail |
|---|---|---|---|
| 0 | Bootstrap | Repo, CI, layer skeletons, catalog template, first ADRs | `steps/00-bootstrap.md` |
| 1 | Charge and readout core | `qdt-sim`'s model and noise migrated; parameter functions; reproducibility contract proven | `steps/01-charge-readout-core.md` |
| 2 | Layout and closed-form electrostatics | A layout + stack gives V(r) and lever-arm fields in microseconds; two-dot lattice emerges from a layout | `steps/02-layout-electrostatics.md` |
| 3 | Disorder and the fault models | Catalog draws perturb the potential; unintentional dots, leaky barriers, switchers, drift appear for physical reasons | `steps/03-disorder-faults.md` |
| 4 | Readout physics | Lead rates, latching, sweep-rate dependence, sensor and RF response, 1/f and telegraph noise from fluctuator ensembles | `steps/04-readout.md` |
| 5 | Valley | Statistical valley splitting per dot; PSB consequences | `steps/05-valley.md` |
| 6 | Spin | Few-qubit dynamics under pulses; readout fidelity, T₁, T₂* | `steps/06-spin.md` |
| 7 | Compile and optimum | Offline FEM electrostatics wrapped; `compile()`; `optimum()` for every layer | `steps/07-compile-optimum.md` |
| 8 | qdt cutover | `qdt-sim` becomes an adapter pinning phantom; qdt's golden files regenerated once | `steps/08-qdt-cutover.md` |

Steps 1–4 give faithful faults and are the body of the work; 5–6 extend to spin; 7 closes the loop; 8 is
coordinated with qdt. 1 and 2 can start together.

## 8. Open decisions
| # | Decision | State |
|---|---|---|
| 1 | Which external electrostatics solver to wrap at Step 7 | Open; closed-form first |
| 2 | Hubbard vs. full CI as the default charge solver | Open; Hubbard first, CI where PSB needs orbitals |
| 3 | Serialisation of a compiled device (format, size) | Open; must be versioned and diffable |
| 4 | Catalog file format and schema (YAML proposed) | Open; ADR 0006 proposes |
| 5 | Whether holes (Ge) join the material classes | Open; electrons in Si first |

## 9. Success criteria
1. A tune-up loop runs against a phantom at ≥ 10⁴ gate points per second per core.
2. The fault catalog's entries each produce their signature in a stability diagram, with a P0 or P1 test.
3. Reproducibility holds across platforms in CI.
4. qdt's adapter is under 300 lines and contains no physics.
5. A P1 suite of at least five published results reproduces within stated tolerance.
