# Architecture Decision Records

Immutable once merged (except a `Superseded by` status line). New decision = new ADR via `/adr <title>`.

| # | Title | Status |
|---|---|---|
| 0001 | Two tiers: compile once per device instance, run cheap at every gate point | Accepted |
| 0002 | Statistical realism, never prediction of a particular device | Accepted |
| 0003 | Dependency-free library with a reproducibility contract; qdt's adapter lives in qdt | Accepted |
| 0004 | Nine layers that import only downward; `compile` imports all | Accepted |
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
| 0017 | Readout dynamics: WKB lead rates from the landscape, a master equation with detailed balance, latching as lag | Accepted |
