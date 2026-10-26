# 0015 — Disorder generators: interface charges as image-charge energy offsets, lithography faults as polygon transforms, one seeded scenario draw

**Status:** Accepted
**Date:** 2026-10-26
**Step:** 3

## Context
Step 3 turns catalog entries into generators that perturb the device for physical reasons. The catalog
entries themselves (#1–#14) are the owner's to write, and only the draft `interface-fixed-charge` exists,
with placeholder densities. The generators must be buildable now, on parameters, and must fit the
closed-form electrostatics of Step 2 (ADR 0013) without a solver.

## Decision
- **Interface charges** (`disorder.charges`): a Poisson number of point charges ±e at the interface depth
  over the layout's bounding box, density × severity. Each adds to the electron energy in the electron
  plane the potential of the charge plus its image in the pinned surface (the same half-space model as
  the gate potentials), with a smoothing length for the finite thickness of the electron layer so that
  charges in the plane itself (MOS) stay finite. No screening by the electrons; the perturbation is an
  upper bound and the entry is `calibration: none` until a confirmed reference sets density and sign.
- **Gate misalignment and width error** (`disorder.geometry_faults`): a common shift of the gate layer
  and an independent outward edge displacement per gate, applied to the polygons. The perturbed
  lever-arm basis follows from the electrostatics with no new physics.
- **The landscape carries disorder as an offset**: `GateBasis.offset_ev` is a gate-independent energy
  added to the electron energy; everything above (dot finding, derived CI parameters, sensor coupling)
  works unchanged on a perturbed basis.
- **One scenario draw** (`disorder.scenario`): material-class parameters + seed + severity → all
  generators from one generator stream, in a fixed order; `to_dict()`/`from_dict()` make the draw
  plain data a consumer can store. Severity 0 draws nothing and returns the clean layout and basis,
  bit for bit (P0 test).
- **Constants** live in `phantom.constants`, below every layer, so `disorder` imports nothing from
  `charge` (rule 2).
- **Protocols**: the landscape and the coarse stability map, clean beside three draws; and the mean
  number of unintentional minima against charge density over many draws, which the Step 3 P1 row will
  compare with a published MOS/SiGe comparison.

## Consequences
- Easier: disorder is a few hundred microseconds of compile-tier work per draw; every Step 1 and
  Step 2 protocol can be rendered for a draw by passing a perturbed basis.
- Harder: without screening, charge densities at the MOS placeholder give large perturbations; the
  numbers are placeholders until the catalog is written, and nothing here claims a rate.
- Must now: strain and g-factor draws wait for the spin layer; switchers and drift remain noise
  processes in `readout` (Step 4 ties their rates to the catalog); the catalog entry gains
  `calibration: none` and `tier: runtime`.

## Alternatives considered
- Random potential fields with a chosen correlation length — rejected: not traceable to a catalog quantity.
- Screened (Thomas–Fermi) charge potentials — deferred to Step 7's solver tier (ADR 0001).
