# 0016 — Wiring faults as a gate matrix and a leak network; sensor asymmetry as coupling draws; trapping as a history state

**Status:** Accepted
**Date:** 2026-10-26
**Step:** 3

## Context
Three catalog entries describe faults with trivial physics but common occurrence: an open or shorted gate
(#4), gate-to-gate leakage (#3) and sensor–dot coupling asymmetry (#11). A fourth signature, the
hysteresis loop of a gate sweep, comes from charge trapping with gate history (#1 at its slowest).
Step 3 needs generators for all four that fit the existing chain without new electrostatics.

## Decision
- **Wiring faults** (`disorder.wiring`) act on the gate vector before the basis: the voltages the gates
  actually carry are M·V + c for programmed V. A shorted pair copies the kept gate's row onto the
  follower; an open gate zeroes its row and sets a stuck potential. Protocols sweep *programmed*
  voltages, so the signature appears where a tune-up would see it.
- **Leakage** is a resistance network between adjacent gate pairs; the current out of each gate is the
  sum over its leaks, equal and opposite on the two ends. It adds to the dielectric leakage of
  `readout.transport` in the leakage protocols.
- **Sensor asymmetry** scales each dot's per-electron sensor coupling by a log-normal factor.
- **Adjacent pairs** are defined geometrically from the layout (polygons within 30 nm in the same band),
  so the fault model follows the design.
- **Trapping** (`disorder.trapping`) is a trapped fraction with first-order kinetics per sample along the
  gate sequence: charging above a capture voltage, releasing below a release voltage, holding between.
  The effective gate voltage is V − shift × f. It is static and history-dependent, the slow partner of
  the telegraph and 1/f processes in `readout.noise`; Step 4 ties the rates of both to the catalog.
- All four are part of the scenario draw (`wiring`), seeded, serialisable, and clean at severity 0;
  probabilities and parameters are placeholders (`calibration: none`) until the entries are written.

## Consequences
- Easier: the fault signatures (faint dot, gates moving together, floating barrier, leak current pair,
  hysteresis loop) render from the same layout and models as the clean device.
- Harder: a shorted or open gate can leave the intended dots unformed at the clean operating point;
  protocols that derive dots must handle the failure (they re-derive at the actual voltages).
- Must now: `readout.transport` turn-on protocols accept a trapping draw; Step 4 connects the trap
  kinetics to the fluctuator rates.

## Alternatives considered
- Faults as perturbations of the derived CI parameters — rejected: the signature of a short is in how
  *programmed* voltages map to the device, which the gate matrix captures and a CI tweak would not.
- A full trap ensemble with capture and emission time constants — deferred to Step 4.
