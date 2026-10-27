# 0023 — A compiled device is a parameter record in its own layer; the basis is recomputed, not stored; JSON serialisation; the FEM backend is a hook; `optimum` reads the record

**Status:** Accepted
**Date:** 2026-10-26
**Step:** 7

## Context
ADR 0001 promised a `compile()` that runs every expensive draw and derivation once and a `Device` the
runtime layers consume without a solver. Plan §8 left two decisions open: which external solver to
wrap (#1) and how to serialise a device (#3). Rule 2 says nothing imports `compile`, yet `optimum`
must read a compiled device.

## Decision
- **A `device` layer** (`phantom.device`) between `spin` and `optimum` holds the `Device` record and
  the builders of every runtime object (CI model, sensor, barrier, Hubbard model, rate equation,
  sensed double dot with the noise draw, spin draw). `compile` writes it; `optimum` and consumers
  read it; the layer order in rule 2 and `scripts/check_imports.py` gains `device`.
- **What `compile_device` does, in order**: scenario draw (charges, misalignment, wiring, noise) →
  basis on the clean grid with the draw applied → programmed point mapped through the wiring faults →
  dots, sensor and CI model from the basis with the sensor gate moved onto its peak's flank → barriers
  to the leads and between the dots at the tunnelling level, their WKB transmissions by the action
  integral along the cut (`electrostatics.barrier_action`, which follows the barrier's shape instead
  of the rectangular height × width estimate, so t(V_B) is smooth) → the interdot t ≈ ħω₀ √T / π and
  the barrier gate's exponential scale from two nearby evaluations → valley splittings from the dots'
  radii → g-factors. Every draw takes a sub-seed derived from the device seed.
- **Serialisation** (plan §8 #3): JSON, sorted keys, one entry per line, floats at full precision,
  with a `format_version` and the model version of every layer that produced the device; the basis
  grid is not stored, `Device.basis()` recomputes it deterministically from the layout and the draw.
  A device file is a few kilobytes and diffs line by line.
- **Backends** (plan §8 #1 stays open): `closed_form` is the only one; `fem` is a hook behind the
  `[fem]` extra that raises with the open decision's name until a solver is chosen. Everything
  downstream of the basis is backend-agnostic.
- **`optimum`** reads the record and answers per layer: sensor flank by a coarse scan and a fine scan of |dI/dV|, the answer a grid
  point (a golden-section search and a parabolic interpolation both differed between platforms); the virtual-gate matrix from the inverse of the plunger block of the
  lever-arm matrix scaled to the raw plungers' lever arms, with cross-talk compensation for the other
  gates and the sensor gate holding the sensor potential fixed; the barrier voltage for a target
  t from the compiled exponential; the Pauli-blockade read point from the receiving dot's valley draw
  and the Elzerman peak time from the compiled lead rate; the exchange sweet spot in the operating
  cell (the symmetric point shifted by whole electrons to the cell the device sits in), its T₂* and
  SWAP fidelity with the device's noise. Brute force attains each (P0).
- **Screening**: still absent from the closed-form basis; it is the FEM backend's job (plan §8 #1),
  not a runtime approximation.

## Consequences
- Easier: one call compiles a device; one file carries it; every runtime object comes from it; the
  layer rule stays checkable.
- Harder: a tenth layer. The rectangular WKB estimate of ADR 0017 remains only as the bound in
  `readout.rates` and its P0 limits; the Step 4 dynamics protocols moved to the action integral
  (`readout.wkb_rate_from_action_hz`) in the second Step 7 PR, with two operating points shifted so
  their pictures keep the same rate-times-dwell regimes.
- Must now: P0 round trip (text-identical after reload; identical stability map from the reloaded
  device), the clean device bit-identical to the Step 2 derivation, `optimum` attained by brute force;
  protocols `compiled-device-roundtrip`, `virtual-gates-map`, `optimum-sensor-flank`.

## Alternatives considered
- Storing the basis grid in the device (hundreds of megabytes per device) — rejected; recomputation
  is deterministic and takes well under a second for the closed form.
- Pickle or NPZ — rejected: not diffable, not versioned by construction.
- Letting `optimum` import `compile` — rejected: rule 2's "nothing imports compile" is what keeps
  the offline tier out of every consumer's runtime.
