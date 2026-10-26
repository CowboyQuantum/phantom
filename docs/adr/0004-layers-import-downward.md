# 0004 — Nine layers that import only downward; `compile` imports all, nothing imports `compile`

**Status:** Accepted
**Date:** 2026-10-24
**Step:** 0

## Context
End-to-end realism needs layout, electrostatics, disorder, valley, charge, readout and spin, each with a
different cost and a different source of truth. Without a dependency order the layers would entangle and
the two-tier split (ADR 0001) would erode.

## Decision
`layout` ← `electrostatics` ← `disorder` ← `valley` ← `charge` ← `readout` ← `spin` ← `device` ← `optimum` (`device`, the compiled-device record, added by ADR 0023), each
importing only layers below it; `compile` may import everything; nothing imports `compile`.
`scripts/check_imports.py` enforces it, and the no-qdt rule, in CI. Each layer's `__init__` docstring states
what it owns and its runtime form.

## Consequences
- Easier: a consumer can use `charge` + `readout` without spin; layers are testable alone; the seam to
  extract or replace a layer (FEM electrostatics) is explicit.
- Harder: cross-layer couplings (valley affecting readout, charge noise affecting spin) flow through
  parameters passed downward-to-upward, never by upward import.
- Must now: keep the order in `CLAUDE.md`; a new layer is an ADR.

## Alternatives considered
- One flat `model` module — rejected: the entanglement this ADR exists to prevent.
