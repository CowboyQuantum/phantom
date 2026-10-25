# 0002 — phantom claims statistical realism, never prediction of a particular device

**Status:** Accepted
**Date:** 2026-10-24
**Step:** all

## Context
What makes a device hard to tune is what is unknown about it: its disorder, defects and valley splittings.
A model cannot know them, and the tuning loop exists to discover them. Overclaiming here would discredit the
library with the people it serves.

## Decision
A phantom instance behaves like a device of its class — the right kinds of faults at the right rates, the
right sweep-rate dependence, the right readout behaviour — because its disorder is a *draw* from catalog
distributions. phantom predicts classes, not devices. Documentation, APIs and any grading built on
`optimum` state this; nothing in the library is named or described as predicting a specific device.

## Consequences
- Easier: the fidelity target is clear (P1 reproductions and P2 statistics), not a device-matching task.
- Harder: consumers wanting a digital twin of one device must fit phantom's parameters to it themselves.
- Must now: the P2 level compares distributions, not single devices; `truth` is labelled as the draw.

## Alternatives considered
- A digital twin fitted per device — out of scope; a consumer may build one on top.
