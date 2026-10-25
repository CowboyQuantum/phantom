# 0007 — Validate against published results only; every claim names a reference entry

**Status:** Accepted (supersedes 0005)
**Date:** 2026-10-25
**Step:** all

## Context
ADR 0005 defined three rungs: P0 analytic invariants, P1 published reproductions, P2 comparisons with
internal device statistics. P2 was never an exit criterion and never entered the repo, but its existence
let a quantity look "validated" by data nobody outside could see. For now the owner does not want the
simulator matched against internal device data at all. The remaining standard is the literature, and the
comparison with it must be explicit enough that a reader can open the paper.

## Decision
- **Two rungs.** P0 (analytic limits and invariants) and P1 (reproductions of published results) are the
  whole ladder. P2 is removed: the `internal` marker, the `-m internal` command and every mention of
  device statistics as a check go. A later ADR may bring it back; nothing here anticipates it.
- **A reference library**, `docs/references/`, is the only admissible standard. One entry per paper
  (`entries/<id>.yaml`, schema in `docs/references/SCHEMA.md`): DOI, the layers it informs, and the
  concrete items it provides (a figure, a table, a value, a range) with the paper's own tolerance.
  Digitized curves live beside the entry under `data/<id>/` with a note on how they were digitized.
- **Every claim cites an id.** A P1 test is `@pytest.mark.literature("<id>", item="<label>")` and loads
  its reference curve from the library, never from a literal in the test. A catalog entry's `sources` and
  `rate` fields point to ids. Benchmarks and docs that claim agreement with physics name the id.
  `scripts/check_references.py` (in `check_all.sh`) fails on an id that does not exist, an entry whose
  id or DOI is malformed, or a `literature` test without an id.
- **Agent-proposed, human-confirmed.** An entry is `proposed` when a session writes it and `confirmed`
  only after a human has opened the paper and checked the DOI, the figure and the digitized values.
  A P1 test against a `proposed` entry is skipped, not passed; a step cannot exit on it.
- **Uncalibrated is said out loud.** A quantity with no published standard (some fault rates will be)
  is marked `calibration: none` in its catalog entry and its docstring. Nothing downstream may describe
  it as realistic. This is the cost of dropping P2 and it is stated rather than hidden.

## Consequences
- Easier: the repo's fidelity claim is auditable from public sources alone; provenance (rule 6) is
  trivially satisfied; a reader can reproduce every number.
- Harder: fault rates and other device-population quantities are uncalibrated until a paper supplies
  them; a model can be fully "validated" and still be unconvincing, which ADR 0008 addresses.
- Must now: `check_references.py` exists before the first P1 test; a session that writes an entry
  cannot flip it to `confirmed`; the catalog schema's `sources` field takes ids, not free text.

## Alternatives considered
- Keep P2 as documented in 0005 — rejected by the owner for now.
- Free-text citations in tests — rejected: not checkable, and the digitized data would scatter.
