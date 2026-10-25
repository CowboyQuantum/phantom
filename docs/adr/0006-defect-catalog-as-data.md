# 0006 — The defect and material catalog is data a human fills in; generators read it

**Status:** Proposed
**Date:** 2026-10-24
**Step:** 0, 3

## Context
The map from fabrication defects and material traits to model perturbations is the knowledge an experienced
tuner has and a model does not. Written as code it would be scattered across generators and hard to review;
written as data it can be reviewed, cited, extended and argued about entry by entry.

## Decision
`docs/catalog/` holds YAML entries (schema in `docs/catalog/SCHEMA.md`): id, name, material classes it
applies to, the layer it perturbs, the parameter(s), the distribution and its defaults, a severity knob,
the signature it leaves in a tune-up, a source, and a status (`draft|reviewed|runtime|offline|rejected`).
`disorder` and `compile` read entries whose status is `runtime` or `offline`. Entries are filed as issues
for a human to write (`scripts/issues_catalog.py`); a session may draft an entry but a human marks it
`reviewed`. A catalog entry never contains a real device's values (CLAUDE.md rule 6).

## Consequences
- Easier: review per entry; the fault library grows without code changes once generators exist.
- Harder: a schema to maintain; generators must handle unknown fields gracefully.
- Must now: the template and schema in Step 0; the issue set filed.

## Alternatives considered
- Faults as code only — rejected: unreviewable by the person who knows them.
