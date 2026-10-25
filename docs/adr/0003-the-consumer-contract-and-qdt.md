# 0003 — phantom is a dependency-free library with a reproducibility contract; qdt consumes it through an adapter that lives in qdt

**Status:** Accepted (rule 3 amended by 0010)
**Date:** 2026-10-24
**Step:** 1, 8

## Context
qdt's `qdt-sim` holds a double-dot model and noise models that implement qdt's instrument protocols. qdt
decided (its ADR 0078) that the physics core stays a dependency-free module and is extracted once a second
consumer exists; agent evaluation and hosted scenarios are those consumers, and both need ground truth and
an optimum that must never travel through qdt's API.

## Decision
- phantom imports nothing from qdt, measurement software or instrument SDKs; it is numpy/scipy in, numpy
  out. Its public API is plan §4: `compile(...) -> Device`, `Device.at`, `Device.measure`, `Device.truth`,
  `optimum(Device)`, all versioned.
- **Reproducibility contract:** parameters + seed + `model_version` → bit-identical output on every
  platform CI runs. A change that alters output bumps `model_version`.
- **qdt's adapter lives in qdt** (`qdt-sim`): it maps `StandConfig`/`ChannelMap` onto a `Device`, implements
  qdt's capability protocols over it, pins a phantom release, and contains no physics (target: < 300 lines).
  phantom ships a contract test the adapter also runs.
- **Truth stays with the consumer.** phantom exposes `truth` and `optimum`; a consumer that must hide them
  from a user (a tuning loop, a product, an agent) does not forward them. phantom never has a "hidden" mode.
- Changes that break the adapter are breaking releases; qdt pins and upgrades deliberately (Step 8).

## Consequences
- Easier: three consumers share one physics with one version to cite; qdt's golden files depend on a
  pinned release, not a moving module.
- Harder: physics changes reach qdt only through a release and a regeneration of its fixtures.
- Must now: Step 8 cutover in qdt; a changelog; the contract test.

## Alternatives considered
- Keep physics in qdt-sim and have consumers import qdt — rejected: pulls qdt's schemas and protocols into
  every consumer; ground truth would sit behind qdt's API.
- phantom exposes qdt's protocols directly — rejected: couples the library to one consumer's interface.
