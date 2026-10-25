# 0008 — A gallery of measurement protocols, rendered by CI, signed off per protocol by a human

**Status:** Accepted (publishing amended by 0012: previews per pull request)
**Date:** 2026-10-25
**Step:** all, from Step 1

## Context
qdt has a front end, so a human could look at what the simulator produced. phantom has none and will
not grow one (rule 7, ADR 0003). A passing P0 or P1 test says a number is right within tolerance; it
does not let a human notice that a stability diagram looks wrong in a way nobody thought to assert.
That judgement is the owner's, exercised in conversation with the agent, and it needs a durable record
without an interactive UI in this repo.

## Decision
- **Protocols, not tasks.** The gallery is stratified by *measurement protocol*: a sweep definition in
  physics terms (gates, ranges, resolution, readout channel, time base, seed), one YAML file per
  protocol under `docs/gallery/protocols/`. Turn-on curves, barrier pinch-off, coarse and fine charge
  stability maps, Coulomb diamonds, sensor peaks, detuning cuts, noise time traces and spectra, later
  the spin-readout sweeps. These are the scans a tune-up uses, described as physics; a consumer maps
  them onto its own stages. Nothing here is named for a tuning step.
- **CI renders; git stores no images.** `scripts/gallery.py` runs every protocol with its fixed seed and
  writes, outside git, the swept arrays (`.npz`), a PNG per protocol, the budget ledger (ADR 0009) and
  a mkdocs page. A `gallery` CI job uploads this as a workflow artifact on every PR and, on `main`,
  deploys the docs site with the gallery to GitHub Pages. Locally `uv run python scripts/gallery.py`
  then `uv run mkdocs serve`. matplotlib lives in the `gallery` dependency group, never in runtime deps.
- **Sign-off is on the data.** A human signs a protocol by a record in `docs/gallery/signoff/<id>.yaml`:
  the canonical hash of the swept arrays, the generator commit, the `model_version`s involved, the seed,
  the date, the signer's GitHub handle, and the notes that came out of the review (what the sweep must
  show). CI recomputes the hash and fails if it differs from a signed record. Changing the physics
  therefore requires a new sign-off, which is the point.
- **Who writes the record.** The human reviews the rendered gallery and discusses it with the agent in
  chat; the agent adjusts the protocol or the model and re-renders until the human says it is right.
  At the human's instruction the agent writes the sign-off record naming them; the human's approval of
  that PR is the signature. A session never signs on its own initiative.
- **Step exit.** A step is settled only when every protocol its layers touch is rendered and signed off,
  in addition to P0 and P1 green (ADR 0007). The step files carry a fixed heading `Gallery protocols`.
- **Determinism.** The arrays are covered by the reproducibility contract (rule 3) and must hash
  identically on both CI platforms. The PNG is for eyes only and is never compared.

## Consequences
- Easier: a human gate on realism without a UI in this repo; review happens where code review happens;
  the record of *why a sweep looks the way it does* survives the chat that produced it.
- Harder: one more CI job and a Pages deployment; every physics change that moves a signed array needs
  a human in the loop before merge.
- Must now: `scripts/gallery.py`, the protocol and sign-off schemas (`docs/gallery/README.md`), the CI
  job and the mkdocs nav exist in Step 1; each step lists its protocols before it starts.

## Alternatives considered
- Commit the PNGs — rejected by the owner: binaries in git; CI can render them on demand.
- Pixel comparison of PNGs — rejected: matplotlib output is not platform-stable; the arrays are.
- A notebook or interactive viewer in this repo — rejected: rule 7, and a maintenance burden phantom
  does not want.
