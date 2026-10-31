# Gallery (ADR 0008, 0009)

A stratified set of *measurement protocols*: sweeps a device of this class is put through, described
as physics (gates, ranges, resolution, readout channel, time base, seed). CI renders every protocol on
every PR; nothing rendered is stored in git. A human signs each protocol off; a step is settled only when
all protocols its layers touch are signed.

**View it.** On a pull request: the preview link in the sticky comment on the PR, which points at
`https://cowboyquantum.github.io/phantom/pr-<n>/gallery/rendered/index.html` (ADR 0012); it updates a minute
or two after each push and is removed when the PR closes. On `main`: the docs site, "Gallery" in the nav, then
[the rendered gallery](rendered/index.html). The CI artifact `gallery-ubuntu-latest` holds the same files. Locally:

```
uv run python scripts/gallery.py          # renders to build/gallery/ (gitignored)
uv run mkdocs serve                       # browse it
```

The candidate roster for every step is `PROTOCOLS.md`; rows become protocol files when their step opens.

## Protocol file — `protocols/<id>.yaml`

| Field | Meaning |
|---|---|
| `id` | kebab-case, stable; the file name, e.g. `stability-coarse-2dot` |
| `name`, `description` | what the sweep is and what a device of this class should show |
| `layers` | phantom layers exercised; the protocol is owed by the step that delivers them |
| `instance` | which model and material class, and the `seed` that draws the instance (disorder is a draw) |
| `sweep` | ordered axes: `{gate_or_time, from, to, points, unit}`; one axis is a trace, two a map |
| `readout` | channel: `current`, `sensor_conductance`, `sensor_lockin` (demodulated response to a square-wave modulation of a named gate; amplitude and frequency are parameters), `rf_response`, `charge_state`, `spin_outcome`, with its unit |
| `render` | plot kind and axes labels; cosmetic only, never compared |
| `budget` | required, ADR 0009: `compile: {wall_s, solver_calls}`, `run: {wall_ms_per_point, points, model_evals_per_point, solver_calls: 0}` |
| `must_show` | plain-language checklist a reviewer looks for; filled in during review |

## Sign-off record — `signoff/<id>.yaml`

| Field | Meaning |
|---|---|
| `protocol` | the protocol id |
| `array_sha256` | canonical hash of the swept arrays: each array's `.tobytes()` in sorted key order, dtype and shape included |
| `generator_commit` | git SHA of the tree that produced the arrays |
| `model_versions` | `{layer: model_version}` for every layer involved |
| `seed` | the instance seed |
| `signed_by`, `signed_on` | GitHub handle and date |
| `notes` | what the review concluded the sweep must show, and what was changed to get there |

The runtime protocols run on a device compiled once and stored under `devices/` (ADR 0033), so every
platform runs the same compiled numbers. After a change to compile, `uv run python scripts/gallery.py
--record-devices` stores them again; a stale or missing device stops the gallery, and
`tests/test_gallery_devices.py` checks each stored device against a fresh compile within 1e-9.

CI computes `array_sha256` on `ubuntu-latest`, the reference platform (ADR 0010), and fails if it differs from a signed record; the Windows job compares its arrays to the reference within 1e-9 relative (ADR 0014). A protocol
with no record renders and is labelled *unsigned*. A session writes a record only when the human says
so, naming them; their approval of the PR is the signature.

## Budget ledger — `build/gallery/budget.json`

Per protocol and platform: declared and measured `compile` and `run` figures. Counts must match
exactly; wall time fails above 3× budget in CI (`--strict` locally for factor 1). The table is rendered
next to the pictures.
