# Catalog schema (ADR 0006)

One YAML file per entry under `docs/catalog/entries/`, named `<id>.yaml`. Fields:

| Field | Meaning |
|---|---|
| `id` | kebab-case, stable; the file name |
| `name` | short human name |
| `kind` | `defect` (something went wrong in fabrication) or `trait` (a property of the material class or design) |
| `materials` | list of material classes it applies to: `SiGe`, `MOS`, later `Ge` |
| `layer` | the phantom layer it perturbs: `layout`, `electrostatics`, `disorder`, `valley`, `charge`, `readout`, `spin` |
| `parameters` | list of `{name, unit, distribution, defaults}`; distributions name their family and parameters |
| `severity` | how the knob scales the distribution (e.g. "multiplies density", "widens spread") |
| `signature` | what a tuner sees: free text, plus `in: [stability, diamonds, sensor, psb, coherence]` |
| `rate` | how often a device of this class shows it, as an order of magnitude with a source, or `unknown` |
| `tier` | `runtime` (a generator in `disorder`/`readout`) or `offline` (needs `compile`) |
| `sources` | citations (DOI or arXiv); textbooks allowed; never internal data |
| `status` | `draft` → `reviewed` (human) → `runtime` / `offline` (generator exists) / `rejected` |
| `notes` | anything else, including what the model deliberately ignores |

Rules: no real device, wafer or site values; published values cite their paper; a session may write a
`draft`, only a human sets `reviewed`.
