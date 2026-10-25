# Reference entry schema (ADR 0007)

One YAML file per paper, `docs/references/entries/<id>.yaml`.

| Field | Meaning |
|---|---|
| `id` | `<firstauthor><year>-<short>` in kebab-case, e.g. `hanson2007-spins-rmp`; the file name |
| `title`, `authors`, `year`, `venue` | as published |
| `doi` | the DOI (`10.xxxx/...`); `arxiv` allowed in addition, never instead, unless the paper has no DOI |
| `layers` | phantom layers the paper informs: `layout`, `electrostatics`, `disorder`, `valley`, `charge`, `readout`, `spin` |
| `provides` | list of items, see below |
| `status` | `proposed` (written by a session) or `confirmed` (a human opened the paper and checked every item) |
| `confirmed_by`, `confirmed_on` | GitHub handle and date; present only when `confirmed` |
| `notes` | what the paper does not say, approximations in the comparison, anything a reader should know |

Each item in `provides`:

| Field | Meaning |
|---|---|
| `label` | stable short label a test names, e.g. `fig3b`, `table1-row2`, `eq12` |
| `kind` | `figure` (a digitized curve), `table`, `value`, `range`, `formula` |
| `quantity` | what it is, with unit in the name or the text, e.g. `valley_splitting_ueV` |
| `conditions` | material class, temperature, field, dot size: whatever the paper states |
| `value`, `tolerance` | for `value`/`range`; the paper's own uncertainty, or the stated reading precision |
| `data` | for `figure`/`table`: path under `data/<id>/`, a CSV with a header row naming columns and units |
| `digitized` | `{method, by, on}` for curves read off a figure; `extracted` for data the paper publishes |

Rules: no internal data; a `confirmed` entry is never edited except to add items, which revert it to
`proposed` until confirmed again; a test may not hard-code a value that lives in the library.
