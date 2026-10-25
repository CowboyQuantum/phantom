# Defect and material catalog

The knowledge only an experienced tuner has: which defects and material traits exist, what they perturb,
and what they look like during tune-up. Data, not code (ADR 0006). Schema: `SCHEMA.md`. Entries:
`entries/*.yaml`. One example entry is provided as a `draft` with every value marked TODO; the rest are
filed as GitHub issues by `scripts/issues_catalog.py` for a human to write.

Suggested entries (the issue set): interface fixed charge · oxide traps and two-level fluctuators · gate
misalignment and width error · gate-to-gate overlap and leakage · open or shorted gate · Ge-concentration
fluctuation (valley) · interface steps and roughness · strain variation · residual dopants · dielectric
thickness variation · 2DEG density and mobility by material class · sensor-dot coupling asymmetry ·
charge drift (slow) · material classes SiGe vs MOS as summary traits · (later) Ge holes.
