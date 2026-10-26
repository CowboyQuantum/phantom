"""phantom.disorder — Static perturbations of the potential and parameters drawn from the catalog (`docs/catalog/`): interface charges, misalignment, strain. Runtime: precomputed fields; the draw itself happens in `compile`.

Modules: `charges` (random interface charges as image-charge energy offsets), `geometry_faults` (gate
misalignment and width error as polygon transforms), `scenario` (seeded draw of all generators,
serialisable; severity 0 reproduces the clean device bit for bit).
"""

from phantom.disorder.charges import (
    FixedChargeDraw,
    FixedChargeParams,
    charge_energy_ev,
    draw_fixed_charges,
)
from phantom.disorder.geometry_faults import (
    MisalignmentDraw,
    MisalignmentParams,
    apply_misalignment,
    draw_misalignment,
)
from phantom.disorder.scenario import ScenarioDraw, ScenarioParams, apply_scenario, draw_scenario

MODEL_VERSION = 1

__all__ = [
    "MODEL_VERSION",
    "FixedChargeDraw",
    "FixedChargeParams",
    "MisalignmentDraw",
    "MisalignmentParams",
    "ScenarioDraw",
    "ScenarioParams",
    "apply_misalignment",
    "apply_scenario",
    "charge_energy_ev",
    "draw_fixed_charges",
    "draw_misalignment",
    "draw_scenario",
]
