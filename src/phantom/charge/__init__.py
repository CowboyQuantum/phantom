"""phantom.charge — Dots, barriers, few-electron Hubbard/CI; parameters as functions of gate voltage from the layers below. Runtime: small diagonalisations per gate point. Seeded from qdt-sim's DoubleDotModel (Step 1).

Modules: `ci` (constant-interaction double dot and its exact derived quantities), `dot` (`Dot`,
`Barrier`), `hubbard` (two-site extended Hubbard model whose t → 0 limit is the CI model).
"""

from phantom.charge.ci import (
    E_CHARGE_C,
    E_UNIT_V_PER_E_AF,
    K_B_EV_PER_K,
    ConstantInteractionParams,
)
from phantom.charge.dot import Barrier, Dot
from phantom.charge.hubbard import HubbardParams

MODEL_VERSION = 1
"""Layer model version: bumped whenever any output of the layer changes for the same inputs (ADR 0003)."""

__all__ = [
    "E_CHARGE_C",
    "E_UNIT_V_PER_E_AF",
    "K_B_EV_PER_K",
    "MODEL_VERSION",
    "Barrier",
    "ConstantInteractionParams",
    "Dot",
    "HubbardParams",
]
