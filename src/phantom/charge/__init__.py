"""phantom.charge — Dots, barriers, few-electron Hubbard/CI; parameters as functions of gate voltage from the layers below. Runtime: small diagonalisations per gate point. Seeded from qdt-sim's DoubleDotModel (Step 1).

Modules: `ci` (constant-interaction double dot and its exact derived quantities), `dot` (`Dot`,
`Barrier`), `hubbard` (two-site extended Hubbard model whose t → 0 limit is the CI model), `hybrid` (the
interdot coupling at the gates and the CI window's hybridized eigenstates, #58), `derived` (CI parameters
from a layout through `electrostatics`, Step 2).
"""

from phantom.charge.ci import ConstantInteractionParams
from phantom.charge.derived import DerivedDoubleDot, derive_double_dot, mutual_capacitance_af
from phantom.charge.dot import Barrier, Dot
from phantom.charge.hubbard import HubbardParams
from phantom.charge.hybrid import InterdotCoupling
from phantom.constants import E_CHARGE_C, E_UNIT_V_PER_E_AF, K_B_EV_PER_K

MODEL_VERSION = 4
"""Layer model version: bumped whenever any output of the layer changes for the same inputs (ADR 0003)."""

__all__ = [
    "E_CHARGE_C",
    "E_UNIT_V_PER_E_AF",
    "K_B_EV_PER_K",
    "MODEL_VERSION",
    "Barrier",
    "ConstantInteractionParams",
    "DerivedDoubleDot",
    "Dot",
    "HubbardParams",
    "InterdotCoupling",
    "derive_double_dot",
    "mutual_capacitance_af",
]
