"""phantom.electrostatics — V(r; gate vector) from a layout: a per-gate basis by closed-form polygonal-gate potentials over a 2DEG (runtime, microseconds); FEM with screening arrives through `compile` (ADR 0001). Also the 2D effective-mass solve for dot positions and level spacings.

Modules: `closed_form` (solid-angle potential of polygon gates, Davies–Larkin–Stern; `GateBasis`),
`dots` (minima, harmonic radii, lever arms, barrier heights from the landscape), `schrodinger`
(offline finite-difference levels).
"""

from phantom.electrostatics.closed_form import (
    GateBasis,
    polygon_solid_angle,
    rectangle_potential_dls,
)
from phantom.electrostatics.dots import (
    DerivedBarrier,
    DerivedDot,
    barrier_action,
    barrier_between,
    barrier_width_nm,
    find_dots,
    line_cut,
)

MODEL_VERSION = 1

__all__ = [
    "MODEL_VERSION",
    "DerivedBarrier",
    "DerivedDot",
    "GateBasis",
    "barrier_action",
    "barrier_between",
    "barrier_width_nm",
    "find_dots",
    "line_cut",
    "polygon_solid_angle",
    "rectangle_potential_dls",
]
