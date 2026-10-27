"""Electrostatics backends for `compile` (ADR 0023).

`closed_form` is the Davies–Larkin–Stern solid-angle basis of `electrostatics.closed_form` (ADR 0013)
evaluated on the clean layout's grid with the scenario draw applied (`disorder.apply_scenario`). `fem`
is the hook for a wrapped external Poisson solver behind the `[fem]` extra (plan §8 #1, open): no
solver is written here, and until one is chosen the backend raises with that message. Both return a
`GateBasis`, so everything downstream of the basis is backend-agnostic.
"""

from __future__ import annotations

from typing import Literal

from phantom.disorder.scenario import ScenarioDraw, apply_scenario
from phantom.electrostatics.closed_form import GateBasis
from phantom.layout.geometry import Layout

Backend = Literal["closed_form", "fem"]


def closed_form_basis(layout: Layout, draw: ScenarioDraw, grid_pitch_nm: float) -> tuple[Layout, GateBasis]:
    """The perturbed layout and its basis with the interface charges' energy offset."""
    return apply_scenario(layout, draw, grid_pitch_nm)


def fem_basis(layout: Layout, draw: ScenarioDraw, grid_pitch_nm: float) -> tuple[Layout, GateBasis]:
    raise NotImplementedError(
        "the FEM backend wraps an external Poisson solver behind the `phantom[fem]` extra; no solver has been "
        "chosen yet (DEVELOPMENT_PLAN §8 #1). Use backend='closed_form'."
    )


def basis_for(backend: Backend, layout: Layout, draw: ScenarioDraw, grid_pitch_nm: float) -> tuple[Layout, GateBasis]:
    if backend == "closed_form":
        return closed_form_basis(layout, draw, grid_pitch_nm)
    if backend == "fem":
        return fem_basis(layout, draw, grid_pitch_nm)
    raise ValueError(f"unknown electrostatics backend {backend!r}")
