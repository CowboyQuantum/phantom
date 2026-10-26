"""A scenario draw: material + catalog parameters + seed → a perturbed device, serialisable.

`draw_scenario` runs every generator with one seeded generator stream, in a fixed order, so the same
inputs give the same draw bit for bit (ADR 0003), and `severity = 0` gives a draw that leaves the layout
and the landscape unchanged. `apply_scenario` turns a draw into the perturbed layout and the perturbed
gate basis the charge layer derives from. The draw is a plain dict in `to_dict()`, so a consumer can
store it and `from_dict()` it back without phantom's classes.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field

import numpy as np

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
from phantom.disorder.wiring import WiringDraw, WiringFaultParams, draw_wiring_faults
from phantom.electrostatics.closed_form import GateBasis
from phantom.layout.geometry import Layout

MODEL_VERSION = 2


@dataclass(frozen=True)
class ScenarioParams:
    fixed_charges: FixedChargeParams = field(default_factory=FixedChargeParams)
    misalignment: MisalignmentParams = field(default_factory=MisalignmentParams)
    wiring: WiringFaultParams = field(default_factory=WiringFaultParams)

    @classmethod
    def for_material(cls, material: str) -> ScenarioParams:
        """Placeholder densities per material class from the catalog draft (`calibration: none`): MOS
        charges sit in the electron plane itself; SiGe charges at the oxide/cap interface above it."""
        if material == "MOS":
            return cls(FixedChargeParams(density_cm2=3e10, depth_nm=8.0), MisalignmentParams(), WiringFaultParams())
        return cls(FixedChargeParams(density_cm2=3e9, depth_nm=5.0), MisalignmentParams(), WiringFaultParams())


@dataclass(frozen=True)
class ScenarioDraw:
    seed: int
    severity: float
    charges: FixedChargeDraw
    misalignment: MisalignmentDraw
    wiring: WiringDraw | None = None
    model_version: int = MODEL_VERSION

    def to_dict(self) -> dict:
        return {
            "seed": self.seed,
            "severity": self.severity,
            "model_version": self.model_version,
            "charges": self.charges.to_dict(),
            "misalignment": self.misalignment.to_dict(),
            "wiring": self.wiring.to_dict() if self.wiring is not None else None,
        }

    @classmethod
    def from_dict(cls, d: dict) -> ScenarioDraw:
        wiring = WiringDraw.from_dict(d["wiring"]) if d.get("wiring") else None
        return cls(int(d["seed"]), float(d["severity"]), FixedChargeDraw.from_dict(d["charges"]), MisalignmentDraw.from_dict(d["misalignment"]), wiring, int(d.get("model_version", MODEL_VERSION)))


def adjacent_gate_pairs(layout: Layout) -> tuple[tuple[str, str], ...]:
    """Gate pairs whose polygons come within 30 nm of each other, by centroid distance along x within the
    same y band: the pairs that can short or leak."""
    pairs = []
    gates = layout.gates
    for i, a in enumerate(gates):
        for b in gates[i + 1 :]:
            _, ay = a.centroid_nm
            _, by = b.centroid_nm
            va, vb = a.vertices, b.vertices
            gap_x = max(va[:, 0].min(), vb[:, 0].min()) - min(va[:, 0].max(), vb[:, 0].max())
            same_band = abs(ay - by) < 0.5 * (np.ptp(va[:, 1]) + np.ptp(vb[:, 1]))
            if same_band and gap_x < 30.0:
                pairs.append((a.name, b.name))
    return tuple(pairs)


def draw_scenario(layout: Layout, params: ScenarioParams, seed: int, severity: float = 1.0, wiring_faults: bool = True) -> ScenarioDraw:
    rng = np.random.default_rng(seed)
    region = layout.bounding_box_nm()
    charges = draw_fixed_charges(params.fixed_charges, region, rng, severity)
    misalignment = draw_misalignment(params.misalignment, layout, rng, severity)
    wiring = draw_wiring_faults(params.wiring, layout.gate_names, adjacent_gate_pairs(layout), rng, severity if wiring_faults else 0.0)
    return ScenarioDraw(seed, severity, charges, misalignment, wiring)


def apply_scenario(layout: Layout, draw: ScenarioDraw, pitch_nm: float = 2.0) -> tuple[Layout, GateBasis]:
    """The perturbed layout and the gate basis with the charges' energy added as an offset. The basis is
    evaluated on the *clean* layout's grid, so draws of one layout share axes and the design's channel
    and sensor regions keep their meaning."""
    perturbed = apply_misalignment(layout, draw.misalignment)
    x0, x1, y0, y1 = layout.bounding_box_nm()
    grid = (np.arange(x0, x1 + pitch_nm / 2, pitch_nm), np.arange(y0, y1 + pitch_nm / 2, pitch_nm))
    basis = GateBasis.from_layout(perturbed, pitch_nm, grid=grid)
    if draw.charges.count == 0:
        return perturbed, basis
    offset = charge_energy_ev(draw.charges, basis.x_nm, basis.y_nm, layout.stack.effective_depth_nm, layout.stack.epsilon_well)
    return perturbed, dataclasses.replace(basis, offset_ev=offset)
