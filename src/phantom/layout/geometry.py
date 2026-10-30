"""Gate polygons, the dielectric stack and the material class: data only, no physics.

Lengths in nm, voltages in V, permittivities relative to vacuum. A `Layout` is serialisable to JSON
(standard library only; rule 1 allows no other runtime dependency), so "a layout file is the only
input needed" (Step 2 exit criterion) means a JSON file read by `Layout.from_json`.
"""

from __future__ import annotations

import json
import pathlib
from dataclasses import asdict, dataclass

import numpy as np

MODEL_VERSION = 2
"""Bumped whenever any output derived from these records changes for the same inputs (ADR 0003).
2: `two_dot_layout` gains the sensor's barriers and bath (SB_L, SB_R, SR_L, SR_R; ADR 0026)."""


@dataclass(frozen=True)
class Gate:
    """A metal gate as a planar polygon at the surface, held at a voltage."""

    name: str
    polygon_nm: tuple[tuple[float, float], ...]
    """Vertices (x, y) in nm, in order around the polygon; convex or not, no self-intersection."""
    role: str = "gate"
    """`plunger`, `barrier`, `sensor`, `reservoir`, `screening` or `gate`; informs defaults, not physics."""

    def __post_init__(self) -> None:
        if len(self.polygon_nm) < 3:
            raise ValueError(f"gate {self.name!r} needs at least three vertices")

    @property
    def vertices(self) -> np.ndarray:
        return np.asarray(self.polygon_nm, dtype=float)

    @property
    def area_nm2(self) -> float:
        """Polygon area by the shoelace formula."""
        v = self.vertices
        x, y = v[:, 0], v[:, 1]
        return float(0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))

    @property
    def centroid_nm(self) -> tuple[float, float]:
        c = self.vertices.mean(axis=0)
        return float(c[0]), float(c[1])

    @classmethod
    def rectangle(cls, name: str, x0: float, x1: float, y0: float, y1: float, role: str = "gate") -> Gate:
        return cls(name, ((x0, y0), (x1, y0), (x1, y1), (x0, y1)), role)


@dataclass(frozen=True)
class Layer:
    """One dielectric layer between the gate plane and the electron plane."""

    name: str
    thickness_nm: float
    epsilon_r: float

    def __post_init__(self) -> None:
        if self.thickness_nm <= 0 or self.epsilon_r <= 0:
            raise ValueError(f"layer {self.name!r}: thickness and permittivity must be positive")


@dataclass(frozen=True)
class Stack:
    """Layers from the gate plane down to the plane of the two-dimensional electron gas (2DEG).

    The closed-form potential (electrostatics) assumes one uniform dielectric. A layered stack is
    mapped onto it by an *effective depth*: each layer counts with its thickness scaled by
    epsilon_well / epsilon_layer, so a thin low-permittivity oxide acts like a thicker layer of the well
    material. This is the standard thin-film approximation and is exact for the normal field far from
    gate edges; Step 7's solver replaces it (ADR 0013).
    """

    layers: tuple[Layer, ...]
    epsilon_well: float = 11.7
    """Relative permittivity of the layer holding the electrons (silicon: 11.7)."""

    def __post_init__(self) -> None:
        if not self.layers:
            raise ValueError("a stack needs at least one layer")

    @property
    def depth_nm(self) -> float:
        return float(sum(layer.thickness_nm for layer in self.layers))

    @property
    def effective_depth_nm(self) -> float:
        return float(sum(layer.thickness_nm * self.epsilon_well / layer.epsilon_r for layer in self.layers))


@dataclass(frozen=True)
class MaterialClass:
    """Material-class defaults (ADR 0002: a class, not a device). Values are textbook ranges, not a
    measured device; a reference entry replaces each default as it is confirmed (ADR 0007)."""

    name: str
    effective_mass_me: float
    """Electron effective mass in the plane, in units of the free-electron mass (silicon valleys: 0.19)."""
    epsilon_well: float
    default_stack: Stack

    @classmethod
    def sige(cls) -> MaterialClass:
        """Si/SiGe: a silicon quantum well under a SiGe spacer and a thin silicon cap."""
        return cls(
            "SiGe",
            effective_mass_me=0.19,
            epsilon_well=11.7,
            default_stack=Stack(
                (Layer("oxide", 5.0, 3.9), Layer("Si cap", 2.0, 11.7), Layer("SiGe spacer", 30.0, 13.0)),
                epsilon_well=11.7,
            ),
        )

    @classmethod
    def mos(cls) -> MaterialClass:
        """Si MOS: electrons at the silicon–oxide interface under a thin oxide."""
        return cls(
            "MOS",
            effective_mass_me=0.19,
            epsilon_well=11.7,
            default_stack=Stack((Layer("oxide", 8.0, 3.9),), epsilon_well=11.7),
        )

    @classmethod
    def by_name(cls, name: str) -> MaterialClass:
        return {"SiGe": cls.sige, "MOS": cls.mos}[name]()


@dataclass(frozen=True)
class Layout:
    """Gates, stack, material and the regions where dots are looked for."""

    gates: tuple[Gate, ...]
    stack: Stack
    material: str = "SiGe"
    channel_nm: tuple[float, float, float, float] = (-150.0, 150.0, -40.0, 40.0)
    """(x0, x1, y0, y1): where the double dot forms; minima of the potential energy are searched here."""
    sensor_nm: tuple[float, float, float, float] | None = None
    """Where the sensor dot forms, if the layout has one."""
    margin_nm: float = 60.0
    """Grid margin around the gates' bounding box."""

    def __post_init__(self) -> None:
        names = [g.name for g in self.gates]
        if len(set(names)) != len(names):
            raise ValueError(f"gate names must be unique: {names}")

    @property
    def gate_names(self) -> tuple[str, ...]:
        return tuple(g.name for g in self.gates)

    def reservoir_points_nm(self) -> dict[str, tuple[float, float]]:
        """Centroids of the reservoir gates, by name: where the leads are for tunnel-rate estimates."""
        return {g.name: g.centroid_nm for g in self.gates if g.role == "reservoir"}

    def gate(self, name: str) -> Gate:
        return self.gates[self.gate_names.index(name)]

    def bounding_box_nm(self) -> tuple[float, float, float, float]:
        v = np.concatenate([g.vertices for g in self.gates])
        m = self.margin_nm
        return float(v[:, 0].min() - m), float(v[:, 0].max() + m), float(v[:, 1].min() - m), float(v[:, 1].max() + m)

    # --- files -----------------------------------------------------------------------------------

    def to_json(self, path: str | pathlib.Path) -> None:
        pathlib.Path(path).write_text(json.dumps(asdict(self), indent=1), encoding="utf-8")

    @classmethod
    def from_json(cls, path: str | pathlib.Path) -> Layout:
        d = json.loads(pathlib.Path(path).read_text(encoding="utf-8"))
        gates = tuple(Gate(g["name"], tuple(tuple(p) for p in g["polygon_nm"]), g.get("role", "gate")) for g in d["gates"])
        stack = Stack(tuple(Layer(**layer) for layer in d["stack"]["layers"]), d["stack"]["epsilon_well"])
        sensor = tuple(d["sensor_nm"]) if d.get("sensor_nm") else None
        return cls(gates, stack, d.get("material", "SiGe"), tuple(d["channel_nm"]), sensor, d.get("margin_nm", 60.0))


def two_dot_layout(material: str = "SiGe", pitch_nm: float = 80.0) -> Layout:
    """A textbook single-layer double-dot layout with a charge sensor: two plungers P1, P2 over a channel
    along x, three barriers B_L, B_M, B_R across it, and two reservoir accumulation gates R_L, R_R beyond
    the outer barriers that hold the leads populated, so the lead barriers are the outer barrier gates and
    not the ungated region. The sensor is a dot of the same kind in a parallel channel at positive y: a
    sensor plunger SG between sensor barriers SB_L and SB_R, and the sensor's bath, two reservoir gates
    SR_L and SR_R beyond them that accumulate its source and drain (ADR 0026). The sensor gates follow the
    qubit channel's pattern (plunger and barrier widths, gaps) and are listed last, so the first eight
    gates keep their order."""
    half = pitch_nm / 2
    pw, bw, gl = 0.55 * pitch_nm, 0.3 * pitch_nm, 120.0  # plunger width, barrier width, gate length (y)
    gap = (pitch_nm - pw - bw) / 2  # between a plunger and the next barrier
    ys = 150.0  # the sensor gates' lower edge
    gates = (
        Gate.rectangle("B_L", -half - pitch_nm / 2 - bw / 2, -half - pitch_nm / 2 + bw / 2, -gl / 2, gl / 2, "barrier"),
        Gate.rectangle("P1", -half - pw / 2, -half + pw / 2, -gl / 2, gl / 2, "plunger"),
        Gate.rectangle("B_M", -bw / 2, bw / 2, -gl / 2, gl / 2, "barrier"),
        Gate.rectangle("P2", half - pw / 2, half + pw / 2, -gl / 2, gl / 2, "plunger"),
        Gate.rectangle("B_R", half + pitch_nm / 2 - bw / 2, half + pitch_nm / 2 + bw / 2, -gl / 2, gl / 2, "barrier"),
        Gate.rectangle("SG", -pw / 2, pw / 2, ys, ys + gl, "sensor"),
        Gate.rectangle("R_L", -4 * pitch_nm, -half - pitch_nm / 2 - bw / 2 - 10.0, -gl / 2, gl / 2, "reservoir"),
        Gate.rectangle("R_R", half + pitch_nm / 2 + bw / 2 + 10.0, 4 * pitch_nm, -gl / 2, gl / 2, "reservoir"),
        Gate.rectangle("SB_L", -pw / 2 - gap - bw, -pw / 2 - gap, ys, ys + gl, "barrier"),
        Gate.rectangle("SB_R", pw / 2 + gap, pw / 2 + gap + bw, ys, ys + gl, "barrier"),
        Gate.rectangle("SR_L", -4 * pitch_nm, -pw / 2 - gap - bw - 10.0, ys, ys + gl, "reservoir"),
        Gate.rectangle("SR_R", pw / 2 + gap + bw + 10.0, 4 * pitch_nm, ys, ys + gl, "reservoir"),
    )
    mat = MaterialClass.by_name(material)
    return Layout(
        gates,
        mat.default_stack,
        material,
        channel_nm=(-2 * pitch_nm, 2 * pitch_nm, -0.3 * gl, 0.3 * gl),
        sensor_nm=(-pw / 2 - gap - bw, pw / 2 + gap + bw, ys + 0.2 * gl, ys + 0.8 * gl),
    )
