"""phantom.device — The compiled-device record: parameters only, serialisable, versioned, and the runtime objects every layer consumes, built from it without a solver (ADR 0001, 0023).

Sits between `spin` and `optimum`: `compile` writes a `Device`, `optimum` and consumers read one.
Modules: `session` (`State` at a gate point, `Session` with a clock: `sense`, `hold`, `wait`,
`spin_readout`), `truth` (`Truth`, `SensorTruth`, and `Knowable`, what a consumer may show a user).
The consumer surface over all of it is `phantom.api` (ADR 0024).
The electrostatics basis is not stored; `Device.basis()` recomputes it from the layout and the draw
(compile tier, deterministic), so a device file stays small and diffable.
"""

from __future__ import annotations

import json
import pathlib
from dataclasses import asdict, dataclass
from typing import Any

import numpy as np

from phantom.charge.ci import ConstantInteractionParams
from phantom.charge.dot import Barrier
from phantom.charge.hubbard import HubbardParams
from phantom.device.session import ChargeMode, Measurement, Session, SpinOutcome, State, state_at
from phantom.device.truth import Knowable, Truth, knowable, truth
from phantom.disorder.scenario import ScenarioDraw, apply_scenario
from phantom.electrostatics.closed_form import GateBasis
from phantom.layout.geometry import Gate, Layer, Layout, Stack
from phantom.readout.master_equation import RateEquation
from phantom.readout.measure import SensedDoubleDot
from phantom.readout.sensor import SensorSETParams
from phantom.spin.dynamics import SpinDraw

MODEL_VERSION = 2
"""2: the compiled sensor carries the wiring draw's sensor-coupling asymmetry (`compiled_sensor`, #43)."""
DEVICE_FORMAT_VERSION = 1
"""Bumped when the serialised layout of a `Device` changes incompatibly."""


def compiled_sensor(alpha_v_per_v: np.ndarray, beta_factors: tuple[float, float] = (1.0, 1.0)) -> SensorSETParams:
    """The SET sensor of a compiled device: its lever arms on the gates from the basis, and the nominal
    per-electron coupling to each dot scaled by the wiring draw's sensor-coupling asymmetry (ADR 0016,
    catalog #11). Compile places the set point with it and every runtime object reads through it."""
    beta = SensorSETParams().beta_v_per_e * np.asarray(beta_factors, dtype=float)
    return SensorSETParams(alpha_v_per_v=np.asarray(alpha_v_per_v, dtype=float), beta_v_per_e=beta)


@dataclass(frozen=True)
class CompiledDot:
    name: str
    x_nm: float
    y_nm: float
    energy_ev: float
    curvature_ev_per_nm2: float
    radius_nm: float
    lever_arms: tuple[float, ...]
    """Dimensionless potential per volt of every gate at the minimum."""
    valley_splitting_uev: float


@dataclass(frozen=True)
class CompiledBarrier:
    between: tuple[str, str]
    height_ev: float
    width_nm: float
    rate_hz: float


@dataclass(frozen=True)
class Device:
    """A compiled device: parameters only, serialisable, versioned."""

    seed: int
    material: str
    layout: dict[str, Any]
    """The design layout as data (`Layout` fields)."""
    scenario: dict[str, Any]
    """The scenario draw (`ScenarioDraw.to_dict`)."""
    programmed_point_v: tuple[float, ...]
    operating_point_v: tuple[float, ...]
    """The voltages the gates carry at the programmed point after wiring faults and the sensor set-point."""
    gate_names: tuple[str, ...]
    dots: tuple[CompiledDot, ...]
    sensor_dot: CompiledDot | None
    c_dot_af: tuple[tuple[float, ...], ...]
    c_gate_af: tuple[tuple[float, ...], ...]
    temperature_k: float
    sensor_alpha_v_per_v: tuple[float, ...]
    barriers: tuple[CompiledBarrier, ...]
    tunnel_coupling_ev: float
    """Interdot t at the operating point (estimate ħω₀√T/π)."""
    barrier_gate: str
    barrier_v_scale_v: float
    """Gate voltage that changes t by a factor e."""
    g_factors: tuple[float, ...]
    params: dict[str, Any]
    """The `CompileParams` used, as data."""
    provenance: dict[str, Any]
    """phantom version, format version and the model versions of every layer that produced this device."""

    # --- files -----------------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> Device:
        if int(d["provenance"]["format_version"]) != DEVICE_FORMAT_VERSION:
            raise ValueError(f"device format {d['provenance']['format_version']} != {DEVICE_FORMAT_VERSION}")

        def dot(x: dict | None) -> CompiledDot | None:
            return None if x is None else CompiledDot(**{**x, "lever_arms": tuple(x["lever_arms"])})

        return cls(
            seed=int(d["seed"]),
            material=d["material"],
            layout=d["layout"],
            scenario=d["scenario"],
            programmed_point_v=tuple(d["programmed_point_v"]),
            operating_point_v=tuple(d["operating_point_v"]),
            gate_names=tuple(d["gate_names"]),
            dots=tuple(dot(x) for x in d["dots"]),
            sensor_dot=dot(d["sensor_dot"]),
            c_dot_af=tuple(tuple(r) for r in d["c_dot_af"]),
            c_gate_af=tuple(tuple(r) for r in d["c_gate_af"]),
            temperature_k=float(d["temperature_k"]),
            sensor_alpha_v_per_v=tuple(d["sensor_alpha_v_per_v"]),
            barriers=tuple(CompiledBarrier(tuple(b["between"]), b["height_ev"], b["width_nm"], b["rate_hz"]) for b in d["barriers"]),
            tunnel_coupling_ev=float(d["tunnel_coupling_ev"]),
            barrier_gate=d["barrier_gate"],
            barrier_v_scale_v=float(d["barrier_v_scale_v"]),
            g_factors=tuple(d["g_factors"]),
            params=d["params"],
            provenance=d["provenance"],
        )

    def to_json(self, path: str | pathlib.Path) -> None:
        pathlib.Path(path).write_text(self.dumps(), encoding="utf-8")

    def dumps(self) -> str:
        """Canonical JSON: sorted keys, one entry per line, floats at full precision (diffable)."""
        return json.dumps(self.to_dict(), indent=1, sort_keys=True)

    @classmethod
    def from_json(cls, path: str | pathlib.Path) -> Device:
        return cls.from_dict(json.loads(pathlib.Path(path).read_text(encoding="utf-8")))

    @classmethod
    def loads(cls, text: str) -> Device:
        return cls.from_dict(json.loads(text))

    # --- runtime objects (never a solver) --------------------------------------------------------

    def layout_object(self) -> Layout:
        d = self.layout
        gates = tuple(Gate(g["name"], tuple(tuple(p) for p in g["polygon_nm"]), g.get("role", "gate")) for g in d["gates"])
        stack = Stack(tuple(Layer(**layer) for layer in d["stack"]["layers"]), d["stack"]["epsilon_well"])
        sensor = tuple(d["sensor_nm"]) if d.get("sensor_nm") else None
        return Layout(gates, stack, d.get("material", "SiGe"), tuple(d["channel_nm"]), sensor, d.get("margin_nm", 60.0))

    def scenario_draw(self) -> ScenarioDraw:
        return ScenarioDraw.from_dict(self.scenario)

    def basis(self) -> tuple[Layout, GateBasis]:
        """Recompute the perturbed layout and the closed-form basis (compile tier; deterministic for
        the draw). A device compiled with another backend records it in `params["backend"]`; its basis
        is rebuilt by `compile.backends.basis_for`."""
        if self.params.get("backend", "closed_form") != "closed_form":
            raise ValueError("rebuild a non-closed-form basis with compile.backends.basis_for")
        return apply_scenario(self.layout_object(), self.scenario_draw(), float(self.params["grid_pitch_nm"]))

    def operating_point(self) -> np.ndarray:
        return np.array(self.operating_point_v, dtype=float)

    def effective_voltages(self, v_programmed: np.ndarray) -> np.ndarray:
        """Programmed (..., gates) → carried, through the wiring draw."""
        draw = self.scenario_draw()
        return draw.wiring.apply(v_programmed) if draw.wiring is not None else np.asarray(v_programmed, dtype=float)

    def runtime_ci(self) -> ConstantInteractionParams:
        return ConstantInteractionParams(gate_names=self.gate_names, c_dot_af=np.array(self.c_dot_af), c_gate_af=np.array(self.c_gate_af), temperature_k=self.temperature_k)

    def sensor_beta_factors(self) -> tuple[float, float]:
        """The wiring draw's sensor-coupling factors per dot ((1, 1) without a wiring draw)."""
        wiring = self.scenario_draw().wiring
        return tuple(wiring.sensor_beta_factors) if wiring is not None else (1.0, 1.0)

    def runtime_sensor(self) -> SensorSETParams:
        return compiled_sensor(np.array(self.sensor_alpha_v_per_v), self.sensor_beta_factors())

    def runtime_barrier(self) -> Barrier:
        j = self.gate_names.index(self.barrier_gate)
        return Barrier(self.barrier_gate, t0_ev=self.tunnel_coupling_ev, v_ref_v=self.operating_point_v[j], v_scale_v=self.barrier_v_scale_v)

    def runtime_hubbard(self) -> HubbardParams:
        return HubbardParams.from_ci(self.runtime_ci(), self.runtime_barrier(), barrier_gate=self.barrier_gate)

    def runtime_rates(self) -> RateEquation:
        rates = {b.between: b.rate_hz for b in self.barriers}
        return RateEquation(self.runtime_ci(), rates[("dot1", "lead_left")], rates[("dot2", "lead_right")], rates[("dot1", "dot2")])

    def runtime_sensed(self, sample_dt_s: float = 1e-3, noise: bool = True) -> SensedDoubleDot:
        draw = self.scenario_draw()
        if noise and draw.noise is not None:
            return SensedDoubleDot.from_noise_draw(draw.noise, sample_dt_s, charge=self.runtime_ci(), sensor=self.runtime_sensor())
        return SensedDoubleDot(charge=self.runtime_ci(), sensor=self.runtime_sensor())

    def runtime_spin(self) -> SpinDraw:
        return SpinDraw(self.g_factors, self.seed)

    def valley_splittings_uev(self) -> np.ndarray:
        return np.array([d.valley_splitting_uev for d in self.dots])

    def noise_sigma_gate_v(self) -> float:
        draw = self.scenario_draw()
        return float(draw.noise.one_over_f_amplitude_v) if draw.noise is not None else 0.0

    # --- the consumer surface (plan §4, ADR 0024) ----------------------------------------------------

    def at(self, v_programmed: np.ndarray, v_sd: float | np.ndarray | None = None) -> State:
        """The noise-free `State` at programmed voltages (..., gates); `v_sd` is the sensor bias (V)."""
        return state_at(self, v_programmed, v_sd)

    def session(
        self,
        noise: bool = True,
        dynamics: bool = True,
        charge: ChargeMode = "sampled",
        noise_seed: int | None = None,
    ) -> Session:
        """A fresh run with its own clock (`device.session`): `noise_seed=None` replays the device's
        noise realisation; an integer gives a new one with the same noise parameters."""
        return Session(self, noise=noise, dynamics=dynamics, charge=charge, noise_seed=noise_seed)

    def measure(
        self,
        v_programmed: np.ndarray,
        dwell_s: float,
        samples: int | None = None,
        v_sd: float | np.ndarray | None = None,
        noise: bool = True,
        dynamics: bool = True,
        charge: ChargeMode = "sampled",
        noise_seed: int | None = None,
    ) -> Measurement:
        """`sense` on a fresh session."""
        return self.session(noise, dynamics, charge, noise_seed).sense(v_programmed, dwell_s, samples, v_sd)

    def truth(
        self,
        lattice_gates: tuple[str, str] = ("P1", "P2"),
        near_v: float | None = None,
        v_sd: float | None = None,
        v_programmed: np.ndarray | None = None,
    ) -> Truth:
        """The hidden parameters under the quantities a tune-up recovers (`device.truth`)."""
        return truth(self, lattice_gates, near_v, v_sd, v_programmed)

    def knowable(self) -> Knowable:
        """The fields a consumer may show a user; everything else on the device is truth."""
        return knowable(self)



__all__ = [
    "DEVICE_FORMAT_VERSION",
    "MODEL_VERSION",
    "CompiledBarrier",
    "CompiledDot",
    "Device",
    "Knowable",
    "Measurement",
    "Session",
    "SpinOutcome",
    "State",
    "Truth",
    "compiled_sensor",
]
