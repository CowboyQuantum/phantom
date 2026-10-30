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
from phantom.constants import K_B_EV_PER_K
from phantom.device.session import ChargeMode, Measurement, Session, SpinOutcome, State, state_at
from phantom.device.truth import Knowable, Truth, knowable, truth
from phantom.disorder.scenario import ScenarioDraw, apply_scenario
from phantom.electrostatics.closed_form import GateBasis
from phantom.layout.geometry import Gate, Layer, Layout, Stack
from phantom.readout.master_equation import RateEquation
from phantom.readout.measure import SensedDoubleDot
from phantom.readout.rates import GateRates
from phantom.readout.sensor import SensorLeads, SensorSETParams
from phantom.spin.dynamics import SpinDraw

MODEL_VERSION = 4
"""2: the compiled sensor carries the wiring draw's sensor-coupling asymmetry (`compiled_sensor`, #43).
3: a sensor with leads is a dot like the qubit dots, read through its barriers and baths (ADR 0026).
4: the dots' lead and interdot rates follow the gates at runtime (`runtime_rates`, #53, ADR 0027)."""
DEVICE_FORMAT_VERSION = 3
"""Bumped when the serialised layout of a `Device` changes incompatibly. 2: the sensor's charging
energy, per-electron coupling and leads (ADR 0026). 3: each barrier's attempt rate, tunnelling level
and gate slopes (#53)."""


def compiled_sensor(
    alpha_v_per_v: np.ndarray,
    beta_factors: tuple[float, float] = (1.0, 1.0),
    charging_energy_ev: float | None = None,
    beta_v_per_e: tuple[float, float] | None = None,
    temperature_k: float | None = None,
    leads: SensorLeads | None = None,
) -> SensorSETParams:
    """The SET sensor of a compiled device: its lever arms on the gates from the basis, and its
    per-electron coupling to each dot scaled by the wiring draw's sensor-coupling asymmetry (ADR 0016,
    catalog #11). Compile places the set point with it and every runtime object reads through it.

    With `leads` (a layout that gives the sensor barriers and baths) the sensor is a dot like the qubit
    dots (ADR 0026): φ_s is the potential at the dot in volts (lever arm 1 eV/V), its peaks are its
    charging energy apart, its thermal width is the device's 2kT, its coupling β the layout's, and its
    current comes through the leads. Without leads it is the typed Step 1 SET with the typed β."""
    alpha = np.asarray(alpha_v_per_v, dtype=float)
    factors = np.asarray(beta_factors, dtype=float)
    if leads is None or charging_energy_ev is None or beta_v_per_e is None or temperature_k is None:
        return SensorSETParams(alpha_v_per_v=alpha, beta_v_per_e=SensorSETParams().beta_v_per_e * factors)
    e_c = float(charging_energy_ev)
    return SensorSETParams(
        alpha_v_per_v=alpha,
        beta_v_per_e=np.asarray(beta_v_per_e, dtype=float) * factors,
        period_v=e_c,
        width_v=2 * K_B_EV_PER_K * float(temperature_k),
        charging_energy_ev=e_c,
        leads=leads,
    )


def switcher_scale(beta_v_per_e: tuple[float, float] | None) -> float:
    """Factor from the noise draw's switcher amplitudes, stated in the typed sensor's V-equivalent
    units, to a sensor dot's potential in volts: the ratio of the two sensors' mean per-electron
    coupling, so that a switcher keeps its size against a dot electron's signal (ADR 0026). 1 for the
    typed sensor."""
    if beta_v_per_e is None:
        return 1.0
    return float(np.mean(np.abs(beta_v_per_e)) / np.mean(np.abs(SensorSETParams().beta_v_per_e)))


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
    """f_attempt·exp(−2S) at the operating point, Hz."""
    attempt_rate_hz: float
    """Attempt frequency of the dot the barrier is seen from, Hz: the ceiling of the rate at runtime."""
    tunnelling_level_ev: float
    """The level the rate is taken at, placed at the operating point and held fixed, eV."""
    log_rate_slopes_per_v: tuple[float, ...]
    """d ln Γ/dV against every gate at that level, 1/V: the runtime rate follows the gates (#53)."""


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
    sensor_charging_energy_ev: float | None = None
    """The sensor dot's charging energy e²(C⁻¹)_ss, eV (ADR 0026); None without a sensor dot."""
    sensor_beta_v_per_e: tuple[float, float] | None = None
    """Shift of the sensor potential per electron on (dot 1, dot 2) from the layout, V, before the
    wiring draw's asymmetry; None without a sensor dot."""
    sensor_leads: SensorLeads | None = None
    """The sensor's barriers and baths compiled at the operating point; None when the layout gives the
    sensor no reservoirs (the typed sensor then reads it)."""

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
            barriers=tuple(CompiledBarrier(**{**b, "between": tuple(b["between"]), "log_rate_slopes_per_v": tuple(b["log_rate_slopes_per_v"])}) for b in d["barriers"]),
            tunnel_coupling_ev=float(d["tunnel_coupling_ev"]),
            barrier_gate=d["barrier_gate"],
            barrier_v_scale_v=float(d["barrier_v_scale_v"]),
            g_factors=tuple(d["g_factors"]),
            params=d["params"],
            provenance=d["provenance"],
            sensor_charging_energy_ev=None if d["sensor_charging_energy_ev"] is None else float(d["sensor_charging_energy_ev"]),
            sensor_beta_v_per_e=None if d["sensor_beta_v_per_e"] is None else tuple(float(x) for x in d["sensor_beta_v_per_e"]),
            sensor_leads=None if d["sensor_leads"] is None else SensorLeads.from_dict(d["sensor_leads"]),
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
        return compiled_sensor(
            np.array(self.sensor_alpha_v_per_v),
            self.sensor_beta_factors(),
            self.sensor_charging_energy_ev,
            self.sensor_beta_v_per_e,
            self.temperature_k,
            self.sensor_leads,
        )

    def switcher_scale(self) -> float:
        """The noise draw's switcher amplitudes in this device's sensor units (`switcher_scale`)."""
        return switcher_scale(self.sensor_beta_v_per_e if self.sensor_leads is not None else None)

    def runtime_barrier(self) -> Barrier:
        j = self.gate_names.index(self.barrier_gate)
        return Barrier(self.barrier_gate, t0_ev=self.tunnel_coupling_ev, v_ref_v=self.operating_point_v[j], v_scale_v=self.barrier_v_scale_v)

    def runtime_hubbard(self) -> HubbardParams:
        return HubbardParams.from_ci(self.runtime_ci(), self.runtime_barrier(), barrier_gate=self.barrier_gate)

    def runtime_rates(self) -> RateEquation:
        """The dots' rate equation, its lead and interdot rates following the carried gate voltages
        from their values at the operating point (#53)."""
        by = {b.between: b for b in self.barriers}
        left, right, inter = (by[k] for k in (("dot1", "lead_left"), ("dot2", "lead_right"), ("dot1", "dot2")))
        gated = GateRates(
            v_ref_v=tuple(self.operating_point_v),
            gamma_ref_hz=(left.rate_hz, right.rate_hz, inter.rate_hz),
            attempt_rate_hz=(left.attempt_rate_hz, right.attempt_rate_hz, inter.attempt_rate_hz),
            log_rate_slopes_per_v=(left.log_rate_slopes_per_v, right.log_rate_slopes_per_v, inter.log_rate_slopes_per_v),
        )
        return RateEquation(self.runtime_ci(), left.rate_hz, right.rate_hz, inter.rate_hz, gated)

    def runtime_sensed(self, sample_dt_s: float = 1e-3, noise: bool = True) -> SensedDoubleDot:
        draw = self.scenario_draw()
        if noise and draw.noise is not None:
            return SensedDoubleDot.from_noise_draw(draw.noise, sample_dt_s, self.switcher_scale(), charge=self.runtime_ci(), sensor=self.runtime_sensor())
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
