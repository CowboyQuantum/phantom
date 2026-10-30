"""What a compiled device knows about itself, and what a consumer may show a user (ADR 0003, 0024).

`Knowable` is the part of a `Device` that a person running the device would know before measuring:
the material, the layout and its gate names and roles, the electron temperature and the programmed
operating point. Everything else on a `Device` is truth: the seed (recompiling from it reproduces the
draw), the scenario draw, the dots, capacitances and rates, the valley and spin draws, `params`, and
`operating_point_v`, which holds the voltages the gates carry after the wiring draw and the sensor
set-point compilation chose. A consumer that hides truth forwards `knowable()` and the device's
measurements, and nothing else.

`Truth` gathers the hidden parameters under the quantities a tune-up recovers. Two are derived here:

- `charge_lattice_e_per_v`: electrons per volt on (dot 1, dot 2) from two programmed gates, the rows
  of C_gate·M / e with M the wiring matrix (`charge.ci.charge_lattice` when the wiring is clean);
- `sensor` (`SensorTruth`): the sensor's Coulomb peak in programmed sensor-gate volts, with the
  other gates at a reference point and the dot charge held at its thermal occupation there. For the
  typed sensor a peak in linear response is I(V) = h · sech²((V − V_peak)/w) above the background,
  h = peak_a · V_sd / bias_v, and the wing is its inflection, w · atanh(1/√3) above the peak, where
  |dI/dV| = 4h / (3√3 w). A sensor dot with leads (ADR 0026) has the lifetime-broadened lineshape of
  `readout.broadening`: h, the wing (the steepest point) and the slope there are those of the thermal
  peak convolved with a Lorentzian of full width ħ(Γ_S + Γ_D) (`lifetime_width_v`), and reduce to the
  sech² statements when that width is small against w. These are qdt's `sensor-operating-point` and
  `sensor-wing` quantities (qdt ADR 0062), stated by phantom.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from phantom.constants import E_UNIT_V_PER_E_AF

if TYPE_CHECKING:
    from phantom.device import CompiledDot, Device

MODEL_VERSION = 2
"""2: the sensor's leads in `SensorTruth` and `Knowable`; the wing and slope of a broadened peak (ADR 0026)."""


@dataclass(frozen=True)
class Knowable:
    """What a consumer may show a user: known before any measurement."""

    material: str
    gate_names: tuple[str, ...]
    layout: dict[str, Any]
    """The design layout as data (`Layout` fields): gate polygons, roles, stack."""
    temperature_k: float
    """Electron temperature, K."""
    programmed_point_v: tuple[float, ...]
    """The voltages the consumer programs, one per gate, V."""
    sensor_gate: str
    barrier_gate: str
    reservoir_gates: tuple[str, ...]
    sensor_barrier_gates: tuple[str, ...] = ()
    """The sensor's barrier gates, source side first (empty when the layout has none)."""
    sensor_reservoir_gates: tuple[str, ...] = ()
    """The sensor's bath gates, source side first (empty when the layout has none)."""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class SensorTruth:
    """The sensor's Coulomb peak against the programmed sensor-gate voltage."""

    gate: str
    peak_v: float
    """Programmed sensor-gate voltage of the peak nearest the requested voltage, V."""
    period_v: float
    """Sensor-gate volts between peaks, V."""
    charging_energy_ev: float
    """Charging energy of the sensor SET, eV."""
    lever_arm_ev_per_v: float
    """Sensor level shift per programmed sensor-gate volt, eV/V (charging energy / period)."""
    width_v: float
    """Thermal width w of a peak in sensor-gate volts, V (I ∝ sech²((V − V_peak)/w) without lifetime broadening)."""
    wing_v: float
    """The wing, where the linear-response slope is steepest above the peak, V: V_peak + w·atanh(1/√3)
    without lifetime broadening, further out with it."""
    sensitivity_a_per_v: float
    """|dI/dV| at the wing at bias `v_sd_v`, A/V (linear response)."""
    v_sd_v: float
    """Source–drain bias the sensitivity is stated at, V."""
    lifetime_width_v: float = 0.0
    """ħ(Γ_S + Γ_D) in sensor-gate volts: the Lorentzian full width of a level, V; 0 for the typed sensor."""
    gamma_source_hz: float | None = None
    """Tunnel rate through the sensor's source barrier into a full bath, Hz; None for the typed sensor."""
    gamma_drain_hz: float | None = None
    """Tunnel rate through the sensor's drain barrier into a full bath, Hz; None for the typed sensor."""


@dataclass(frozen=True)
class Truth:
    """The hidden parameters of a compiled device (`Device.truth`)."""

    seed: int
    material: str
    scenario: dict[str, Any]
    """The scenario draw: interface charges, misalignment, wiring faults, noise draw."""
    operating_point_v: tuple[float, ...]
    """Carried voltages at the programmed point, with the compiled sensor set-point, V."""
    dots: tuple[CompiledDot, ...]
    sensor_dot: CompiledDot | None
    addition_energies_ev: tuple[float, float]
    mutual_charging_energy_ev: float
    lever_arms_ev_per_v: tuple[tuple[float, ...], ...]
    """[dot, gate] on carried voltages, eV/V."""
    lattice_gates: tuple[str, str]
    charge_lattice_e_per_v: tuple[tuple[float, float], tuple[float, float]]
    """Electrons per programmed volt on (dot 1, dot 2) [row] from `lattice_gates` [column]."""
    tunnel_coupling_ev: float
    """Interdot t at the operating point, eV."""
    barrier_gate: str
    barrier_v_scale_v: float
    """Barrier-gate voltage that changes t by a factor e, V."""
    gamma_left_hz: float
    """Dot 1 to the left lead at the operating point, Hz; `Device.at` gives it at any point (#53)."""
    gamma_right_hz: float
    """Dot 2 to the right lead at the operating point, Hz."""
    gamma_interdot_hz: float
    """Electron transfer between the dots at the operating point, Hz."""
    valley_splittings_uev: tuple[float, ...]
    g_factors: tuple[float, ...]
    sensor: SensorTruth | None
    """None when the sensor gate does not move the sensor (an open or missing sensor gate)."""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def knowable(device: Device) -> Knowable:
    p = device.params
    leads = device.sensor_leads
    return Knowable(
        material=device.material,
        gate_names=tuple(device.gate_names),
        layout=device.layout,
        temperature_k=float(device.temperature_k),
        programmed_point_v=tuple(device.programmed_point_v),
        sensor_gate=str(p.get("sensor_gate", "SG")),
        barrier_gate=device.barrier_gate,
        reservoir_gates=tuple(p.get("reservoir_gates", ())),
        sensor_barrier_gates=tuple(leads.barrier_gates) if leads is not None else (),
        sensor_reservoir_gates=tuple(leads.reservoir_gates) if leads is not None else (),
    )


def _wiring_matrix(device: Device) -> np.ndarray:
    wiring = device.scenario_draw().wiring
    return wiring.wiring_matrix()[0] if wiring is not None else np.eye(len(device.gate_names))


def sensor_truth(device: Device, near_v: float | None = None, v_sd: float | None = None, v_programmed: np.ndarray | None = None) -> SensorTruth | None:
    """The sensor peak nearest programmed sensor-gate voltage `near_v` (the reference point's when
    None), with the other gates at programmed `v_programmed` (the device's programmed point when None)
    and the dot charge held at its thermal occupation there; the sensitivity at bias `v_sd` (V; the
    sensor's default bias when None)."""
    gate = str(device.params.get("sensor_gate", "SG"))
    if gate not in device.gate_names:
        return None
    j = device.gate_names.index(gate)
    sensor = device.runtime_sensor()
    a = float((sensor.alpha_v_per_v @ _wiring_matrix(device))[j])
    if not a > 0:
        return None
    v_ref = np.array(device.programmed_point_v if v_programmed is None else v_programmed, dtype=float)
    vc = device.effective_voltages(v_ref)
    phi0 = float(sensor.potential(vc, device.runtime_ci().mean_occupation(vc)))
    near = float(v_ref[j]) if near_v is None else float(near_v)
    k = np.rint((phi0 + a * (near - v_ref[j]) - sensor.phase_v) / sensor.period_v)
    peak = float(v_ref[j] + (sensor.phase_v + k * sensor.period_v - phi0) / a)
    period, width = sensor.period_v / a, sensor.width_v / a
    bias = sensor.bias_v if v_sd is None else float(v_sd)
    # the leads' rates at the peak: the reference point with the sensor gate moved there
    vp = v_ref.copy()
    vp[j] = peak
    v_peak = device.effective_voltages(vp)
    rates = sensor.rates_hz(v_peak)
    return SensorTruth(
        gate=gate,
        peak_v=peak,
        period_v=float(period),
        charging_energy_ev=float(sensor.charging_energy_ev),
        lever_arm_ev_per_v=float(sensor.charging_energy_ev / period),
        width_v=float(width),
        wing_v=float(peak + sensor.wing_offset_v(v_peak) / a),
        sensitivity_a_per_v=float(sensor.wing_slope_a_per_v(bias, v_peak) * a),
        v_sd_v=float(bias),
        lifetime_width_v=float(sensor.lifetime_width_v(v_peak) / a),
        gamma_source_hz=None if rates is None else float(rates[0]),
        gamma_drain_hz=None if rates is None else float(rates[1]),
    )


def truth(
    device: Device,
    lattice_gates: tuple[str, str] = ("P1", "P2"),
    near_v: float | None = None,
    v_sd: float | None = None,
    v_programmed: np.ndarray | None = None,
) -> Truth:
    """The full record; `near_v`, `v_sd` and `v_programmed` as in `sensor_truth`."""
    ci = device.runtime_ci()
    i, j = (device.gate_names.index(g) for g in lattice_gates)
    lattice = (np.asarray(ci.c_gate_af) @ _wiring_matrix(device))[:, [i, j]] / E_UNIT_V_PER_E_AF
    rates = device.runtime_rates()
    add = ci.addition_energies_ev()
    return Truth(
        seed=int(device.seed),
        material=device.material,
        scenario=device.scenario,
        operating_point_v=tuple(device.operating_point_v),
        dots=tuple(device.dots),
        sensor_dot=device.sensor_dot,
        addition_energies_ev=(float(add[0]), float(add[1])),
        mutual_charging_energy_ev=ci.mutual_charging_energy_ev(),
        lever_arms_ev_per_v=tuple(tuple(float(x) for x in row) for row in ci.lever_arms_ev_per_v()),
        lattice_gates=(lattice_gates[0], lattice_gates[1]),
        charge_lattice_e_per_v=((float(lattice[0, 0]), float(lattice[0, 1])), (float(lattice[1, 0]), float(lattice[1, 1]))),
        tunnel_coupling_ev=float(device.tunnel_coupling_ev),
        barrier_gate=device.barrier_gate,
        barrier_v_scale_v=float(device.barrier_v_scale_v),
        gamma_left_hz=float(rates.gamma_left_hz),
        gamma_right_hz=float(rates.gamma_right_hz),
        gamma_interdot_hz=float(rates.gamma_interdot_hz),
        valley_splittings_uev=tuple(float(x) for x in device.valley_splittings_uev()),
        g_factors=tuple(float(g) for g in device.g_factors),
        sensor=sensor_truth(device, near_v, v_sd, v_programmed),
    )


__all__ = ["MODEL_VERSION", "Knowable", "SensorTruth", "Truth", "knowable", "sensor_truth", "truth"]
