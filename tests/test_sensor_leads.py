"""P0 for the sensor dot's barriers and bath (ADR 0026, #47): the lifetime-broadened lineshape against
its closed-form limits, the compiled leads against the compile tier, and the sensor dot of the textbook
layout as a dot like the qubit dots. Analytic limits and invariants only; no claim of agreement with a
measured device."""

import dataclasses
import math

import numpy as np
import pytest
from scipy.special import expit, psi
from scipy.special import polygamma as real_polygamma

from phantom.charge import derive_double_dot
from phantom.charge.derived import sensor_coupling
from phantom.compile import CompileParams, compile_device
from phantom.constants import E_CHARGE_C, HBAR_EV_S
from phantom.contract import REFERENCE_OPERATING_POINT_V, reference_device
from phantom.device import Device, switcher_scale
from phantom.electrostatics.dots import DerivedDot, barrier_action
from phantom.layout import Layout, two_dot_layout
from phantom.readout.broadening import (
    lineshape_per_ev,
    lineshape_slope_per_ev2,
    occupation,
    peak_conductance_quanta,
    polygamma,
    wing_ev,
)
from phantom.readout.sensor import SensorSETParams

KT = 12.926e-6  # k_B · 150 mK, eV
H_EV_S = 2 * math.pi * HBAR_EV_S


@pytest.fixture(scope="module")
def device() -> Device:
    return reference_device()


# --- the lineshape --------------------------------------------------------------------------------


@pytest.mark.parametrize("m", [1, 2, 3])
def test_polygamma_agrees_with_scipy_on_the_real_axis_and_with_its_own_derivative_off_it(m: int):
    x = np.array([0.5, 0.71, 1.3, 5.0, 33.0])
    np.testing.assert_allclose(polygamma(m, x).real, real_polygamma(m, x), rtol=1e-13)
    z = 0.8 + 1j * np.linspace(-40.0, 40.0, 9)
    h = 1e-5
    below = psi if m == 1 else (lambda w: polygamma(m - 1, w))
    np.testing.assert_allclose((below(z + h) - below(z - h)) / (2 * h), polygamma(m, z), rtol=1e-7)


def test_occupation_is_the_fermi_function_without_broadening_and_the_lorentzian_step_without_temperature():
    eps = np.linspace(-1e-3, 1e-3, 41)
    np.testing.assert_allclose(occupation(eps, 0.0, KT), 1 / (1 + np.exp(eps / KT)), atol=1e-14)
    gamma = 2e-4
    np.testing.assert_allclose(occupation(eps, gamma, 1e-10), 0.5 - np.arctan(2 * eps / gamma) / np.pi, atol=1e-9)


@pytest.mark.parametrize("gamma_ev", [0.0, 1e-5, 2e-4])
def test_the_lineshape_is_normalised_and_peaks_at_the_level(gamma_ev: float):
    e = np.linspace(-0.05, 0.05, 400001)
    lineshape = lineshape_per_ev(e, gamma_ev, KT)
    tails = gamma_ev / (np.pi * 0.05)  # the Lorentzian's weight beyond ±50 meV
    assert np.trapezoid(lineshape, e) == pytest.approx(1.0 - tails, abs=1e-6)
    assert abs(e[np.argmax(lineshape)]) < 0.5 * (e[1] - e[0])
    if gamma_ev == 0.0:
        assert lineshape_per_ev(0.0, 0.0, KT) == pytest.approx(1 / (4 * KT), rel=1e-13)


def test_the_wing_reduces_to_the_inflections_of_sech2_and_of_the_lorentzian_and_is_the_steepest_point():
    assert wing_ev(0.0, KT) == pytest.approx(2 * KT * math.atanh(1 / math.sqrt(3)), rel=1e-12)
    gamma = 1e-3
    assert wing_ev(gamma, gamma * 1e-6) == pytest.approx(gamma / (2 * math.sqrt(3)), rel=1e-5)
    for gamma in (3e-5, 3e-4):
        e = np.linspace(0.0, 10 * (KT + gamma), 400001)
        k = int(np.argmin(lineshape_slope_per_ev2(e, gamma, KT)))
        assert wing_ev(gamma, KT) == pytest.approx(e[k], abs=2 * (e[1] - e[0]))


def test_peak_conductance_is_sequential_tunnelling_when_closed_and_bounded_by_the_quantum_when_open():
    g_small = peak_conductance_quanta(1e5, 3e5, KT)  # γ/kT ≈ 2e-5: the first-order correction is 3.4·γ/(4πkT)
    assert g_small == pytest.approx(H_EV_S * (1e5 * 3e5 / 4e5) / (4 * KT), rel=1e-5), "Γ_eff·h/(4kT)"
    assert peak_conductance_quanta(1e14, 1e14, KT) == pytest.approx(1.0, rel=1e-3), "e²/h for symmetric barriers"
    assert peak_conductance_quanta(1e14, 4e14, KT) == pytest.approx(4 * 4 / 25, rel=1e-3), "4Γ_SΓ_D/(Γ_S + Γ_D)²"
    rates = np.logspace(6, 15, 400)
    assert np.all(peak_conductance_quanta(rates, rates, KT) <= 1.0 + 1e-12)
    assert np.all(peak_conductance_quanta(rates, 0.1 * rates, KT) <= 4 * 0.1 / 1.21 + 1e-12)


def test_the_sensor_current_through_its_leads_is_sequential_tunnelling_in_the_thermal_limit(device: Device):
    """With rates far below kT/ħ the current is e·Γ_SΓ_D/(Γ_S + Γ_D) times the Fermi window."""
    sensor = device.runtime_sensor()
    closed = dataclasses.replace(sensor, leads=dataclasses.replace(sensor.leads, gamma_ref_hz=(1e4, 3e4)))
    phi = sensor.phase_v + np.linspace(-0.3, 0.3, 61) * sensor.period_v
    mu, half, kt = sensor.levels_ev(phi), sensor.bias_v / 2, sensor.kt_ev
    fermi = (expit(-(mu - half) / kt) - expit(-(mu + half) / kt)).sum(axis=-1)
    expected = E_CHARGE_C * (1e4 * 3e4 / 4e4) * fermi
    # far from the peaks the Lorentzian tails (elastic cotunnelling) outlast the thermal exponential
    np.testing.assert_allclose(closed.current_a(phi), expected, rtol=1e-5, atol=1e-6 * expected.max())


# --- the compiled leads ---------------------------------------------------------------------------


def _compile_tier_rates_hz(device: Device, v: np.ndarray) -> np.ndarray:
    """Each lead's rate recomputed from the WKB action at the leads' fixed tunnelling level."""
    leads = device.sensor_leads
    lay, basis = device.basis()
    s = device.sensor_dot
    dot = DerivedDot(s.name, s.x_nm, s.y_nm, s.energy_ev, s.curvature_ev_per_nm2, s.radius_nm, np.array(s.lever_arms))
    points = [lay.gate(n).centroid_nm for n in leads.reservoir_gates]
    actions = [barrier_action(basis, v, dot, p, 0.19, level_ev=lv)[2] for p, lv in zip(points, leads.tunnelling_level_ev, strict=True)]
    return leads.attempt_rate_hz * np.exp(-2.0 * np.array(actions))


def test_the_runtime_rates_are_the_compiled_ones_at_the_anchor_and_follow_the_compile_tier_nearby(device: Device):
    leads = device.sensor_leads
    v0 = np.array(leads.v_ref_v)
    np.testing.assert_allclose(leads.barrier_rates_hz(v0), leads.gamma_ref_hz, rtol=1e-12)
    np.testing.assert_allclose(_compile_tier_rates_hz(device, v0), leads.gamma_ref_hz, rtol=1e-12)
    for gate in ("SB_L", "SB_R", "SG", "SR_L", "SR_R"):
        for dv in (-0.01, 0.01):
            v = v0.copy()
            v[device.gate_names.index(gate)] += dv
            # a few e-folds over ±10 mV; the linearised exponent is within 5 % of the action's change
            np.testing.assert_allclose(leads.barrier_rates_hz(v), _compile_tier_rates_hz(device, v), rtol=0.2, err_msg=f"{gate} {dv:+}")


def test_each_sensor_barrier_controls_its_own_lead_and_a_rate_saturates_at_the_attempt_rate(device: Device):
    leads = device.sensor_leads
    slopes = np.array(leads.log_rate_slopes_per_v)
    jl, jr = device.gate_names.index("SB_L"), device.gate_names.index("SB_R")
    assert slopes[0, jl] > slopes[1, jl] > 0 and slopes[1, jr] > slopes[0, jr] > 0
    v0 = np.array(leads.v_ref_v)
    shut = v0.copy()
    shut[jl] -= 0.05
    ratio = leads.barrier_rates_hz(shut) / leads.barrier_rates_hz(v0)
    assert ratio[0] == pytest.approx(np.exp(-0.05 * slopes[0, jl]), rel=1e-2), "exponential while closed"
    wide = v0.copy()
    wide[[jl, jr]] += 1.0
    np.testing.assert_allclose(leads.barrier_rates_hz(wide), leads.attempt_rate_hz, rtol=1e-6)


def test_a_bath_conducts_above_its_turn_on_and_a_depleted_bath_or_closed_barrier_stops_the_sensor(device: Device):
    sensed = device.runtime_sensed(noise=False)
    v = device.operating_point()
    i_op = float(sensed.current_a(v))
    assert i_op > 1e-10, "the reference device's sensor conducts at its set point"
    np.testing.assert_allclose(device.sensor_leads.bath_occupation(v), 1.0, atol=1e-9)
    for gate, value in (("SR_L", 0.0), ("SR_R", 0.0)):
        w = v.copy()
        w[device.gate_names.index(gate)] = value
        assert float(sensed.current_a(w)) < 1e-6 * i_op, f"{gate} depleted"
    for gate in ("SB_L", "SB_R"):
        w = v.copy()
        w[device.gate_names.index(gate)] -= 0.1
        assert float(sensed.current_a(w)) < 1e-6 * i_op, f"{gate} closed"
    sr = np.linspace(0.0, 0.3, 301)
    w = np.broadcast_to(v, sr.shape + v.shape).copy()
    w[:, device.gate_names.index("SR_L")] = sr
    occ = device.sensor_leads.bath_occupation(w)[:, 0]
    assert np.all(np.diff(occ) >= 0) and occ[0] < 0.01 and occ[-1] > 0.99, "a monotone turn-on inside 0.3 V"


def test_the_current_never_exceeds_the_conductance_quantum(device: Device):
    sensor = device.runtime_sensor()
    v = device.operating_point()
    for opening in (0.0, 0.05, 0.1, 0.3):
        w = v.copy()
        w[[device.gate_names.index("SB_L"), device.gate_names.index("SB_R")]] += opening
        g = sensor.linear_response_peak_a(1e-7, w) / 1e-7  # peak conductance, A/V
        assert g <= E_CHARGE_C / H_EV_S * (1 + 1e-9)


# --- the sensor dot -------------------------------------------------------------------------------


def test_the_textbook_sensor_dot_forms_between_its_barriers_with_baths_beyond():
    lay = two_dot_layout()
    roles = {g.name: g.role for g in lay.gates}
    assert roles["SB_L"] == roles["SB_R"] == "barrier" and roles["SR_L"] == roles["SR_R"] == "reservoir"
    assert lay.gate_names[:8] == ("B_L", "P1", "B_M", "P2", "B_R", "SG", "R_L", "R_R"), "the first eight keep their order"
    xs = {n: lay.gate(n).centroid_nm[0] for n in ("SR_L", "SB_L", "SG", "SB_R", "SR_R")}
    assert xs["SR_L"] < xs["SB_L"] < xs["SG"] < xs["SB_R"] < xs["SR_R"]
    op = dict(REFERENCE_OPERATING_POINT_V, SB_L=0.35, SB_R=0.35)
    d = derive_double_dot(lay, np.array([op[g] for g in lay.gate_names]))
    s = d.sensor_dot
    assert s.x_nm == 0.0 and lay.sensor_nm[2] <= s.y_nm <= lay.sensor_nm[3]
    assert d.sensor_beta_v_per_e[0] == pytest.approx(d.sensor_beta_v_per_e[1], rel=1e-9), "symmetric coupling"


def test_the_sensor_charging_energy_and_coupling_come_from_the_three_dot_capacitance_matrix(device: Device):
    lay, basis = device.basis()
    d = derive_double_dot(lay, device.operating_point(), basis=basis)
    e_c, beta = sensor_coupling(d.ci.c_dot_af, d.dots, d.sensor_dot, lay.stack.epsilon_well)
    assert e_c == d.sensor_charging_energy_ev and np.array_equal(beta, d.sensor_beta_v_per_e)
    assert 1e-3 < e_c < 30e-3, "meV-scale charging energy"
    assert np.all(beta < 0) and np.all(np.abs(beta) < 0.2 * e_c), "an electron raises the sensor level by a fraction of a period"
    dist = [math.hypot(x.x_nm - d.sensor_dot.x_nm, x.y_nm - d.sensor_dot.y_nm) for x in d.dots]
    assert (abs(beta[0]) > abs(beta[1])) == (dist[0] < dist[1]), "the nearer dot couples more strongly"
    sensor = device.runtime_sensor()
    assert sensor.lever_arm_ev_per_v == pytest.approx(1.0, rel=1e-12) and sensor.kt_ev == pytest.approx(KT, rel=1e-3)
    j = device.gate_names.index("SG")
    assert 10e-3 < device.truth().sensor.period_v < 100e-3, "peaks tens of mV apart on the sensor plunger"
    assert device.truth().sensor.period_v == pytest.approx(e_c / sensor.alpha_v_per_v[j], rel=0.05)


def test_the_noise_draws_switchers_keep_their_size_against_a_dot_electron(device: Device):
    amps = device.scenario["noise"]["switcher_amplitudes_v"]
    scale = device.switcher_scale()
    assert scale == pytest.approx(np.mean(np.abs(device.sensor_beta_v_per_e)) / np.mean(np.abs(SensorSETParams().beta_v_per_e)), rel=1e-12)
    session = device.session()
    got = [p.amplitude for p in session._sensed.sensor_noise.processes] if amps else []
    np.testing.assert_allclose(got, np.array(amps) * scale, rtol=1e-12)
    assert switcher_scale(None) == 1.0


def test_a_layout_without_sensor_reservoirs_keeps_the_typed_sensor():
    full = two_dot_layout()
    lay = Layout(tuple(g for g in full.gates if g.name not in ("SB_L", "SB_R", "SR_L", "SR_R")), full.stack, full.material, full.channel_nm, (-80.0, 80.0, 174.0, 246.0))
    op = {g: REFERENCE_OPERATING_POINT_V[g] for g in lay.gate_names}
    dev = compile_device(lay, CompileParams(operating_point_v=op), seed=2)
    assert dev.sensor_leads is None and dev.knowable().sensor_barrier_gates == ()
    sensor = dev.runtime_sensor()
    assert sensor.leads is None and sensor.period_v == SensorSETParams().period_v
    np.testing.assert_array_equal(sensor.beta_v_per_e, SensorSETParams().beta_v_per_e * np.array(dev.sensor_beta_factors()))
    assert dev.at(np.array(dev.programmed_point_v)).sensor_gamma_hz is None
    s = dev.truth().sensor
    assert s.lifetime_width_v == 0.0 and s.gamma_source_hz is None
    assert s.wing_v - s.peak_v == pytest.approx(math.atanh(1 / math.sqrt(3)) * s.width_v, rel=1e-9)
