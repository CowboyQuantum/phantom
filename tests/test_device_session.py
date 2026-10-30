"""P0 for Step 8: the runtime session and the truth record against closed forms. A sampled charge is a
Boltzmann draw without dynamics and a Markov trajectory with them, whose switching matches the
two-state propagator; white readout noise keeps its density as the dwell changes; a noise seed keeps the
noise parameters; the sensor truth is the peak, width, wing and slope of the sensor's own current, thermal
and lifetime-broadened; the charge lattice follows the wiring draw; the spin readout reproduces the PSB
occupations."""

import dataclasses
import math

import numpy as np
import pytest

from phantom.api import Device
from phantom.constants import K_B_EV_PER_K
from phantom.contract import reference_device
from phantom.spin import psb_occupation


@pytest.fixture(scope="module")
def device() -> Device:
    return reference_device()


def _dot1_degeneracy(dev: Device) -> np.ndarray:
    """A programmed point where dot 1 sits half-way through a transition to its lead."""
    v0 = np.array(dev.programmed_point_v)
    j = dev.gate_names.index("P1")
    grid = v0 + np.linspace(-0.03, 0.03, 4001)[:, None] * np.eye(len(v0))[j]
    n = dev.at(grid).n_mean
    frac = n[:, 0] - np.floor(n[:, 0])
    flat2 = np.abs(n[:, 1] - np.rint(n[:, 1])) < 1e-6
    k = int(np.argmin(np.where(flat2, np.abs(frac - 0.5), np.inf)))
    assert abs(frac[k] - 0.5) < 0.05
    return grid[k]


def test_sampled_charge_without_dynamics_is_a_boltzmann_draw(device: Device):
    v = _dot1_degeneracy(device)
    p = device.at(v).n_mean[0] - np.floor(device.at(v).n_mean[0])
    n = device.session(noise=False, dynamics=False).sense(v, 1e-3, samples=4000).n_present[:, 0]
    upper = np.mean(n > np.floor(device.at(v).n_mean[0]))
    assert abs(upper - p) < 5 * math.sqrt(p * (1 - p) / 4000)


@pytest.mark.parametrize("gamma_dwell", [0.02, 0.5, 50.0])
def test_the_trajectory_switches_as_the_two_state_propagator(device: Device, gamma_dwell: float):
    """At a dot-lead degeneracy the charge is a two-state Markov process with total rate Γ_L: between
    consecutive samples it changes with probability 2p(1−p)(1 − e^{−Γ_L·dwell}) on average."""
    v = _dot1_degeneracy(device)
    gamma = device.at(v).gamma_left_hz
    dwell = gamma_dwell / gamma
    n = device.session(noise=False).sense(v, dwell, samples=6000).n_present[:, 0]
    p = float(np.mean(n > n.min()))
    changes = float(np.mean(n[1:] != n[:-1]))
    expected = 2 * p * (1 - p) * (1 - math.exp(-gamma_dwell))
    assert set(np.unique(n)) <= {n.min(), n.min() + 1}
    assert abs(changes - expected) < 6 * math.sqrt(expected / 6000) + 0.002


def test_white_readout_noise_keeps_its_density_as_the_dwell_changes(device: Device):
    noise = {**device.scenario["noise"], "one_over_f_amplitude_v": 0.0, "switcher_amplitudes_v": [], "switcher_rates_up_hz": [], "switcher_rates_down_hz": []}
    quiet = dataclasses.replace(device, scenario={**device.scenario, "noise": noise})
    v = np.array(device.programmed_point_v)
    floor = device.at(v).sensor_current_a
    sigma_1ms = noise["white_sigma_a_at_1ms"]
    for dwell in (1e-3, 1e-4):
        i = quiet.session(dynamics=False, charge="mean").sense(v, dwell, samples=4000).current_a
        assert np.std(i - floor) == pytest.approx(sigma_1ms * math.sqrt(1e-3 / dwell), rel=0.06)


def test_a_noise_seed_changes_the_realisation_not_the_parameters(device: Device):
    a, b = device.session(noise_seed=1), device.session(noise_seed=2)
    assert a._sensed.gate_noise.amplitude == b._sensed.gate_noise.amplitude == device.noise_sigma_gate_v()
    amps = [[p.amplitude for p in s._sensed.sensor_noise.processes] for s in (a, b)] if device.scenario["noise"]["switcher_amplitudes_v"] else [[], []]
    assert amps[0] == amps[1]
    assert a._sensed.readout_noise.sigma == b._sensed.readout_noise.sigma
    v = np.array(device.programmed_point_v)
    assert not np.array_equal(a.sense(v, 1e-3, samples=50).current_a, b.sense(v, 1e-3, samples=50).current_a)


@pytest.mark.parametrize("at_set_point", [False, True], ids=["programmed-thermal", "set-point-broadened"])
def test_sensor_truth_is_the_peak_wing_and_slope_of_the_sensor_current(device: Device, at_set_point: bool):
    """The peak nearest the programmed sensor gate is nearly thermal on the reference device; the one at
    the compiled set point is lifetime-broadened (ħΓ about four thermal widths). The truth states each
    peak's lineshape with the leads' rates at its peak, so the current is read with the rates held there."""
    bias = 1e-7  # linear response
    near = device.operating_point_v[device.gate_names.index("SG")] if at_set_point else None
    s = device.truth(near_v=near, v_sd=bias).sensor
    assert (s.lifetime_width_v > 2 * s.width_v) is at_set_point
    j = device.gate_names.index(s.gate)
    sensor = device.runtime_sensor()
    v_ref = np.array(device.programmed_point_v)
    n_ref = device.runtime_ci().mean_occupation(device.effective_voltages(v_ref))
    v_peak = v_ref.copy()
    v_peak[j] = s.peak_v
    sg = s.peak_v + np.linspace(-0.05, 0.05, 200001) * s.period_v
    v = np.broadcast_to(v_ref, sg.shape + v_ref.shape).copy()
    v[:, j] = sg
    i = sensor.current_a(sensor.potential(device.effective_voltages(v), n_ref), bias, device.effective_voltages(v_peak))
    step = sg[1] - sg[0]
    assert abs(sg[np.argmax(i)] - s.peak_v) <= step
    slope = -np.gradient(i, step)
    k = int(np.argmax(np.where(sg > s.peak_v, slope, -np.inf)))
    assert abs(sg[k] - s.wing_v) < 2 * step
    assert slope[k] == pytest.approx(s.sensitivity_a_per_v, rel=2e-3)
    assert i.max() == pytest.approx(sensor.linear_response_peak_a(bias, device.effective_voltages(v_peak)), rel=2e-3)
    assert s.width_v == pytest.approx(sensor.width_v / sensor.alpha_v_per_v[j], rel=1e-12)
    rates = sensor.rates_hz(device.effective_voltages(v_peak))
    assert (s.gamma_source_hz, s.gamma_drain_hz) == (float(rates[0]), float(rates[1]))


def test_sensor_truth_agrees_with_the_optimum_flank(device: Device):
    """The optimum searches the actual current at the default bias, whose rates change along the sensor
    gate (a factor e per few mV, ADR 0026); the truth states the linear-response lineshape with the rates
    at the peak. Across the 0.2 mV from the peak to its wing that is a few per cent of slope and a tenth
    of the wing's offset."""
    from phantom.api import optimum

    o = optimum(device).sensor
    s = device.truth(near_v=o.v_sensor_gate_v).sensor
    assert o.slope_a_per_v == pytest.approx(s.sensitivity_a_per_v, rel=0.08)
    assert abs(abs(o.v_sensor_gate_v - s.peak_v) - (s.wing_v - s.peak_v)) < 0.2 * (s.wing_v - s.peak_v)


def test_charge_lattice_and_sensor_truth_follow_the_wiring_draw(device: Device):
    names = device.gate_names
    wiring = {**device.scenario["wiring"], "shorted_pairs": [["P1", "B_M"]], "open_gates": [names[names.index("SG")]]}
    faulty = dataclasses.replace(device, scenario={**device.scenario, "wiring": wiring})
    c = np.array(device.c_gate_af)
    ip, ib, i2 = names.index("P1"), names.index("B_M"), names.index("P2")
    from phantom.constants import E_UNIT_V_PER_E_AF

    expected = np.stack([c[:, ip] + c[:, ib], c[:, i2]], axis=1) / E_UNIT_V_PER_E_AF
    np.testing.assert_allclose(np.array(faulty.truth().charge_lattice_e_per_v), expected, rtol=1e-12)
    assert faulty.truth().sensor is None, "an open sensor gate has no peak to find"
    clean = dataclasses.replace(device, scenario={**device.scenario, "wiring": None})
    np.testing.assert_allclose(np.array(clean.truth().charge_lattice_e_per_v), clean.runtime_ci().charge_lattice(), rtol=1e-12)


def test_knowable_hides_the_compiled_set_point(device: Device):
    k = device.knowable()
    j = device.gate_names.index(k.sensor_gate)
    assert k.programmed_point_v[j] != device.operating_point_v[j], "the sensor set point is truth"
    assert not {"seed", "scenario", "operating_point_v", "params"} & set(k.to_dict())


@pytest.mark.parametrize("p_triplet", [0.0, 1.0])
def test_spin_readout_reproduces_the_psb_occupations(device: Device, p_triplet: float):
    eps = 0.5 * device.dots[1].valley_splitting_uev * 1e-6
    out = device.session().spin_readout(p_triplet, eps, 4000, 1e-6, 1.0)
    p_s, p_t = psb_occupation(np.array([eps]), device.dots[1].valley_splitting_uev * 1e-6, device.tunnel_coupling_ev, K_B_EV_PER_K * device.temperature_k)
    expected = float(p_t[0] if p_triplet else p_s[0]) * (1 - 1e-6) + float(p_s[0]) * 1e-6 * p_triplet
    assert out.prepared_triplet.all() == bool(p_triplet) and out.prepared_triplet.any() == bool(p_triplet)
    assert abs(out.in_02.mean() - expected) < 5 * math.sqrt(max(expected * (1 - expected), 1e-4) / 4000)


def test_a_relaxed_triplet_reads_as_a_singlet(device: Device):
    e_v = device.dots[1].valley_splitting_uev * 1e-6
    eps = 0.5 * e_v
    out = device.session().spin_readout(1.0, eps, 4000, 1e-3, 1e-3)
    p_s, p_t = psb_occupation(np.array([eps]), e_v, device.tunnel_coupling_ev, K_B_EV_PER_K * device.temperature_k)
    survive = math.exp(-1.0)
    expected = survive * float(p_t[0]) + (1 - survive) * float(p_s[0])
    assert abs(out.in_02.mean() - expected) < 5 * math.sqrt(expected * (1 - expected) / 4000)
