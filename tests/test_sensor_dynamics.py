"""P0 for the sensor dot's charge in time (#54, ADR 0028): the counting statistics against the two-state
closed forms and exact trajectories, the long-dwell mean against the stationary current of ADR 0026, the
shot noise shrinking with the dwell, a slow sensor switching as a two-state process, no current through a
depleted bath or a closed barrier, nothing changed without dynamics, and the back-action on the dots.
Analytic limits and invariants only."""

import math

import numpy as np
import pytest

from phantom.constants import E_CHARGE_C
from phantom.contract import reference_device
from phantom.device import Device
from phantom.readout import sensor_dynamics as sd


@pytest.fixture(scope="module")
def device() -> Device:
    return reference_device()


def _point(device: Device, **shift: float) -> np.ndarray:
    """The programmed point with the sensor gate at its compiled set point, shifted per gate."""
    v = np.array(device.programmed_point_v)
    v[device.gate_names.index("SG")] = device.operating_point_v[device.gate_names.index("SG")]
    for gate, dv in shift.items():
        v[device.gate_names.index(gate)] += dv
    return v


def _on_peak(device: Device, v: np.ndarray, offset: float = 0.0) -> np.ndarray:
    """`v` with the sensor gate moved so the sensor's nearest level sits `offset` kT above the leads,
    with the dots in their ground state (the reference device's SG is not faulted)."""
    sensor, j = device.runtime_sensor(), device.gate_names.index("SG")
    out = v.copy()
    for _ in range(2):
        s = device.at(out)
        x = (float(sensor.potential(s.v_carried_v, s.n_ground)) - sensor.phase_v) / sensor.period_v
        out[j] += (np.rint(x) - offset * sensor.kt_ev / sensor.charging_energy_ev - x) * sensor.period_v / sensor.alpha_v_per_v[j]
    return out


def _chain_at(device: Device, v: np.ndarray) -> tuple[sd.SensorChain, float, np.ndarray]:
    s = device.at(v)
    sensor = device.runtime_sensor()
    phi = float(sensor.potential(s.v_carried_v, s.n_ground))
    lo, hi = sd.window_states(sensor, phi, sensor.bias_v)
    return sd.chain(sensor, phi, s.v_carried_v, sensor.bias_v, lo, hi), phi, s.v_carried_v


SLOW = {"SB_L": -0.085, "SB_R": -0.06}
"""Both sensor barriers closed so that each rate is in the 10³–10⁵ Hz range on the reference device."""


def test_the_counting_statistics_are_the_two_state_closed_forms(device: Device):
    c, _, _ = _chain_at(device, _on_peak(device, _point(device), 1.0))
    assert len(c.states) == 2
    a_s, a_d = c.w[1, 0] - (-c.j1[1, 0]), -c.j1[1, 0]  # in from the source, in from the drain
    b_d = c.j1[0, 1]
    b_s = c.w[0, 1] - b_d
    total = a_s + a_d + b_s + b_d
    mean, var = sd.counting_cumulants(c)
    assert mean == pytest.approx((a_s * b_d - a_d * b_s) / total, rel=1e-9)
    assert var == pytest.approx((a_s * b_d + a_d * b_s) / total - 2 * (a_s * b_d - a_d * b_s) ** 2 / total**3, rel=1e-9)
    s = device.at(_on_peak(device, _point(device), 1.0))
    assert E_CHARGE_C * mean + c.background_a == pytest.approx(float(s.sensor_current_a), rel=1e-9), "the stationary current is ADR 0026's"
    np.testing.assert_allclose(c.w @ sd.stationary(c), 0.0, atol=1e-12 * total)


def test_exact_trajectories_reproduce_the_counting_statistics(device: Device):
    v = _on_peak(device, _point(device, **SLOW), 0.5)
    c, phi, vc = _chain_at(device, v)
    sensor = device.runtime_sensor()
    total = float(np.max(-np.diag(c.w)))
    assert 1e3 < total < 1e6
    dwell = 50.0 / total
    rng = np.random.default_rng(5)
    p = sd.stationary(c)
    counts = np.empty(4000)
    for k in range(len(counts)):
        m0 = int(c.states[int(rng.random() > p[0])])
        counts[k] = sd.trajectory(sensor, phi, vc, sensor.bias_v, m0, dwell, rng)[0]
    mean, var = sd.counting_cumulants(c)
    assert counts.mean() / dwell == pytest.approx(mean, rel=0.05)
    assert counts.var() / dwell == pytest.approx(var, rel=0.12), "the long-time variance, within its 1/(ΓT) correction"


def test_with_fast_rates_the_long_dwell_mean_is_the_stationary_current_and_the_noise_shrinks(device: Device):
    v = _on_peak(device, _point(device), 1.0)
    expected = float(device.at(v).sensor_current_a)
    run = device.session(noise=False, charge="mean").sense(v, 1e-3, samples=3)
    # the integral of e^{Ws} by one matrix exponential at ΓT ≈ 10⁸ keeps about eight digits
    np.testing.assert_allclose(run.current_a, expected, rtol=1e-6)
    c, _, _ = _chain_at(device, v)
    _, var = sd.counting_cumulants(c)
    spreads = []
    for dwell in (1e-6, 1e-4):
        i = device.session(noise=False).sense(v, dwell, samples=2000).current_a
        sigma = E_CHARGE_C * math.sqrt(var / dwell)
        assert abs(i.mean() - expected) < 5 * sigma / math.sqrt(len(i))
        assert np.std(i) == pytest.approx(sigma, rel=0.08)
        spreads.append(np.std(i))
    assert spreads[0] / spreads[1] == pytest.approx(10.0, rel=0.12), "σ ∝ 1/√dwell"


def test_a_slow_sensor_switches_as_the_two_state_process_and_its_current_comes_in_electrons(device: Device):
    v = _on_peak(device, _point(device, **SLOW), 0.5)
    c, _, _ = _chain_at(device, v)
    up, down = c.w[1, 0], c.w[0, 1]
    p = up / (up + down)
    assert 0.1 < p < 0.9
    for rate_dwell in (0.05, 0.5):
        dwell = rate_dwell / (up + down)
        m = device.session(noise=False).sense(v, dwell, samples=6000)
        charge = m.sensor_charge
        assert set(np.unique(charge)) <= set(c.states.astype(float).tolist())
        changes = float(np.mean(charge[1:] != charge[:-1]))
        expected = 2 * p * (1 - p) * (1 - math.exp(-rate_dwell))
        assert abs(changes - expected) < 6 * math.sqrt(expected / 6000) + 0.002
        electrons = (m.current_a - c.background_a) * dwell / E_CHARGE_C
        np.testing.assert_allclose(electrons, np.rint(electrons), atol=1e-6)
        assert np.all(np.abs(np.rint(electrons)) <= 3), "a few electrons per dwell: bursts"


@pytest.mark.parametrize("charge", ["sampled", "mean"])
def test_a_depleted_bath_or_a_closed_barrier_carries_no_current(device: Device, charge: str):
    v = _on_peak(device, _point(device), 1.0)
    i_open = float(device.at(v).sensor_current_a)
    for shift in ({"SR_L": -0.7}, {"SR_R": -0.7}, {"SB_L": -0.3}, {"SB_R": -0.3}):
        w = v.copy()
        for gate, dv in shift.items():
            w[device.gate_names.index(gate)] += dv
        i = device.session(noise=False, charge=charge).sense(w, 1e-3, samples=20).current_a
        assert np.all(np.abs(i) < 1e-6 * i_open), shift


def test_without_dynamics_nothing_changes(device: Device):
    v = np.stack([_point(device), _on_peak(device, _point(device), 1.0)])
    for charge in ("sampled", "mean"):
        m = device.session(noise=False, dynamics=False, charge=charge).sense(v, 1e-3)
        assert m.sensor_charge is None
        sensor = device.runtime_sensor()
        vc = device.effective_voltages(v)
        np.testing.assert_array_equal(m.current_a, sensor.current_a(sensor.potential(vc, m.n_present), sensor.bias_v, vc))


def test_the_back_action_is_the_sensor_charge_against_its_mean(device: Device):
    """In `mean` mode a sensor in its stationary state does not act back; one electron above its mean
    shifts dot i's levels by −β_i (the coupling by which dot i shifts the sensor), and the dots' rates
    keep detailed balance in the shifted energies."""
    v = _on_peak(device, _point(device), 1.0)
    vc = device.effective_voltages(v)
    run = device.session(noise=False, charge="mean")
    run.hold(v, 1e-3)
    assert np.all(np.abs(run._back_action_ev(vc, 0.0)) < 1e-12)
    run._m_p = run._m_p.copy()
    m_bar = float(run._m_p @ run._m_states)
    run._m_states, run._m_p = np.array([round(m_bar) + 1]), np.ones(1)
    shift = run._back_action_ev(vc, 0.0)
    beta = device.runtime_sensor().beta_v_per_e
    np.testing.assert_allclose(shift, -beta * (round(m_bar) + 1 - m_bar), rtol=1e-9)
    assert np.all(shift > 0), "a sensor electron raises the dots' levels"
    rates = device.runtime_rates()
    states = rates.window_states(vc)
    u, w = rates.rate_matrix(vc, states, shift)
    boltzmann = np.exp(-(u - u.min()) / rates.kt_ev)
    np.testing.assert_allclose(u, device.runtime_ci().state_energies_ev(vc, states) + states @ shift, rtol=1e-12)
    assert np.abs(w @ boltzmann).max() <= 1e-9 * np.abs(w).max()


def test_a_typed_sensor_has_no_charge_process():
    from phantom.compile import CompileParams, compile_device
    from phantom.contract import REFERENCE_OPERATING_POINT_V
    from phantom.layout import Layout, two_dot_layout

    full = two_dot_layout()
    lay = Layout(tuple(g for g in full.gates if g.name not in ("SB_L", "SB_R", "SR_L", "SR_R")), full.stack, full.material, full.channel_nm, (-80.0, 80.0, 174.0, 246.0))
    dev = compile_device(lay, CompileParams(operating_point_v={g: REFERENCE_OPERATING_POINT_V[g] for g in lay.gate_names}), seed=2)
    assert dev.session().sense(np.array(dev.programmed_point_v), 1e-3, samples=3).sensor_charge is None
