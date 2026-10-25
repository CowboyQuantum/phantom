"""P0: the SET sensor's Coulomb peaks and diamonds, and the charge-transition step it reads."""

import numpy as np
import pytest

from phantom.charge import ConstantInteractionParams
from phantom.readout import SensedDoubleDot, SensorSETParams


def test_sensor_is_periodic_and_peaked():
    s = SensorSETParams()
    phi = np.linspace(-0.6, 0.6, 12001)
    i = s.current_a(phi)
    assert i.max() == pytest.approx(s.peak_a + s.background_a, rel=0.03)
    valley = s.background_a + 2 * s.peak_a / np.cosh(s.period_v / 2 / s.width_v) ** 2
    assert i.min() == pytest.approx(valley, rel=0.05)
    tiny = (s.current_a(phi, v_sd=1e-9) - s.background_a) * s.bias_v / 1e-9
    k = np.arange(-3, 4)
    d = phi[:, None] - s.phase_v - k * s.period_v
    sech = (s.peak_a / np.cosh(d / s.width_v) ** 2).sum(axis=1)
    np.testing.assert_allclose(tiny, sech, rtol=1e-4)
    peaks = np.flatnonzero((i[1:-1] > i[:-2]) & (i[1:-1] > i[2:])) + 1
    np.testing.assert_allclose(np.diff(phi[peaks]), s.period_v, rtol=1e-2)
    np.testing.assert_allclose(phi[peaks], s.peak_positions_v(-0.6, 0.6), atol=2e-4)


def test_sensor_bias_opens_coulomb_diamonds_with_slopes_from_the_lever_arm():
    s = SensorSETParams()
    phi = s.phase_v + s.period_v / 2
    e_c = s.charging_energy_ev
    assert s.current_a(phi, 0.2 * e_c) < s.background_a + 0.1 * s.peak_a
    assert s.current_a(phi, 1.3 * e_c) > s.background_a + s.peak_a
    assert s.current_a(phi, -1.3 * e_c) < s.background_a - s.peak_a
    (phi_0, s_up), (_, s_down) = s.diamond_edges(0)
    assert phi_0 == s.phase_v and s_down == -s_up
    assert s_up == pytest.approx(2 * s.lever_arm_ev_per_v)
    assert (s.period_v / 2) * s_up == pytest.approx(e_c, rel=1e-12), "apex at V_sd = E_c / e"
    for d in (-0.08, 0.05):
        v = s_up * d
        i_on, i_off = s.current_a(phi_0 + d, v * 1.3), s.current_a(phi_0 + d, v * 0.7)
        assert min(i_on, i_off) < s.current_a(phi_0 + d, v) < max(i_on, i_off)


def test_charge_transition_steps_the_sensor():
    d = SensedDoubleDot(charge=ConstantInteractionParams(temperature_k=0.01))
    v = np.zeros((3000, 3))
    v[:, 0] = np.linspace(0.0, 0.02, 3000)
    n1 = d.charge.mean_occupation(v)[:, 0]
    k = int(np.argmax(n1 > 0.5))
    i = d.current_a(v)
    assert abs(i[k + 2] - i[k - 2]) > 20 * abs(i[k + 200] - i[k + 196])


def test_parameters_are_validated():
    with pytest.raises(ValueError):
        SensorSETParams(width_v=0.0)
    with pytest.raises(ValueError):
        SensedDoubleDot(sensor=SensorSETParams(alpha_v_per_v=np.array([0.01, 0.01])))
