"""P0: the layout-derived transport shapes behave as their definitions say."""

import numpy as np
import pytest

from phantom.layout import two_dot_layout
from phantom.readout import channel_transmission_current_a, gate_leakage_current_a, ohmic_current_a


def test_channel_transmission_is_a_monotonic_step_between_zero_and_saturation():
    u = np.linspace(-0.2, 0.1, 601)
    i = channel_transmission_current_a(u, threshold_ev=-0.05, width_ev=0.004, i_sat_a=2e-9)
    assert np.all(np.diff(i) <= 0), "the channel closes as the barrier energy rises"
    assert i[0] == pytest.approx(2e-9, rel=1e-9) and i[-1] < 1e-20
    assert channel_transmission_current_a(np.array([-0.05]), -0.05, 0.004, 2e-9)[0] == pytest.approx(1e-9)


def test_gate_leakage_is_antisymmetric_scales_with_area_and_is_tiny_in_the_window():
    lay = two_dot_layout()
    t = lay.stack.depth_nm
    v = np.linspace(-3.0, 3.0, 61)
    big = gate_leakage_current_a(v, lay.gate("P1").area_nm2, t)
    small = gate_leakage_current_a(v, lay.gate("B_L").area_nm2, t)
    np.testing.assert_allclose(big, -big[::-1], atol=1e-30)
    np.testing.assert_allclose(big, small * lay.gate("P1").area_nm2 / lay.gate("B_L").area_nm2)
    inside = np.abs(v) < 1.0
    assert np.abs(big[inside]).max() < 1e-3 * np.abs(big).max(), "negligible inside the safe window"
    v_onset = 0.08 * t
    at_onset = gate_leakage_current_a(np.array([v_onset]), 1e6, t)[0]  # one square micrometre
    assert at_onset == pytest.approx(1e-9, rel=1e-9)


def test_ohmic_current_is_linear():
    v = np.linspace(-1e-3, 1e-3, 11)
    np.testing.assert_allclose(ohmic_current_a(v, 20e3), v / 20e3)
