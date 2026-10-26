"""P0 for the wiring faults and the trapping hysteresis: documented effect, clean when severity is zero,
seeded and serialisable."""

import numpy as np
import pytest

from phantom.disorder import (
    ScenarioDraw,
    ScenarioParams,
    TrappingParams,
    WiringDraw,
    WiringFaultParams,
    draw_scenario,
    draw_wiring_faults,
    effective_gate_voltage,
    trapped_fraction,
)
from phantom.disorder.scenario import adjacent_gate_pairs
from phantom.layout import two_dot_layout

NAMES = ("B_L", "P1", "B_M", "P2", "B_R", "SG", "R_L", "R_R")


def test_clean_wiring_is_the_identity():
    w = WiringDraw(NAMES)
    v = np.random.default_rng(0).uniform(-1, 1, (7, 8))
    assert w.is_clean
    np.testing.assert_array_equal(w.apply(v), v)
    np.testing.assert_array_equal(w.leakage_current_a(v), 0.0)


def test_shorted_follower_takes_the_kept_voltage_and_open_gate_floats():
    w = WiringDraw(NAMES, open_gates=("B_L",), shorted_pairs=(("P2", "B_M"),), stuck_potential_v=0.1)
    v = np.array([0.2, 0.6, 0.2, 0.7, 0.2, 0.6, 0.7, 0.7])
    a = w.apply(v)
    assert a[0] == 0.1, "open gate floats at the stuck potential whatever is programmed"
    assert a[2] == a[3] == 0.7, "B_M follows P2"
    assert a[1] == 0.6 and a[4] == 0.2 and a[5] == 0.6
    m, c = w.wiring_matrix()
    assert m.shape == (8, 8) and np.all(m[0] == 0) and c[0] == 0.1


def test_leakage_pair_currents_are_equal_and_opposite_and_linear():
    w = WiringDraw(NAMES, leaks=(("P1", "B_M", 1e9),))
    v = np.zeros((5, 8))
    v[:, 1] = np.linspace(-1, 1, 5)
    i = w.leakage_current_a(v)
    np.testing.assert_allclose(i[:, 1], v[:, 1] / 1e9)
    np.testing.assert_allclose(i[:, 2], -v[:, 1] / 1e9)
    assert np.all(i[:, [0, 3, 4, 5, 6, 7]] == 0)


def test_wiring_draw_is_clean_at_zero_severity_and_seeded_otherwise():
    lay = two_dot_layout()
    pairs = adjacent_gate_pairs(lay)
    assert pairs == (("B_L", "P1"), ("B_L", "R_L"), ("P1", "B_M"), ("B_M", "P2"), ("P2", "B_R"), ("B_R", "R_R"))
    assert draw_wiring_faults(WiringFaultParams(), lay.gate_names, pairs, np.random.default_rng(0), severity=0.0).is_clean
    a = draw_wiring_faults(WiringFaultParams(), lay.gate_names, pairs, np.random.default_rng(4), severity=3.0)
    b = draw_wiring_faults(WiringFaultParams(), lay.gate_names, pairs, np.random.default_rng(4), severity=3.0)
    assert a.to_dict() == b.to_dict()
    for kept, follower in a.shorted_pairs:
        assert kept not in a.open_gates and follower not in a.open_gates
    # over many draws the open-gate rate is the declared probability
    n_open = sum(len(draw_wiring_faults(WiringFaultParams(p_open_per_gate=0.1), lay.gate_names, pairs, np.random.default_rng(s)).open_gates) for s in range(500))
    assert n_open / (500 * 8) == pytest.approx(0.1, abs=0.02)


def test_scenario_includes_wiring_and_round_trips():
    lay = two_dot_layout()
    d = draw_scenario(lay, ScenarioParams.for_material("SiGe"), seed=2, severity=3.0)
    assert d.wiring is not None
    back = ScenarioDraw.from_dict(d.to_dict())
    assert back.wiring.to_dict() == d.wiring.to_dict()
    clean = draw_scenario(lay, ScenarioParams.for_material("SiGe"), seed=2, severity=0.0)
    assert clean.wiring.is_clean
    off = draw_scenario(lay, ScenarioParams.for_material("SiGe"), seed=2, severity=3.0, wiring_faults=False)
    assert off.wiring.is_clean and off.charges.count == d.charges.count


def test_trapping_opens_a_loop_and_is_clean_at_zero_severity():
    up = np.linspace(-0.6, 1.2, 361)
    seq = np.concatenate([up, up[::-1]])
    params = TrappingParams()
    v_eff, f = effective_gate_voltage(seq, params)
    n = len(up)
    assert f.max() == pytest.approx(1.0, abs=1e-4) and f[0] == 0.0
    assert np.all(np.diff(f[:n]) >= 0), "charging only on the way up"
    hold = (up > params.release_v) & (up < params.capture_v)
    assert np.all(v_eff[n:][::-1][hold] < v_eff[:n][hold]), "on the way down the gate is screened by the trapped charge"
    assert np.all(v_eff[:n][up <= params.capture_v] == up[up <= params.capture_v]), "nothing trapped before capture on the way up"
    v0, f0 = effective_gate_voltage(seq, params, severity=0.0)
    np.testing.assert_array_equal(v0, seq)
    assert np.all(f0 == 0)
    with pytest.raises(ValueError):
        TrappingParams(capture_v=0.0, release_v=0.5)
    assert trapped_fraction(np.array([1.0, 1.0, 1.0]), TrappingParams(capture_rate_per_sample=1.0))[0] == 1.0
