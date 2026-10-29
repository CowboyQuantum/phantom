"""P0 for Step 4's dynamics core: WKB rates, detailed balance, and latching above the rate threshold."""

import numpy as np
import pytest

from phantom.charge import ConstantInteractionParams, derive_double_dot
from phantom.electrostatics import GateBasis, barrier_width_nm
from phantom.layout import MaterialClass, two_dot_layout
from phantom.readout import RateEquation, attempt_rate_hz, wkb_transmission, wkb_tunnel_rate_hz


def test_wkb_transmission_limits_and_monotonicity():
    m = 0.19
    assert wkb_transmission(0.0, 10.0, m) == 1.0 and wkb_transmission(-1.0, 10.0, m) == 1.0
    heights = np.linspace(1e-3, 50e-3, 50)
    t = wkb_transmission(heights, 20.0, m)
    assert np.all(np.diff(t) < 0) and t[-1] < 1e-6
    widths = np.linspace(1.0, 60.0, 60)
    assert np.all(np.diff(wkb_transmission(10e-3, widths, m)) < 0)
    # κ for 10 meV in silicon is 0.22 /nm: a 20 nm barrier transmits exp(-8.9)
    kappa = np.sqrt(10e-3 * m / 0.0380998212)
    assert wkb_transmission(10e-3, 20.0, m) == pytest.approx(np.exp(-2 * kappa * 20.0), rel=1e-12)
    assert 1e9 <= attempt_rate_hz(1e-5, m) <= 1e13


def test_lead_rate_falls_by_decades_as_the_barrier_closes():
    lay = two_dot_layout()
    basis = GateBasis.from_layout(lay)
    m = MaterialClass.sige().effective_mass_me
    rates = []
    for vb in (0.2, 0.1, 0.0):
        v = np.array([vb, 0.6, 0.2, 0.6, vb, 0.6, 0.7, 0.7, 0.4, 0.3, 0.7, 0.7])
        d = derive_double_dot(lay, v, basis=basis)
        h, w = barrier_width_nm(basis, v, d.dots[0], lay.reservoir_points_nm()["R_L"])
        rates.append(wkb_tunnel_rate_hz(h, w, d.dots[0].curvature_ev_per_nm2, m))
    assert rates[0] > rates[1] > rates[2]
    assert rates[0] / rates[2] > 1e3


def test_stationary_distribution_is_boltzmann_detailed_balance():
    ci = ConstantInteractionParams(temperature_k=0.3)
    re = RateEquation(ci, gamma_left_hz=2e4, gamma_right_hz=7e3, gamma_interdot_hz=5e3)
    rng = np.random.default_rng(1)
    for _ in range(10):
        v = np.array([rng.uniform(0, 0.1), rng.uniform(0, 0.1), 0.0])
        states, p = re.stationary(v)
        u_states = ci.state_energies_ev(v, states)
        pb = np.exp(-(u_states - u_states.min()) / re.kt_ev)
        pb /= pb.sum()
        np.testing.assert_allclose(p, pb, atol=1e-9)
        # and the CI model's own thermal occupation agrees with the stationary mean
        np.testing.assert_allclose(p @ states, ci.mean_occupation(v), atol=1e-9)
        _, u, w = re.states_and_rates(v)
        kt = re.kt_ev
        for i in range(len(u)):
            for j in range(len(u)):
                if i != j and w[j, i] > 0:
                    assert w[j, i] / w[i, j] == pytest.approx(np.exp(-(u[j] - u[i]) / kt), rel=1e-9)
        assert np.allclose(w.sum(axis=0), 0.0, atol=1e-6 * np.abs(w).max())


def test_fast_rates_follow_the_sweep_and_slow_rates_latch():
    ci = ConstantInteractionParams(temperature_k=0.05)
    v = np.zeros((300, 3))
    v[:, 0] = np.linspace(0.0, 0.05, 300)  # across three dot-1 transitions
    equilibrium = ci.mean_occupation(v)[:, 0]
    dwell = 1e-3
    # an interdot rate is needed too: with cold leads a (2,1) → (3,0) rearrangement has both lead paths
    # uphill by tens of kT and would stay metastable, which is real physics but not this test's point
    fast = RateEquation(ci, 1e5, 1e5, 1e5).sweep(v, dwell)[0][:, 0]  # Γ·dwell = 100
    np.testing.assert_allclose(fast, equilibrium, atol=0.02)
    slow = RateEquation(ci, 1e-2, 1e-2, 1e-2).sweep(v, dwell)[0][:, 0]  # Γ·dwell = 1e-5
    assert np.all(np.abs(slow - slow[0]) < 0.02), "the charge never follows: latched at its starting value"
    middle = RateEquation(ci, 1e3, 1e3, 1e3).sweep(v, dwell)[0][:, 0]  # Γ·dwell = 1
    assert np.all(middle <= equilibrium + 0.02) and (equilibrium - middle).max() > 0.2, "lags behind"


def test_direction_dependence_of_a_lagging_sweep():
    ci = ConstantInteractionParams(temperature_k=0.05)
    up = np.zeros((200, 3))
    up[:, 0] = np.linspace(0.0, 0.04, 200)
    re = RateEquation(ci, 2e2, 2e2, 2e2)
    n_up = re.sweep(up, 1e-3)[0][:, 0]
    n_down = re.sweep(up[::-1], 1e-3)[0][:, 0][::-1]
    assert np.all(n_down >= n_up - 1e-9) and (n_down - n_up).max() > 0.3, "hysteresis between directions"


def test_without_an_interdot_rate_a_cold_double_dot_can_stay_metastable():
    """Sequential tunnelling only: when a rearrangement (2,1) → (3,0) needs an electron to leave dot 2 or
    enter dot 1 against tens of kT, both lead paths are suppressed and the charge stays where it was,
    although the ground state has moved. An interdot rate releases it."""
    ci = ConstantInteractionParams(temperature_k=0.05)
    v = np.zeros((300, 3))
    v[:, 0] = np.linspace(0.0, 0.05, 300)
    stuck = RateEquation(ci, 1e5, 1e5, 0.0).sweep(v, 1e-3)[0]
    released = RateEquation(ci, 1e5, 1e5, 1e5).sweep(v, 1e-3)[0]
    equilibrium = ci.mean_occupation(v)
    assert np.abs(released - equilibrium).max() < 0.02
    assert np.abs(stuck - equilibrium).max() > 0.5
