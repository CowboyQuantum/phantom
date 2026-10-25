"""P0: rate-limited loading follows equilibrium when the barrier is open and freezes when it is closed."""

import numpy as np
import pytest

from phantom.readout.loading import LeadBarrier, latched_occupation


def test_lead_barrier_rate_is_exponential():
    b = LeadBarrier("BL", rate0_hz=1e3, v_ref_v=-0.1, v_scale_v=0.03)
    assert b.rate_hz(-0.1) == pytest.approx(1e3)
    assert b.rate_hz(-0.1 + 0.03 * np.log(10)) == pytest.approx(1e4, rel=1e-12)
    with pytest.raises(ValueError):
        LeadBarrier("BL", rate0_hz=0.0)


def test_open_barrier_follows_equilibrium_and_closed_barrier_freezes():
    n_eq = np.stack([np.repeat(np.arange(5.0), 20), np.zeros(100)], axis=1)
    rng = np.random.default_rng(0)
    fast = latched_occupation(n_eq, np.full(100, 1e6), dwell_s=1e-3, rng=rng)
    np.testing.assert_array_equal(fast, n_eq)
    slow = latched_occupation(n_eq, np.full(100, 1e-6), dwell_s=1e-3, rng=rng)
    np.testing.assert_array_equal(slow, np.broadcast_to(n_eq[0], n_eq.shape))


def test_latching_is_seeded_and_intermediate_rates_partly_follow():
    n_eq = np.stack([np.repeat(np.arange(10.0), 50), np.zeros(500)], axis=1)
    a = latched_occupation(n_eq, np.full(500, 20.0), 1e-3, np.random.default_rng(3))
    b = latched_occupation(n_eq, np.full(500, 20.0), 1e-3, np.random.default_rng(3))
    np.testing.assert_array_equal(a, b)
    lag = n_eq[:, 0] - a[:, 0]
    assert lag.min() >= 0 and 0 < lag.mean() < 5, "the measured charge trails equilibrium, never leads it"
