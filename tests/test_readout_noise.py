"""P0: noise processes have their declared statistics, are seeded, and `trace` equals `sample`."""

import itertools

import numpy as np
import pytest

from phantom.readout import (
    Composite,
    NoiseProcess,
    NoNoise,
    OneOverFNoise,
    TelegraphNoise,
    WhiteNoise,
)


def test_processes_satisfy_the_protocol():
    for p in (WhiteNoise(1.0), TelegraphNoise(1.0, 1.0, 1.0), OneOverFNoise(1.0), NoNoise()):
        assert isinstance(p, NoiseProcess)


def test_white_noise_has_the_declared_sigma_and_is_reproducible():
    a = WhiteNoise(2.0, seed=3).trace(np.arange(20000))
    b = WhiteNoise(2.0, seed=3).trace(np.arange(20000))
    np.testing.assert_array_equal(a, b)
    assert np.std(a) == pytest.approx(2.0, rel=0.03)
    assert WhiteNoise(0.0).sample(0.0) == 0.0


def test_telegraph_takes_two_values_and_switches_at_its_rates():
    tn = TelegraphNoise(amplitude=1.0, rate_up_hz=10.0, rate_down_hz=10.0, seed=0)
    t = np.arange(0, 200.0, 1e-3)
    x = tn.trace(t)
    assert set(np.unique(x)) == {-0.5, 0.5}
    assert np.count_nonzero(np.diff(x)) == pytest.approx(200 * 10.0, rel=0.15)
    assert x.mean() == pytest.approx(0.0, abs=0.05)


def test_telegraph_is_sampled_exactly_for_any_time_step():
    tn = TelegraphNoise(amplitude=1.0, rate_up_hz=3.0, rate_down_hz=1.0, seed=5)
    x = tn.trace(1e6 * np.arange(4000))
    assert x.mean() == pytest.approx(0.75 - 0.5, abs=0.03)


def test_one_over_f_spectrum_slope():
    n = OneOverFNoise(amplitude=1.0, f_min_hz=1.0, f_max_hz=1e3, n_fluctuators=40, seed=2)
    dt, length = 1e-4, 2**17
    x = n.trace(np.arange(length) * dt)
    assert np.std(x) == pytest.approx(1.0, rel=0.3)
    psd = np.abs(np.fft.rfft(x - x.mean())) ** 2
    f = np.fft.rfftfreq(length, dt)
    fb, pb = [], []
    for lo, hi in itertools.pairwise(np.geomspace(3.0, 300.0, 12)):
        sel = (f >= lo) & (f < hi)
        fb.append(np.sqrt(lo * hi))
        pb.append(psd[sel].mean())
    slope = np.polyfit(np.log(fb), np.log(pb), 1)[0]
    assert -1.4 < slope < -0.6


@pytest.mark.parametrize(
    "make",
    [
        lambda: WhiteNoise(1.5, seed=4),
        lambda: TelegraphNoise(1.0, 2.0, 5.0, seed=4),
        lambda: OneOverFNoise(1.0, 1e-1, 1e3, n_fluctuators=8, seed=4),
        lambda: Composite((WhiteNoise(1.0, seed=1), TelegraphNoise(1.0, 1.0, 1.0, seed=2))),
    ],
)
def test_trace_equals_successive_samples(make):
    t = np.sort(np.random.default_rng(0).uniform(0, 10, 500))
    t[100] = t[99]  # a repeated time consumes no random number
    one = make()
    by_sample = np.array([one.sample(ti) for ti in t])
    by_trace = make().trace(t)
    np.testing.assert_array_equal(by_sample, by_trace)
    # and mixing the two on one object continues the same stream
    mixed = make()
    head = np.array([mixed.sample(ti) for ti in t[:50]])
    tail = mixed.trace(t[50:])
    np.testing.assert_array_equal(np.concatenate([head, tail]), by_trace)


def test_invalid_rates_rejected():
    with pytest.raises(ValueError):
        TelegraphNoise(1.0, 0.0, 0.0)
    with pytest.raises(ValueError):
        OneOverFNoise(1.0, f_min_hz=10.0, f_max_hz=1.0)
