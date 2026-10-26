"""P0 for the RF tank and single-shot readout: resonance and matching, phase flip, square-root law."""

import numpy as np
import pytest

from phantom.readout import (
    NoNoise,
    SensedDoubleDot,
    TankCircuit,
    WhiteNoise,
    noisy_device_for_sampling,
    sensor_resistance_ohm,
    single_shots,
)


def test_tank_resonance_and_matching():
    t = TankCircuit(inductance_h=820e-9, parasitic_capacitance_f=0.5e-12, z0_ohm=50.0)
    assert t.resonance_hz == pytest.approx(1 / (2 * np.pi * np.sqrt(820e-9 * 0.5e-12)))
    f0 = t.resonance_hz
    assert abs(t.reflection(f0, t.matched_resistance_ohm)) < 0.03, "matched: almost no reflection at resonance"
    assert abs(t.reflection(f0, 1e9)) > 0.99, "open load: total reflection"
    f = np.linspace(0.8 * f0, 1.2 * f0, 2001)
    mag = np.abs(t.reflection(f, 1e5))
    assert abs(f[np.argmin(mag)] - f0) < 0.01 * f0, "the dip sits at the resonance"
    below, above = t.reflection(f0, 0.5 * t.matched_resistance_ohm), t.reflection(f0, 2 * t.matched_resistance_ohm)
    assert abs(np.degrees(np.angle(below)) - np.degrees(np.angle(above))) > 150, "phase flips across the match"
    with pytest.raises(ValueError):
        TankCircuit(inductance_h=0.0)


def test_sensor_resistance_is_bias_over_current_and_clipped():
    r = sensor_resistance_ohm(np.array([1e-9, 1e-15, 1.0]), 50e-6)
    assert r[0] == pytest.approx(5e4) and r[1] == 1e9 and r[2] == 1e3


def _white_only(sigma_ref: float, seed: int) -> SensedDoubleDot:
    return SensedDoubleDot(readout_noise=WhiteNoise(sigma_ref, seed=seed), gate_noise=NoNoise(), sensor_noise=NoNoise())


def test_snr_follows_the_square_root_law_with_white_noise_only():
    base = _white_only(15e-12, 0)
    v = np.array([0.0, 0.0, base.sensor.phase_v + 0.9 * base.sensor.width_v])
    n0 = base.charge.mean_occupation(v)
    states = np.array([n0, n0 + [1.0, 0.0]])
    snr = []
    for k, tau in enumerate((1e-4, 1e-3, 1e-2)):
        dt = tau / 20
        dev = noisy_device_for_sampling(base, dt, seed=10 + k)
        res = single_shots(dev, v, states, tau, dt, 600, np.random.default_rng(k))
        snr.append(res.snr)
    assert snr[1] / snr[0] == pytest.approx(np.sqrt(10), rel=0.2)
    assert snr[2] / snr[1] == pytest.approx(np.sqrt(10), rel=0.2)
    res = single_shots(noisy_device_for_sampling(base, 5e-4, 99), v, states, 1e-2, 5e-4, 600, np.random.default_rng(5))
    assert res.fidelity > 0.99


def test_white_noise_rescales_with_the_sampling_interval():
    base = _white_only(15e-12, 0)
    fine = noisy_device_for_sampling(base, 1e-5, 1)
    assert fine.readout_noise.sigma == pytest.approx(15e-12 * np.sqrt(1e-3 / 1e-5))


def test_shots_are_seeded():
    base = SensedDoubleDot.default_noisy(seed=0)
    v = np.array([0.0, 0.0, base.sensor.phase_v + 0.9 * base.sensor.width_v])
    n0 = base.charge.mean_occupation(v)
    states = np.array([n0, n0 + [1.0, 0.0]])
    a = single_shots(noisy_device_for_sampling(base, 5e-5, 3), v, states, 1e-3, 5e-5, 100, np.random.default_rng(2))
    b = single_shots(noisy_device_for_sampling(base, 5e-5, 3), v, states, 1e-3, 5e-5, 100, np.random.default_rng(2))
    np.testing.assert_array_equal(a.signals_a, b.signals_a)
    assert set(np.unique(a.true_state)) <= {0, 1}
