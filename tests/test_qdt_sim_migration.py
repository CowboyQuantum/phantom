"""The migrated core reproduces qdt-sim's defaults bit for bit (Step 1 exit criterion, one-off check).

`tests/fixtures/qdt_sim_defaults.npz` was generated from qdt-sim's `DoubleDotModel()` and noise
defaults in qdt's own environment; it is synthetic and names the qdt commit. phantom never imports
qdt (rule 1): the fixture is the only link, and Step 8 retires it when qdt's adapter pins phantom.
"""

import pathlib

import numpy as np

from phantom.readout import OneOverFNoise, SensedDoubleDot, TelegraphNoise, WhiteNoise

FIXTURE = pathlib.Path(__file__).parent / "fixtures" / "qdt_sim_defaults.npz"


def load():
    return np.load(FIXTURE)


def test_fixture_is_synthetic():
    assert "qdt-sim" in str(load()["provenance"]) and "no device data" in str(load()["provenance"])


def test_stability_map_and_occupation_match():
    f = load()
    d = SensedDoubleDot()
    v = np.zeros((24, 32, 3))
    v[..., 0] = f["p1_v"][:, None]
    v[..., 1] = f["p2_v"][None, :]
    np.testing.assert_allclose(d.charge.mean_occupation(v), f["mean_occupation"], rtol=0, atol=1e-12)
    np.testing.assert_allclose(d.current_a(v), f["sensor_current_a"], rtol=1e-12, atol=0)
    np.testing.assert_allclose(d.charge.charge_lattice(), f["charge_lattice"], rtol=1e-15)


def test_sensor_diamonds_match():
    f = load()
    d = SensedDoubleDot()
    got = d.sensor.current_a(f["phi_v"][:, None], f["v_sd_v"][None, :])
    np.testing.assert_allclose(got, f["sensor_diamonds_a"], rtol=1e-12, atol=0)


def test_noise_sequences_match_exactly():
    f = load()
    t = f["t_s"]
    np.testing.assert_array_equal(WhiteNoise(15e-12, seed=2).trace(t), f["white_a"])
    np.testing.assert_array_equal(TelegraphNoise(4e-3, 0.5, 0.5, seed=1).trace(t), f["telegraph_v"])
    np.testing.assert_array_equal(OneOverFNoise(150e-6, 1e-2, 1e4, seed=0).trace(t), f["one_over_f_v"])


def test_default_noisy_sweep_is_reproducible_and_seed_dependent():
    v = np.zeros((16, 16, 3))
    v[..., 0] = np.linspace(-0.05, 0.1, 16)[:, None]
    v[..., 1] = np.linspace(-0.05, 0.1, 16)[None, :]
    times = np.arange(256).reshape(16, 16) * 1e-3
    a = SensedDoubleDot.default_noisy(seed=3).sweep_a(v, times)
    b = SensedDoubleDot.default_noisy(seed=3).sweep_a(v, times)
    c = SensedDoubleDot.default_noisy(seed=4).sweep_a(v, times)
    assert a.tobytes() == b.tobytes()
    assert a.tobytes() != c.tobytes()
