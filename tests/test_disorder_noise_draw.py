"""P0: per-device noise draws are seeded, serialisable, clean at zero severity, and build the processes."""

import numpy as np
import pytest

from phantom.disorder import NoiseDraw, NoiseParams, ScenarioParams, draw_noise, draw_scenario
from phantom.layout import two_dot_layout
from phantom.readout import Composite, NoNoise, OneOverFNoise, SensedDoubleDot, WhiteNoise


def test_zero_severity_leaves_only_the_white_floor():
    d = draw_noise(NoiseParams(), np.random.default_rng(0), seed=0, severity=0.0)
    assert d.one_over_f_amplitude_v == 0.0 and d.switcher_amplitudes_v == ()
    dev = SensedDoubleDot.from_noise_draw(d)
    assert isinstance(dev.gate_noise, NoNoise) and isinstance(dev.sensor_noise, NoNoise)
    assert isinstance(dev.readout_noise, WhiteNoise)


def test_draws_are_seeded_round_trip_and_follow_the_parameters():
    params = NoiseParams(switchers_mean=2.0)
    a = draw_noise(params, np.random.default_rng(5), seed=5)
    b = draw_noise(params, np.random.default_rng(5), seed=5)
    assert a.to_dict() == b.to_dict() == NoiseDraw.from_dict(a.to_dict()).to_dict()
    counts = [len(draw_noise(params, np.random.default_rng(s), seed=s).switcher_amplitudes_v) for s in range(400)]
    assert np.mean(counts) == pytest.approx(2.0, rel=0.15)
    amps = [draw_noise(NoiseParams(), np.random.default_rng(s), seed=s).one_over_f_amplitude_v for s in range(400)]
    assert np.exp(np.mean(np.log(amps))) == pytest.approx(150e-6, rel=0.1)
    for d in (a, b):
        for up, down in zip(d.switcher_rates_up_hz, d.switcher_rates_down_hz, strict=True):
            assert 0.1 <= up <= 100 and 0.1 <= down <= 100


def test_processes_from_a_draw_have_the_drawn_statistics():
    d = draw_noise(NoiseParams(switchers_mean=3.0), np.random.default_rng(2), seed=2)
    dev = SensedDoubleDot.from_noise_draw(d, sample_dt_s=1e-4)
    assert isinstance(dev.gate_noise, OneOverFNoise) and dev.gate_noise.f_min_hz == 1e-4
    assert isinstance(dev.sensor_noise, Composite | NoNoise)
    assert dev.readout_noise.sigma == pytest.approx(15e-12 * np.sqrt(10.0))
    x = dev.gate_noise.trace(np.arange(200000) * 1e-3)
    assert np.std(x) == pytest.approx(d.one_over_f_amplitude_v, rel=0.4)


def test_scenario_draw_carries_the_noise_draw():
    s = draw_scenario(two_dot_layout(), ScenarioParams.for_material("MOS"), seed=1)
    assert s.noise is not None and s.noise.one_over_f_amplitude_v > 0
    assert draw_scenario(two_dot_layout(), ScenarioParams.for_material("MOS"), seed=1, severity=0.0).noise.one_over_f_amplitude_v == 0.0
