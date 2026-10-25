"""P0: the dynamic lock-in channel has the three limits that define it.

- rate far above the chop frequency, negligible lifetime width: equals the quasi-static channel;
- rate far below the chop frequency: the charge never follows, the signal vanishes (latching);
- rate far above everything: lifetime broadening lowers the derivative peak monotonically.
"""

import numpy as np

from phantom.charge import ConstantInteractionParams
from phantom.readout import LeadRateLandscape, LockinParams, SensedDoubleDot
from phantom.readout.lockin import lockin_sweep_a_per_v, lockin_sweep_dynamic_a_per_v


def _device() -> SensedDoubleDot:
    return SensedDoubleDot(charge=ConstantInteractionParams(temperature_k=0.1))


def _line(d: SensedDoubleDot, n: int = 400) -> np.ndarray:
    v = np.zeros((n, 3))
    v[:, 0] = np.linspace(0.0, 0.04, n)  # three dot-1 transitions
    return v


def _uniform(d: SensedDoubleDot, rate_hz: float) -> LeadRateLandscape:
    return LeadRateLandscape(d.charge.gate_names, rate0_hz=np.array([rate_hz, rate_hz]), weights_per_v=np.zeros((2, 3)))


def test_fast_rate_reproduces_the_quasi_static_channel():
    d, v, lk = _device(), _line(_device()), LockinParams(gate="P1", frequency_hz=1e3, cycles=2)
    static = lockin_sweep_a_per_v(d, v, lk)
    # 1 MHz: settles every half-cycle (p = 1 - exp(-500)), lifetime width 7e-5 of kT at 0.1 K
    dyn = lockin_sweep_dynamic_a_per_v(d, v, lk, _uniform(d, 1e6), np.random.default_rng(0))
    np.testing.assert_allclose(dyn["signal_a_per_v"], static, rtol=1e-6, atol=1e-8 * np.abs(static).max())


def test_slow_rate_latches_and_the_signal_vanishes():
    d, v, lk = _device(), _line(_device()), LockinParams(gate="P1", frequency_hz=1e3, cycles=2)
    static = lockin_sweep_a_per_v(d, v, lk)
    dyn = lockin_sweep_dynamic_a_per_v(d, v, lk, _uniform(d, 1e-3), np.random.default_rng(0))
    assert np.abs(dyn["signal_a_per_v"]).max() < 0.02 * np.abs(static).max()
    n = dyn["n_present_mean"]
    assert np.all(n == n[0]), "the charge present never changes after the first settled step"


def test_lifetime_broadening_lowers_the_peak_monotonically():
    d, v, lk = _device(), _line(_device()), LockinParams(gate="P1", frequency_hz=1e3, cycles=2)
    # a small excitation (0.2 mV, about 30 ueV through the lever arm) so the line width, not the
    # modulation amplitude, sets the derivative peak
    lk = LockinParams(gate="P1", amplitude_v=0.2e-3, frequency_hz=1e3, cycles=2)
    peaks = []
    for rate in (1e7, 3e10, 1e11, 3e11):
        dyn = lockin_sweep_dynamic_a_per_v(d, v, lk, _uniform(d, rate), np.random.default_rng(0))
        peaks.append(np.abs(dyn["signal_a_per_v"]).max())
    assert peaks[0] > peaks[1] > peaks[2] > peaks[3]
    assert peaks[-1] < 0.5 * peaks[0]


def test_direction_vector_and_detuning_sign():
    names = ("P1", "P2", "SG")
    assert list(LockinParams(gate="P2", amplitude_v=2e-3).vector_v(names)) == [0.0, 2e-3, 0.0]
    det = LockinParams(direction={"P1": 1.0, "P2": -1.0}, amplitude_v=1e-3).vector_v(names)
    assert list(det) == [1e-3, -1e-3, 0.0]
    d = _device()
    v = _line(d)
    common = lockin_sweep_a_per_v(d, v, LockinParams(direction={"P1": 1.0, "P2": 1.0}))
    p1 = lockin_sweep_a_per_v(d, v, LockinParams(gate="P1"))
    assert common.shape == p1.shape and not np.array_equal(common, p1)
