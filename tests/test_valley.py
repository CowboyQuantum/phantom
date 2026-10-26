"""P0 for Step 5: the valley-splitting distribution has the stated scalings and limits; the signatures
sit where the splitting says."""

import numpy as np
import pytest

from phantom.valley import (
    MU_B_UEV_PER_T,
    ValleyParams,
    draw_valley_splitting_uev,
    kink_field_t,
    pulsed_spectroscopy_occupation,
    relaxation_rate_hz,
    rice_mean,
    second_electron_shift_uev,
)


def test_sigma_scales_with_dot_size_and_germanium_concentration():
    p = ValleyParams(sigma_ref_uev=45.0, radius_ref_nm=10.0, ge_ref=0.3)
    assert p.sigma_uev(10.0) == pytest.approx(45.0)
    assert p.sigma_uev(20.0) == pytest.approx(22.5), "inverse in the radius"
    assert p.sigma_uev(10.0, ge=0.5) == pytest.approx(45.0 * np.sqrt(0.25 / 0.21))
    assert p.sigma_uev(10.0, interface_width_nm=4.0) == pytest.approx(90.0), "wider interface, more overlap"
    assert p.sigma_uev(10.0, step_factor=1.5) == pytest.approx(67.5)


def test_rayleigh_and_rice_limits_of_the_mean():
    assert rice_mean(0.0, 10.0) == pytest.approx(10.0 * np.sqrt(np.pi / 2), rel=1e-3)
    assert rice_mean(500.0, 10.0) == pytest.approx(500.0 + 10.0**2 / (2 * 500.0), rel=2e-3), "large deterministic part"
    assert rice_mean(50.0, 0.0) == 50.0


def test_draws_follow_the_distribution_and_the_scalings():
    rng = np.random.default_rng(0)
    p = ValleyParams(delta_det_uev=0.0, sigma_ref_uev=45.0)
    draws = draw_valley_splitting_uev(p, np.full(20000, 10.0), rng)
    assert draws.mean() == pytest.approx(45.0 * np.sqrt(np.pi / 2), rel=0.03)
    big = draw_valley_splitting_uev(p, np.full(20000, 20.0), np.random.default_rng(1))
    assert big.mean() == pytest.approx(draws.mean() / 2, rel=0.05), "larger dots: half the splitting"
    mos = draw_valley_splitting_uev(ValleyParams.for_material("MOS"), np.full(20000, 10.0), np.random.default_rng(2))
    assert mos.mean() > 300 and mos.std() / mos.mean() < 0.2, "MOS: large and deterministic-dominated"
    a = draw_valley_splitting_uev(p, np.full(5, 10.0), np.random.default_rng(7))
    b = draw_valley_splitting_uev(p, np.full(5, 10.0), np.random.default_rng(7))
    np.testing.assert_array_equal(a, b)


def test_magnetospectroscopy_kink_sits_at_the_valley_splitting():
    e_v = 100.0
    b = np.linspace(0, 3, 3001)
    shift = second_electron_shift_uev(e_v, b)
    b_star = kink_field_t(e_v)
    assert b_star == pytest.approx(100.0 / (2 * MU_B_UEV_PER_T))
    k = int(np.argmax(shift))
    assert b[k] == pytest.approx(b_star, abs=2e-3)
    assert shift[k] == pytest.approx(e_v / 2, rel=1e-2)
    slope_before = (shift[200] - shift[100]) / (b[200] - b[100])
    slope_after = (shift[-1] - shift[-100]) / (b[-1] - b[-100])
    assert slope_before == pytest.approx(MU_B_UEV_PER_T) and slope_after == pytest.approx(-MU_B_UEV_PER_T)


def test_hot_spot_peaks_at_the_kink_field():
    b = np.linspace(0, 3, 3001)
    rate = relaxation_rate_hz(b, 100.0, background_hz=2.0, peak_hz=1e4, width_uev=4.0)
    assert b[np.argmax(rate)] == pytest.approx(kink_field_t(100.0), abs=2e-3)
    assert rate.max() == pytest.approx(1e4 + 2.0, rel=1e-3)
    assert rate[0] == pytest.approx(2.0 + 1e4 * 16.0 / (100.0**2 + 16.0), rel=1e-6), "background plus the Lorentzian tail at zero field"


def test_pulsed_spectroscopy_lines():
    v = np.linspace(-0.02, 0.02, 801)
    a = np.array([0.0, 0.01])
    n = pulsed_spectroscopy_occupation(v, a, 0.0, 200.0, 0.2, kt_uev=5.0, excited_weight=0.5)
    np.testing.assert_allclose(n[0], 1.0 / (1.0 + np.exp(-v / (5e-6 / 0.2))), atol=1e-12), "no pulse: the equilibrium step"
    # with a pulse: the plateau between -A and 0 holds half an electron, and the excited-state line sits
    # E_V/alpha = 1 mV below the pulsed line with half the pulsed half's weight
    assert n[1][np.argmin(np.abs(v + 0.005))] == pytest.approx(0.5, abs=0.01)
    at_excited = n[1][np.argmin(np.abs(v + 0.0105))]
    assert 0.2 < at_excited < 0.3
