"""P0 for the second Step 6 PR: blockade suppresses the forward current and the valley splitting or
field-quenched mixing lifts it; the funnel sits on J = E_Z; Elzerman shots average to the analytic
curve; the filter function has its Ramsey, echo and CPMG limits; the exchange gate is perfect without
noise and degrades with the draw's noise."""

import numpy as np
import pytest

from phantom.charge import Barrier, ConstantInteractionParams
from phantom.charge.hubbard import HubbardParams
from phantom.constants import HBAR_EV_S
from phantom.spin import (
    BlockadeRates,
    SpinParams,
    bias_triangle_current_a,
    ci_energies_ev,
    coherence,
    cpmg_filter,
    detuning_noise_sigma_ev,
    elzerman_shots,
    empty_probability,
    funnel_singlet_probability,
    log_grid_hz,
    omega_spectrum,
    one_over_f_amplitude_sq,
    one_over_f_psd,
    swap_gate_fidelity,
    symmetric_point_v,
    t2_star_s,
)


def _vertex_window(ci: ConstantInteractionParams, half_v: float, n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    tp1, tp2 = ci.triple_points_ng()
    c = ci.gate_voltages_for_ng((tp1 + tp2) / 2 + np.array([0.0, 1.0]))
    p1 = np.linspace(c[0] - half_v, c[0] + half_v, n)
    p2 = np.linspace(c[1] - half_v, c[1] + half_v, n)
    v = np.zeros((n, n, 3))
    v[..., 0], v[..., 1] = p1[:, None], p2[None, :]
    return p1, p2, v


def test_blockade_suppresses_forward_current_and_is_lifted_by_valley_or_mixing():
    ci = ConstantInteractionParams()
    _, _, v = _vertex_window(ci, 0.02, 81)
    rates = BlockadeRates(gamma_lead_hz=1e8, gamma_interdot_hz=1e7, gamma_mix_hz=1e7, b_mix_t=1e-3)
    e_v = 150e-6
    fwd_blocked = bias_triangle_current_a(ci, v, 1e-3, e_v, 0.1, rates)
    fwd_lifted = bias_triangle_current_a(ci, v, 1e-3, e_v, 0.0, rates)
    rev = bias_triangle_current_a(ci, v, -1e-3, e_v, 0.1, rates)
    e01, e11, e02, _ = np.moveaxis(ci_energies_ev(ci, v, ((0, 1), (1, 1), (0, 2), (1, 2))), -1, 0)
    eps = e11 - e02
    inside = (eps > 20e-6) & (eps < e_v - 65e-6) & (e11 - e01 < 0.5e-3 - 50e-6) & (e02 - e01 > -0.5e-3 + 50e-6)
    assert inside.sum() > 5
    assert np.all(fwd_blocked[inside] < 0.1 * fwd_lifted[inside]), "blocked base of the triangle"
    assert np.all(fwd_lifted[inside] > 1e-13)
    beyond = (eps > e_v + 50e-6) & (e11 - e01 < 0.5e-3 - 50e-6) & (e02 - e01 > -0.5e-3 + 50e-6)
    assert beyond.sum() > 5
    assert np.all(fwd_blocked[beyond] > 0.5 * fwd_lifted[beyond]), "valley-excited triplet carries current beyond E_V"
    assert rev.max() > 1e-12 and fwd_blocked.sum() < 0.8 * rev.sum(), "the reverse bias has no blockade, so it carries more total current"
    assert np.all(fwd_blocked >= 0) and np.all(rev >= 0)


def test_funnel_follows_j_equals_zeeman_energy():
    j = np.linspace(0.5e-6, 10e-6, 20001)
    ez = 4e-6
    p = funnel_singlet_probability(j, ez, delta_mix_ev=20e-9)
    assert j[np.argmin(p)] == pytest.approx(ez, abs=0.1e-6)
    assert p.min() == pytest.approx(0.5, abs=1e-3) and p[0] > 0.99 and p[-1] > 0.99
    p_t0 = funnel_singlet_probability(np.array([1e-9, 1e-6]), 1.0, 0.0, delta_zeeman_ev=2e-9)
    assert p_t0[0] < p_t0[1], "a Zeeman difference costs singlet probability where J is small"


def test_elzerman_shots_average_to_the_analytic_empty_probability():
    rng = np.random.default_rng(0)
    g_out, g_in, t1 = 2e4, 1e4, 1e-2
    traces, spins = elzerman_shots(rng, 6000, 1.0, g_out, g_in, t1, 2e-6, 400)
    assert spins.all()
    t = np.arange(400) * 2e-6
    empty = 1 - traces.mean(axis=0)
    np.testing.assert_allclose(empty, empty_probability(t, g_out, g_in, t1), atol=0.03)
    traces_down, _ = elzerman_shots(np.random.default_rng(1), 50, 0.0, g_out, g_in, t1, 2e-6, 400)
    assert traces_down.all(), "a lower-state electron never leaves"
    a, _ = elzerman_shots(np.random.default_rng(5), 20, 0.5, g_out, g_in, t1, 2e-6, 100)
    b, _ = elzerman_shots(np.random.default_rng(5), 20, 0.5, g_out, g_in, t1, 2e-6, 100)
    np.testing.assert_array_equal(a, b)


def test_filter_function_limits_and_cpmg_recovers_coherence_under_one_over_f():
    x = np.linspace(0.01, 50, 500)
    np.testing.assert_allclose(cpmg_filter(x, 0), 4 * np.sin(x / 2) ** 2, atol=1e-12), "Ramsey: |e^{ix} − 1|²"
    np.testing.assert_allclose(cpmg_filter(x, 1), 16 * np.sin(x / 4) ** 4, atol=1e-12), "Hahn echo"
    # quasi-static limit: a spectrum confined far below 1/τ gives the Gaussian envelope with T2* = √2/σ_ω
    f = log_grid_hz(1e-3, 1e1)
    sigma_e = 1e-9
    s_omega = omega_spectrum(one_over_f_psd(f, one_over_f_amplitude_sq(sigma_e, 1e-3, 1e1), 1e-3, 1e1), 1.0)
    tau = 0.3 * t2_star_s(sigma_e)
    assert coherence(tau, 0, f, s_omega) == pytest.approx(np.exp(-((tau / t2_star_s(sigma_e)) ** 2)), rel=2e-3)
    # 1/f over a wide band: echo beats Ramsey, and CPMG improves with N
    f = log_grid_hz(1e-2, 1e8)
    s_omega = omega_spectrum(one_over_f_psd(f, one_over_f_amplitude_sq(1e-9, 1e-2, 1e8), 1e-2, 1e8), 1.0)
    tau = 2e-6
    w = [coherence(tau, n, f, s_omega) for n in (0, 1, 4, 16)]
    assert w[0] < w[1] < w[2] < w[3]
    # white noise: decoupling does not help (χ ∝ τ for every N)
    s_white = np.full(f.shape, 1e3)
    w_white = [coherence(1e-5, n, f, s_white) for n in (1, 4, 16)]
    np.testing.assert_allclose(w_white, w_white[0], rtol=0.05)


def test_exchange_gate_is_perfect_without_noise_and_degrades_with_the_draw():
    hub = HubbardParams.from_ci(ConstantInteractionParams(), Barrier("B", t0_ev=20e-6, v_scale_v=0.05))
    v0 = symmetric_point_v(hub)
    clean = SpinParams(t1_s=1e9, t_phi_markov_s=1e9)
    f0, t_gate = swap_gate_fidelity(hub, v0, 0.0, 0.0, clean)
    assert f0 == pytest.approx(1.0, abs=1e-9)
    assert t_gate == pytest.approx(np.pi * HBAR_EV_S / 1.07e-6, rel=0.02)
    f_a, _ = swap_gate_fidelity(hub, v0, 1.0e-3, detuning_noise_sigma_ev(hub, 100e-6), clean)
    f_b, _ = swap_gate_fidelity(hub, v0, 1.0e-3, detuning_noise_sigma_ev(hub, 200e-6), clean)
    assert 1.0 > f_a > f_b, "more detuning noise, lower fidelity away from the sweet spot"
    f_sweet, _ = swap_gate_fidelity(hub, v0, 0.0, detuning_noise_sigma_ev(hub, 200e-6), clean)
    assert f_sweet > f_b, "the symmetric point is the sweet spot"
    f_dec, _ = swap_gate_fidelity(hub, v0, 0.0, 0.0, SpinParams(t1_s=1e9, t_phi_markov_s=t_gate * 20))
    assert 0.9 < f_dec < 0.999, "Markovian dephasing during the gate"
