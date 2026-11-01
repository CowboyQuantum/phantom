"""P0 for Step 7: a compiled device round-trips through serialisation text-identically and its runtime
objects reproduce its diagrams bit for bit; the clean device matches the Step 2 derivation; the FEM
hook refuses with the open decision's name; every `optimum` is attained by brute force."""

import numpy as np
import pytest

from phantom.charge import derive_double_dot
from phantom.compile import DEVICE_FORMAT_VERSION, CompileParams, Device, compile_device, fem_basis
from phantom.disorder import ScenarioParams, draw_scenario
from phantom.layout import two_dot_layout
from phantom.optimum import (
    coherence_optimum,
    coupling_optimum,
    optimum,
    readout_optimum,
    sensor_optimum,
    virtual_gate_matrix,
)
from phantom.spin import psb_contrast, swap_gate_fidelity
from phantom.spin.exchange import detuning_noise_sigma_ev

OP = {"B_L": 0.15, "P1": 0.6, "B_M": 0.2, "P2": 0.6, "B_R": 0.15, "SG": 0.6, "R_L": 0.7, "R_R": 0.7, "SB_L": 0.4, "SB_R": 0.3, "SR_L": 0.7, "SR_R": 0.7}


@pytest.fixture(scope="module")
def device() -> Device:
    return compile_device(two_dot_layout("SiGe", 80.0), CompileParams(operating_point_v=OP), seed=2)


def _coarse_map(dev: Device) -> np.ndarray:
    sensed = dev.runtime_sensed(noise=False)
    v0 = dev.operating_point()
    i1, i2 = dev.gate_names.index("P1"), dev.gate_names.index("P2")
    p = np.linspace(-0.03, 0.03, 25)
    vv = np.broadcast_to(v0, (25, 25, len(v0))).copy()
    vv[..., i1] += p[:, None]
    vv[..., i2] += p[None, :]
    return sensed.sweep_a(vv, None)


def test_round_trip_is_text_identical_and_reproduces_the_diagram(device: Device, tmp_path):
    path = tmp_path / "device.json"
    device.to_json(path)
    again = Device.from_json(path)
    assert again.dumps() == device.dumps()
    assert again.provenance["format_version"] == DEVICE_FORMAT_VERSION and "charge.ci" in again.provenance["model_versions"]
    np.testing.assert_array_equal(_coarse_map(again), _coarse_map(device))
    assert path.stat().st_size < 100_000, "a device file stays small: the basis is not stored (the shape grid, #59, is about 40 kB)"
    _lay, basis = again.basis()
    assert basis.basis.shape[0] == len(device.gate_names)


def test_compile_is_deterministic_in_the_seed_and_the_draws_differ_between_seeds():
    lay = two_dot_layout("SiGe", 80.0)
    a = compile_device(lay, CompileParams(operating_point_v=OP), seed=4)
    b = compile_device(lay, CompileParams(operating_point_v=OP), seed=4)
    c = compile_device(lay, CompileParams(operating_point_v=OP), seed=5)
    assert a.dumps() == b.dumps()
    assert a.valley_splittings_uev().tolist() != c.valley_splittings_uev().tolist()
    assert a.g_factors != c.g_factors and a.scenario["charges"] != c.scenario["charges"]


def test_clean_device_matches_the_step_2_derivation_and_carries_every_draw():
    lay = two_dot_layout("SiGe", 80.0)
    dev = compile_device(lay, CompileParams(operating_point_v=OP, severity=0.0, wiring_faults=False), seed=0)
    v = dev.operating_point()
    derived = derive_double_dot(lay, v, temperature_k=0.15)
    np.testing.assert_allclose(np.array(dev.c_dot_af), derived.ci.c_dot_af, rtol=1e-12)
    np.testing.assert_allclose(np.array(dev.c_gate_af), derived.ci.c_gate_af, rtol=1e-12)
    assert dev.programmed_point_v[dev.gate_names.index("P1")] == 0.6
    assert dev.operating_point_v[dev.gate_names.index("SG")] != 0.6, "the sensor gate moved onto the flank"
    assert len(dev.dots) == 2 and dev.sensor_dot is not None
    assert all(d.valley_splitting_uev > 0 for d in dev.dots) and len(dev.g_factors) == 2
    assert {b.between for b in dev.barriers} == {("dot1", "lead_left"), ("dot2", "lead_right"), ("dot1", "dot2")}
    assert dev.tunnel_coupling_ev > 0 and dev.barrier_v_scale_v > 0
    hub = dev.runtime_hubbard()
    assert hub.barrier_gate == "B_M" and hub.barrier.tunnel_coupling_ev(v[dev.gate_names.index("B_M")]) == pytest.approx(dev.tunnel_coupling_ev)
    assert dev.runtime_rates().gamma_left_hz == dev.barriers[0].rate_hz
    assert dev.scenario["noise"] is not None and dev.noise_sigma_gate_v() == pytest.approx(dev.scenario["noise"]["one_over_f_amplitude_v"])


def test_wiring_faults_map_the_programmed_point(device: Device):
    lay = two_dot_layout("SiGe", 80.0)
    draw = draw_scenario(lay, ScenarioParams.for_material("SiGe"), 2, 1.0, True)
    v_prog = np.array([OP[g] for g in lay.gate_names])
    expected = draw.wiring.apply(v_prog)
    actual = device.effective_voltages(v_prog)
    np.testing.assert_allclose(actual, expected)


def test_fem_backend_names_the_open_decision():
    lay = two_dot_layout("SiGe", 80.0)
    with pytest.raises(NotImplementedError, match="fem"):
        fem_basis(lay, draw_scenario(lay, ScenarioParams.for_material("SiGe"), 0, 0.0), 2.0)


def test_sensor_optimum_is_the_steepest_point_of_a_dense_scan(device: Device):
    opt = sensor_optimum(device)
    sensed = device.runtime_sensed(noise=False)
    j = device.gate_names.index("SG")
    v0 = device.operating_point()
    period = device.runtime_sensor().period_v / device.runtime_sensor().alpha_v_per_v[j]
    grid = v0[j] + np.linspace(-period / 2, period / 2, 6001)  # the same peak period the optimum searches
    vv = np.broadcast_to(v0, (len(grid), len(v0))).copy()
    vv[:, j] = grid
    i = sensed.sweep_a(vv, None)
    slope = np.abs(np.gradient(i, grid))
    k = int(np.argmax(slope))
    at_opt = slope[int(np.argmin(np.abs(grid - opt.v_sensor_gate_v)))]
    assert at_opt >= 0.999 * slope[k], "the chosen point is as steep as the steepest point of a dense scan (the two flanks of a peak tie)"
    assert opt.slope_a_per_v >= 0.999 * slope[k]


def test_virtual_gates_decouple_the_dots(device: Device):
    m = virtual_gate_matrix(device)
    lever = device.runtime_ci().lever_arms_ev_per_v()
    ip = [device.gate_names.index("P1"), device.gate_names.index("P2")]
    np.testing.assert_allclose(lever @ m[:, ip], np.diag(np.diag(lever[:, ip])), atol=1e-12), "virtual plunger i moves dot i alone, at the raw scale"
    js = device.gate_names.index("SG")
    alpha = np.array(device.sensor_alpha_v_per_v)
    for j in range(len(device.gate_names)):
        if j not in ip and j != js:
            np.testing.assert_allclose(lever @ m[:, j], 0.0, atol=1e-12), "other gates' cross-talk on the dots is compensated"
            np.testing.assert_allclose(alpha @ m[:, j], 0.0, atol=1e-12), "and on the sensor"
    np.testing.assert_allclose(alpha @ m[:, ip], 0.0, atol=1e-12), "the sensor potential does not move with the virtual plungers"
    m_raw_sensor = virtual_gate_matrix(device, compensate_sensor=False)
    np.testing.assert_allclose(lever @ m_raw_sensor[:, ip], np.diag(np.diag(lever[:, ip])), atol=1e-12)


def test_coupling_optimum_inverts_the_runtime_coupling(device: Device):
    c = coupling_optimum(device, 20e-6)
    v = device.operating_point()
    v[device.gate_names.index(device.barrier_gate)] = c.v_barrier_v
    assert device.runtime_coupling().tunnel_coupling_ev(v) == pytest.approx(20e-6, rel=1e-9)
    assert np.isnan(coupling_optimum(device, 2 * device.runtime_coupling().t_max_ev).v_barrier_v), "beyond the ceiling"


def test_readout_optimum_matches_brute_force_over_detuning(device: Device):
    r = readout_optimum(device)
    kt = 8.617333262e-5 * device.temperature_k
    e_v = device.dots[1].valley_splitting_uev * 1e-6
    eps = np.linspace(-300e-6, 500e-6, 20001)
    c = psb_contrast(eps, e_v, device.tunnel_coupling_ev, kt)
    assert r.psb_contrast == pytest.approx(c.max(), rel=1e-3)
    assert abs(r.psb_detuning_ev - eps[np.argmax(c)]) < 1e-6
    assert 0 < r.elzerman_peak_time_s and 0 < r.elzerman_peak_probability <= 1


def test_coherence_optimum_is_the_sweet_spot_of_a_brute_force_scan(device: Device):
    c = coherence_optimum(device)
    hub = device.runtime_hubbard()
    v_sweet = np.array(c.sweet_spot_v)
    occ = hub.ground_state_occupation(v_sweet)
    np.testing.assert_allclose(occ, np.rint(occ), atol=1e-6)
    sigma = detuning_noise_sigma_ev(hub, device.noise_sigma_gate_v())
    spin = __import__("phantom.spin", fromlist=["SpinParams"]).SpinParams(**device.params["spin"])
    eps = np.linspace(-0.4e-3, 0.4e-3, 9)
    f = np.array([swap_gate_fidelity(hub, v_sweet, float(e), sigma, spin)[0] for e in eps])
    assert f[4] == pytest.approx(c.swap_fidelity_at_sweet_spot, rel=1e-9) and f[4] >= f.max() - 1e-9
    assert c.t2_star_single_spin_s > 0
    o = optimum(device)
    assert o.to_dict()["coherence"]["swap_fidelity_at_sweet_spot"] == c.swap_fidelity_at_sweet_spot


def test_action_rate_is_never_below_the_rectangular_bound():
    from phantom.charge import derive_double_dot as _derive
    from phantom.electrostatics import GateBasis, barrier_action, barrier_width_nm
    from phantom.readout import wkb_rate_from_action_hz, wkb_tunnel_rate_hz
    from phantom.readout.rates import attempt_rate_hz

    lay = two_dot_layout("SiGe", 80.0)
    basis = GateBasis.from_layout(lay)
    res = lay.reservoir_points_nm()
    for b_l in (-0.1, 0.0, 0.1):
        op = {**OP, "B_L": b_l}
        v = np.array([op[g] for g in lay.gate_names])
        d1 = _derive(lay, v, basis=basis).dots[0]
        h, w = barrier_width_nm(basis, v, d1, res["R_L"])
        _, _, s = barrier_action(basis, v, d1, res["R_L"], 0.19)
        rect = float(wkb_tunnel_rate_hz(h, w, d1.curvature_ev_per_nm2, 0.19))
        act = float(wkb_rate_from_action_hz(s, d1.curvature_ev_per_nm2, 0.19))
        assert act >= rect, "a shaped barrier transmits at least as much as the rectangle that bounds it"
        assert act == pytest.approx(attempt_rate_hz(d1.curvature_ev_per_nm2, 0.19) * np.exp(-2 * s), rel=1e-12)


def test_coupling_optimum_agrees_with_the_compile_tier_near_the_operating_point(device: Device):
    """The optimum's B_M, put back through the compile tier, gives the target: t = (ħω₀/π)·√T with T
    Kemble's transmission of the interdot barrier's top recomputed from the potential with the dots held,
    moved from the compiled top by the barrier gate's screening factor (#70, ADR 0034). The runtime's top
    along B_M is the shape grid's, with the dots moved, which differs by about 1 % in slope."""
    from phantom.charge.barrier import kemble_transmission
    from phantom.compile.device import interdot_barrier_top

    lay, basis = device.basis()
    v0 = device.operating_point()
    d1, d2 = derive_double_dot(lay, v0, basis=basis).dots
    jb = device.gate_names.index(device.barrier_gate)
    level = device.barriers[2].tunnelling_level_ev
    sigma = device.interdot_barrier_screening[jb]
    t_max = device.runtime_coupling().t_max_ev
    for target in (10e-6, 30e-6, 100e-6):
        v = v0.copy()
        v[jb] = coupling_optimum(device, target).v_barrier_v
        top = device.interdot_barrier_top_ev + sigma * (interdot_barrier_top(basis, v, d1, d2, level, 0.19)[0] - device.interdot_barrier_top_ev)
        t = t_max * np.sqrt(kemble_transmission(top, device.interdot_barrier_hbar_omega_ev))
        assert t == pytest.approx(target, rel=0.15), "the compiled transmission holds within a decade of the operating point"


def test_the_compiled_sensor_carries_the_coupling_asymmetry_and_sits_on_its_flank(device: Device):
    """#43, ADR 0016, 0026: each dot's per-electron sensor coupling is the layout's times the wiring
    draw's factor, and compile places the set point on the flank of that same sensor."""
    from phantom.readout.sensor import WING_WIDTHS

    factors = np.array(device.scenario["wiring"]["sensor_beta_factors"])
    assert not np.allclose(factors, 1.0), "the gallery device draws an asymmetry"
    sensor = device.runtime_sensor()
    np.testing.assert_array_equal(sensor.beta_v_per_e, np.array(device.sensor_beta_v_per_e) * factors)
    v = device.operating_point()
    phi = float(sensor.potential(v, device.runtime_ci().mean_occupation(v)))
    offset = (phi - sensor.phase_v) / sensor.period_v
    offset = (offset - np.rint(offset)) * sensor.period_v
    # the target is 0.9 of the sech²-equivalent width at the rates one step before; the rates at the
    # stored point differ by the last step's move along the sensor gate, so the width does too
    width = sensor.wing_offset_v(v) / WING_WIDTHS
    assert 0.5 < offset / (device.params["sensor_flank_widths"] * width) < 1.5
    clean = compile_device(two_dot_layout("SiGe", 80.0), CompileParams(operating_point_v=OP, severity=0.0, wiring_faults=False), seed=0)
    np.testing.assert_array_equal(clean.runtime_sensor().beta_v_per_e, np.array(clean.sensor_beta_v_per_e))
