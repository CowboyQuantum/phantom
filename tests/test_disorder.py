"""P0 for Step 3: the disorder generators have the documented sign, scale and extent; zero severity
reproduces the clean device bit for bit; draws are seeded and serialisable."""

import numpy as np
import pytest

from phantom.charge import derive_double_dot
from phantom.disorder import (
    FixedChargeDraw,
    FixedChargeParams,
    MisalignmentDraw,
    MisalignmentParams,
    ScenarioDraw,
    ScenarioParams,
    apply_misalignment,
    apply_scenario,
    charge_energy_ev,
    draw_fixed_charges,
    draw_misalignment,
    draw_scenario,
)
from phantom.electrostatics import GateBasis
from phantom.layout import two_dot_layout

OP = np.array([0.2, 0.6, 0.2, 0.6, 0.2, 0.6, 0.7, 0.7])


def test_zero_severity_reproduces_the_clean_device_bit_for_bit():
    lay = two_dot_layout()
    draw = draw_scenario(lay, ScenarioParams.for_material("SiGe"), seed=11, severity=0.0)
    perturbed, basis = apply_scenario(lay, draw)
    assert perturbed == lay and basis.offset_ev is None
    clean = GateBasis.from_layout(lay)
    assert basis.basis.tobytes() == clean.basis.tobytes()
    a = derive_double_dot(lay, OP, basis=basis).ci
    b = derive_double_dot(lay, OP, basis=clean).ci
    assert a.c_gate_af.tobytes() == b.c_gate_af.tobytes() and a.c_dot_af.tobytes() == b.c_dot_af.tobytes()


def test_a_positive_charge_lowers_the_electron_energy_locally_with_the_image_scale_and_extent():
    draw = FixedChargeDraw(np.array([0.0]), np.array([0.0]), np.array([1.0]), depth_nm=5.0, smoothing_nm=3.0)
    x = np.linspace(-200, 200, 401)
    u = charge_energy_ev(draw, x, np.array([0.0]), electron_depth_nm=44.0, epsilon_r=11.7)[0]
    assert u.min() < 0 and np.argmin(u) == 200, "attractive, centred on the charge"
    prefactor = 1.602176634e-19 / (4 * np.pi * 8.8541878128e-12 * 11.7) * 1e9  # eV·nm
    expected_peak = -prefactor * (1 / np.sqrt(39.0**2 + 3.0**2) - 1 / 49.0)
    assert u[200] == pytest.approx(expected_peak, rel=1e-12)
    half = np.flatnonzero(u < u.min() / 2)
    extent = x[half[-1]] - x[half[0]]
    assert 30 < extent < 150, f"half-depth extent {extent:.0f} nm is of the order of the depth"
    assert abs(u[0]) < 0.05 * abs(u.min()), "the image makes it fall off faster than a bare Coulomb tail"
    negative = FixedChargeDraw(np.array([0.0]), np.array([0.0]), np.array([-1.0]), 5.0, 3.0)
    np.testing.assert_allclose(charge_energy_ev(negative, x, np.array([0.0]), 44.0, 11.7)[0], -u)


def test_charge_count_follows_density_area_and_severity():
    params = FixedChargeParams(density_cm2=3e10)
    region = (-500.0, 500.0, -500.0, 500.0)  # 1 µm² = 1e-8 cm² → 300 charges expected
    counts = [draw_fixed_charges(params, region, np.random.default_rng(s)).count for s in range(200)]
    assert np.mean(counts) == pytest.approx(300.0, rel=0.05)
    halves = [draw_fixed_charges(params, region, np.random.default_rng(s), severity=0.5).count for s in range(200)]
    assert np.mean(halves) == pytest.approx(150.0, rel=0.08)
    assert draw_fixed_charges(params, region, np.random.default_rng(0), severity=0.0).count == 0
    d = draw_fixed_charges(params, region, np.random.default_rng(3))
    assert np.all((d.x_nm >= -500) & (d.x_nm <= 500)) and set(np.unique(d.sign)) <= {-1.0, 1.0}


def test_misalignment_shifts_and_dilates_exactly_for_rectangles():
    lay = two_dot_layout()
    draw = MisalignmentDraw((4.0, -3.0), {g.name: (1.5 if g.name == "P1" else 0.0) for g in lay.gates})
    moved = apply_misalignment(lay, draw)
    for g, m in zip(lay.gates, moved.gates, strict=True):
        c0, c1 = np.array(g.centroid_nm), np.array(m.centroid_nm)
        np.testing.assert_allclose(c1 - c0, [4.0, -3.0], atol=1e-9)
    p1, m1 = lay.gate("P1").vertices, moved.gate("P1").vertices
    width0 = p1[:, 0].max() - p1[:, 0].min()
    width1 = m1[:, 0].max() - m1[:, 0].min()
    assert width1 - width0 == pytest.approx(3.0, abs=1e-9), "each edge moved out by 1.5 nm"
    assert apply_misalignment(lay, draw_misalignment(MisalignmentParams(), lay, np.random.default_rng(0), severity=0.0)) == lay


def test_scenario_draw_is_seeded_and_round_trips_through_a_dict():
    lay = two_dot_layout("MOS")
    params = ScenarioParams.for_material("MOS")
    a1, a2 = (draw_scenario(lay, params, seed=7) for _ in range(2))
    c = draw_scenario(lay, params, seed=8)
    assert a1.to_dict() == a2.to_dict() and a1.to_dict() != c.to_dict()
    back = ScenarioDraw.from_dict(a1.to_dict())
    np.testing.assert_array_equal(back.charges.x_nm, a1.charges.x_nm)
    assert back.misalignment.shift_nm == a1.misalignment.shift_nm
    _, basis_a = apply_scenario(lay, a1)
    _, basis_back = apply_scenario(lay, back)
    assert basis_a.offset_ev.tobytes() == basis_back.offset_ev.tobytes()


def test_disorder_moves_the_dots_but_keeps_them_in_the_channel():
    lay = two_dot_layout()
    clean = derive_double_dot(lay, OP)
    moved = []
    for seed in range(5):
        draw = draw_scenario(lay, ScenarioParams.for_material("SiGe"), seed)
        perturbed, basis = apply_scenario(lay, draw)
        d = derive_double_dot(perturbed, OP, basis=basis)
        moved.append(np.hypot(d.dots[0].x_nm - clean.dots[0].x_nm, d.dots[0].y_nm - clean.dots[0].y_nm))
        assert all(lay.channel_nm[2] <= dot.y_nm <= lay.channel_nm[3] for dot in d.dots)
    assert max(moved) > 1.0, "at least one draw moves a dot by more than a grid pitch"
