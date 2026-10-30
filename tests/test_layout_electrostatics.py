"""P0 for Step 2: the closed-form potential, the derived dots and the CI parameters built from a layout."""

import numpy as np
import pytest

from phantom.charge import ConstantInteractionParams, derive_double_dot
from phantom.electrostatics import GateBasis, polygon_solid_angle, rectangle_potential_dls
from phantom.electrostatics.dots import harmonic_radius_nm
from phantom.electrostatics.schrodinger import dot_levels_ev
from phantom.layout import Gate, Layout, MaterialClass, two_dot_layout

OPERATING_POINT = np.array([0.2, 0.6, 0.2, 0.6, 0.2, 0.6, 0.7, 0.7, 0.35, 0.35, 0.7, 0.7])  # B_L P1 B_M P2 B_R SG R_L R_R SB_L SB_R SR_L SR_R


def test_polygon_solid_angle_matches_the_rectangle_closed_form():
    x = np.linspace(-120, 120, 25)
    y = np.linspace(-150, 150, 31)
    rect = Gate.rectangle("g", -22, 22, -60, 60)
    for depth in (10.0, 44.0, 120.0):
        sa = polygon_solid_angle(rect.vertices, x[None, :], y[:, None], depth) / (2 * np.pi)
        closed = rectangle_potential_dls(x[None, :], y[:, None], depth, -22, 22, -60, 60)
        np.testing.assert_allclose(sa, closed, rtol=1e-12, atol=1e-14)
    # orientation and fan triangulation do not matter: reversed vertices, and a concave L-shape split two ways
    reversed_rect = Gate("r", tuple(reversed(rect.polygon_nm)))
    np.testing.assert_allclose(polygon_solid_angle(reversed_rect.vertices, 5.0, -7.0, 30.0), polygon_solid_angle(rect.vertices, 5.0, -7.0, 30.0))
    ell = np.array([(0, 0), (60, 0), (60, 20), (20, 20), (20, 60), (0, 60)], dtype=float)
    a = polygon_solid_angle(ell, 10.0, 10.0, 25.0)
    b = polygon_solid_angle(np.array([(0, 0), (60, 0), (60, 20), (0, 20)], float), 10.0, 10.0, 25.0) + polygon_solid_angle(np.array([(0, 20), (20, 20), (20, 60), (0, 60)], float), 10.0, 10.0, 25.0)
    assert a == pytest.approx(b, rel=1e-12)


def test_limits_of_the_potential():
    huge = Gate.rectangle("g", -1e6, 1e6, -1e6, 1e6)
    assert polygon_solid_angle(huge.vertices, 0.0, 0.0, 40.0) / (2 * np.pi) == pytest.approx(1.0, abs=1e-4)
    small = Gate.rectangle("g", -10, 10, -10, 10)
    far = polygon_solid_angle(small.vertices, 5000.0, 0.0, 40.0) / (2 * np.pi)
    assert 0 < far < 1e-5
    # a point above the centre of a square at depth d sees the known solid angle 4 arctan(a²/(d √(2a² + d²)))
    a, d = 30.0, 40.0
    square = Gate.rectangle("g", -a, a, -a, a)
    expected = 4 * np.arctan(a * a / (d * np.sqrt(2 * a * a + d * d)))
    assert polygon_solid_angle(square.vertices, 0.0, 0.0, d) == pytest.approx(expected, rel=1e-12)


def test_layout_round_trips_through_json(tmp_path):
    lay = two_dot_layout("MOS")
    path = tmp_path / "layout.json"
    lay.to_json(path)
    back = Layout.from_json(path)
    assert back == lay
    assert back.stack.effective_depth_nm == pytest.approx(8.0 * 11.7 / 3.9)


def test_lever_arms_are_the_derivative_of_the_dot_energy():
    lay = two_dot_layout()
    basis = GateBasis.from_layout(lay)
    d = derive_double_dot(lay, OPERATING_POINT, basis=basis)
    dot = d.dots[0]
    h = 1e-3
    for j in range(len(lay.gates)):
        e = np.zeros(len(lay.gates))
        e[j] = h
        u_plus = basis.energy_ev(OPERATING_POINT + e)
        u_minus = basis.energy_ev(OPERATING_POINT - e)
        ix = int(np.argmin(np.abs(basis.x_nm - dot.x_nm)))
        iy = int(np.argmin(np.abs(basis.y_nm - dot.y_nm)))
        slope = (u_plus[iy, ix] - u_minus[iy, ix]) / (2 * h)  # eV per V
        assert slope == pytest.approx(-dot.lever_arms[j], rel=1e-9, abs=1e-12)
    np.testing.assert_allclose(d.ci.lever_arms_ev_per_v(), np.stack([d.dots[0].lever_arms, d.dots[1].lever_arms]), rtol=1e-12)


def test_derived_double_dot_is_symmetric_and_plausible():
    d = derive_double_dot(two_dot_layout(), OPERATING_POINT)
    d1, d2 = d.dots
    assert d1.x_nm == pytest.approx(-d2.x_nm) and d1.y_nm == d2.y_nm
    assert 0.0 <= d1.y_nm <= 6.0, "the sensor's baths pull the dots a few nm towards the sensor channel"
    assert d1.lever_arms[1] == pytest.approx(d2.lever_arms[3])  # P1 on dot 1 equals P2 on dot 2
    assert d1.lever_arms[1] > d1.lever_arms[2] > d1.lever_arms[3], "own plunger > middle barrier > other plunger"
    assert 1e-3 < d.ci.addition_energies_ev()[0] < 30e-3, "meV-scale charging energy"
    assert d.interdot_barrier.height_ev > 0, "a barrier separates the dots at this operating point"
    assert d.sensor_dot is not None and d.sensor_alpha_v_per_v[5] > 0.1 and d.sensor_alpha_v_per_v[1] < 0.02
    assert isinstance(d.ci, ConstantInteractionParams) and d.ci.gate_names == two_dot_layout().gate_names


def test_layout_honeycomb_period_matches_the_derived_capacitance():
    d = derive_double_dot(two_dot_layout(), OPERATING_POINT, temperature_k=0.02)
    ci = d.ci
    j = ci.gate_names.index("P1")
    v = np.repeat(OPERATING_POINT[None, :], 4000, axis=0)
    v[:, j] = OPERATING_POINT[j] + np.linspace(0.0, 0.2, 4000)
    n1 = ci.mean_occupation(v)[:, 0]
    steps = np.flatnonzero(np.diff(np.rint(n1)) > 0.5)
    spacing = np.diff(v[steps, j])
    assert len(spacing) >= 2
    # along one plunger line the dot-1 transitions are shifted by every dot-2 transition crossed on the
    # way (the interdot coupling), so single spacings scatter; their mean is the capacitive period
    assert spacing.mean() == pytest.approx(ci.transition_period_v(0, "P1"), rel=0.05)
    shift_v = ci.mutual_charging_energy_ev() / ci.lever_arms_ev_per_v()[0, j]
    assert np.all(np.abs(spacing - ci.transition_period_v(0, "P1")) <= shift_v + 1e-4)


def test_barrier_opens_as_its_gate_rises():
    lay = two_dot_layout()
    basis = GateBasis.from_layout(lay)
    heights = []
    for v_bm in (0.0, 0.2, 0.4):
        v = OPERATING_POINT.copy()
        v[2] = v_bm
        heights.append(derive_double_dot(lay, v, basis=basis).interdot_barrier.height_ev)
    assert heights[0] > heights[1] > heights[2]


def test_schrodinger_levels_match_the_harmonic_oscillator():
    m = MaterialClass.sige().effective_mass_me
    k = 1e-4  # eV/nm²: U = ½ k r²
    x = np.arange(-60, 60.5, 1.0)
    u = 0.5 * k * (x[None, :] ** 2 + x[:, None] ** 2)
    levels = dot_levels_ev(u, 1.0, 1.0, m, n_levels=4)
    hbar_omega = np.sqrt(2 * 0.0380998212 * k / m)  # ħω = sqrt(ħ²k/m*) with ħ²/2m_e = 0.0381 eV nm²
    assert levels[1] - levels[0] == pytest.approx(hbar_omega, rel=0.02)
    assert levels[2] - levels[1] == pytest.approx(0.0, abs=0.02 * hbar_omega), "first excited level is doubly degenerate"
    assert harmonic_radius_nm(k, m) == pytest.approx((2 * 0.0380998212 / (m * k)) ** 0.25)
