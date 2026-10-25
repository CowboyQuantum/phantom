"""P0: analytic limits and invariants of the constant-interaction double dot."""

import numpy as np
import pytest

from phantom.charge import E_UNIT_V_PER_E_AF, ConstantInteractionParams


def cold(**kwargs) -> ConstantInteractionParams:
    return ConstantInteractionParams(temperature_k=0.01, **kwargs)


def test_occupation_steps_along_a_plunger_at_the_capacitive_period():
    m = cold()
    v = np.zeros((2000, 3))
    v[:, 0] = np.linspace(0.0, 0.3, 2000)
    n = m.mean_occupation(v)
    assert np.all(np.diff(n[:, 0]) >= -1e-9), "N1 never decreases as P1 rises"
    n_added = round(n[-1, 0] - n[0, 0])
    assert abs(n_added - round(0.3 / m.transition_period_v(0, "P1"))) <= 1
    away = np.abs(n[:, 0] - np.rint(n[:, 0]))
    assert np.mean(away < 0.02) > 0.9, "cold: occupation is integer almost everywhere"


def test_cross_capacitance_moves_dot_2_with_p1():
    m = cold()
    v = np.zeros((500, 3))
    v[:, 0] = np.linspace(0.0, 1.0, 500)
    n = m.mean_occupation(v)
    assert n[-1, 1] > n[0, 1] and n[-1, 1] < n[-1, 0]


def test_temperature_broadens_transitions():
    v = np.zeros((2001, 3))
    v[:, 0] = np.linspace(0.0, 0.02, 2001)
    sharp = lambda m: np.mean(np.abs(m.mean_occupation(v)[:, 0] % 1.0 - 0.5) > 0.4)
    assert sharp(ConstantInteractionParams(temperature_k=2.0)) < sharp(cold())


def test_lever_arms_are_the_derivative_of_the_chemical_potential():
    """μ(1,0 ← 0,0) = U(1,0) − U(0,0) shifts by −(C⁻¹ C_gate)_0j eV per volt of gate j."""
    m = ConstantInteractionParams()
    v0 = m.gate_voltages_for_ng(np.array([0.4, 0.3]))
    v0 = np.array([v0[0], v0[1], 0.0])

    def mu(v):
        states, u = m.energies(m.induced_charge(v))
        i10 = np.flatnonzero((states == [1, 0]).all(axis=1))[0]
        i00 = np.flatnonzero((states == [0, 0]).all(axis=1))[0]
        return u[i10] - u[i00]

    h = 1e-6
    for j in range(3):
        e = np.zeros(3)
        e[j] = h
        slope = (mu(v0 + e) - mu(v0 - e)) / (2 * h)
        assert slope == pytest.approx(-m.lever_arms_ev_per_v()[0, j], rel=1e-6, abs=1e-12)


def test_addition_energy_equals_period_times_lever_arm_without_cross_capacitance():
    m = ConstantInteractionParams(c_gate_af=np.array([[12.0, 0.0, 0.0], [0.0, 11.5, 0.0]]))
    e_add = m.transition_period_v(0, "P1") * m.lever_arms_ev_per_v()[0, 0]
    assert e_add == pytest.approx(m.addition_energies_ev()[0], rel=1e-12)
    assert m.addition_energies_ev()[0] == pytest.approx(
        E_UNIT_V_PER_E_AF * np.linalg.inv(m.c_dot_af)[0, 0], rel=1e-12
    )


def test_triple_points_are_degenerate_and_bound_the_right_regions():
    m = cold()
    tp1, tp2 = m.triple_points_ng()

    def u_of(n_g, state):
        states, u = m.energies(np.asarray(n_g))
        return u[np.flatnonzero((states == state).all(axis=1))[0]]

    for n_g, trio in ((tp1, ([0, 0], [1, 0], [0, 1])), (tp2, ([1, 0], [0, 1], [1, 1]))):
        energies = [u_of(n_g, s) for s in trio]
        assert np.ptp(energies) < 1e-15, f"three states degenerate at {n_g}"
    eps = 0.05
    for direction, expected in (((-1, -1), [0, 0]), ((1, -1), [1, 0]), ((-1, 1), [0, 1])):
        v = m.gate_voltages_for_ng(tp1 + eps * np.array(direction))
        v = np.array([v[0], v[1], 0.0])
        assert list(m.ground_state_occupation(v)) == expected
    v = m.gate_voltages_for_ng(tp2 + eps * np.array((1, 1)))
    assert list(m.ground_state_occupation(np.array([v[0], v[1], 0.0]))) == [1, 1]


def test_triple_point_separation_vanishes_without_mutual_capacitance():
    m = ConstantInteractionParams(c_dot_af=np.array([[85.0, 0.0], [0.0, 90.0]]))
    tp1, tp2 = m.triple_points_ng()
    np.testing.assert_allclose(tp1, tp2, atol=1e-15)
    assert m.mutual_charging_energy_ev() == 0.0
    tp1, tp2 = ConstantInteractionParams().triple_points_ng()
    assert np.linalg.norm(tp2 - tp1) > 0.1


def test_charge_lattice_and_gate_voltages_are_inverse():
    m = ConstantInteractionParams()
    n_g = np.array([1.3, 0.7])
    v = m.gate_voltages_for_ng(n_g)
    np.testing.assert_allclose(m.charge_lattice() @ v, n_g, rtol=1e-12)
    np.testing.assert_allclose(m.induced_charge(np.array([v[0], v[1], 0.0])), n_g, rtol=1e-12)


def test_sample_occupation_matches_mean_in_expectation():
    m = ConstantInteractionParams(temperature_k=1.0)
    v = np.array([0.5 * m.transition_period_v(0, "P1"), 0.0, 0.0])
    rng = np.random.default_rng(1)
    draws = np.stack([m.sample_occupation(v, rng) for _ in range(2000)])
    assert set(np.unique(draws[:, 0])) <= {0, 1}
    assert abs(draws[:, 0].mean() - m.mean_occupation(v)[0]) < 0.05


def test_shapes_are_validated():
    with pytest.raises(ValueError):
        ConstantInteractionParams(c_gate_af=np.zeros((2, 2)))
