"""P0: the two-site Hubbard model reduces to the CI model at t → 0 and anticrosses at t > 0."""

import numpy as np
import pytest

from phantom.charge import E_UNIT_V_PER_E_AF, Barrier, ConstantInteractionParams, Dot, HubbardParams
from phantom.charge.hubbard import HOP, N1, N2


def test_hopping_matrix_is_hermitian_and_conserves_number_and_spin():
    assert np.array_equal(HOP, HOP.T)
    total = N1 + N2
    rows, cols = np.nonzero(HOP)
    assert np.all(total[rows] == total[cols])
    assert HOP.sum() == 0.0 and np.count_nonzero(HOP) == 16  # 8 hops, each with its conjugate


def test_t_to_zero_reproduces_the_ci_ground_state_and_energy():
    # cold: the spinful Hubbard model has spin degeneracies the CI model does not count, so thermal
    # averages agree only as T -> 0; ground states and energies agree exactly
    ci = ConstantInteractionParams(temperature_k=0.005)
    h = HubbardParams.from_ci(ci, Barrier("B", t0_ev=0.0))
    a = np.linalg.inv(ci.c_dot_af)
    rng = np.random.default_rng(0)
    for _ in range(40):
        n_g = rng.uniform(0.0, 1.8, 2)
        v2 = ci.gate_voltages_for_ng(n_g)
        v_ci = np.array([v2[0], v2[1], 0.0])
        v_h = np.array([v2[0], v2[1], 0.0, 0.0])
        np.testing.assert_array_equal(h.ground_state_occupation(v_h), ci.ground_state_occupation(v_ci))
        np.testing.assert_allclose(h.mean_occupation(v_h), ci.mean_occupation(v_ci), atol=1e-6)
        e0 = h.spectrum(v_h)[0][0] + 0.5 * E_UNIT_V_PER_E_AF * n_g @ a @ n_g
        _, u = ci.energies(n_g)
        assert e0 == pytest.approx(u.min(), abs=1e-12)


def test_tunnel_coupling_anticrosses_the_interdot_transition():
    ci = ConstantInteractionParams(temperature_k=0.01)
    tp1, tp2 = ci.triple_points_ng()
    mid = (tp1 + tp2) / 2  # on the (1,0)–(0,1) interdot line
    v2 = ci.gate_voltages_for_ng(mid)
    widths = []
    for t0 in (0.0, 50e-6, 200e-6):
        h = HubbardParams.from_ci(ci, Barrier("B", t0_ev=t0, v_ref_v=0.0))
        v = np.array([v2[0], v2[1], 0.0, 0.0])
        n = h.ground_state_occupation(v)
        assert n.sum() == pytest.approx(1.0, abs=1e-9)
        if t0 > 0:  # at t = 0 the midpoint is exactly degenerate and the split is undefined
            assert n[0] == pytest.approx(0.5, abs=0.02), "charge is shared at the midpoint"
        # width of the interdot step: move along the detuning direction until N1 reaches 0.9
        direction = (tp2 - tp1) / np.linalg.norm(tp2 - tp1)
        perp = np.array([direction[1], -direction[0]])
        ds = np.linspace(0, 0.5, 2001)
        pts = mid[None, :] + ds[:, None] * perp[None, :]
        vv = np.stack([ci.gate_voltages_for_ng(p) for p in pts])
        vv = np.concatenate([vv, np.zeros((len(ds), 2))], axis=1)
        n1 = h.ground_state_occupation(vv)[:, 0]
        if n1[-1] < 0.5:
            n1 = 1 - n1
        widths.append(ds[np.argmax(n1 > 0.9)])
    assert widths[0] < widths[1] < widths[2], "the step broadens with t"


def test_barrier_is_exponential_in_its_gate():
    b = Barrier("B", t0_ev=20e-6, v_ref_v=-0.1, v_scale_v=0.05)
    assert b.tunnel_coupling_ev(-0.1) == pytest.approx(20e-6)
    assert b.tunnel_coupling_ev(-0.1 + 0.05 * np.log(2)) == pytest.approx(40e-6, rel=1e-12)
    with pytest.raises(ValueError):
        Barrier("B", v_scale_v=0.0)


def test_dot_charging_energy_is_a_few_mev_in_silicon():
    e_c = Dot("d1", radius_nm=25.0, epsilon_r=11.7).charging_energy_ev
    assert 3e-3 < e_c < 15e-3
    assert Dot("big", radius_nm=50.0).charging_energy_ev == pytest.approx(e_c / 2, rel=1e-12)
