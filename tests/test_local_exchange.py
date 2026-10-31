"""P0 for the exchange of the dots' valence electrons (#60, ADR 0031): the Hubbard model's J at its
symmetric point, the second-order limit away from the anticrossings, J following B_M as t² and
saturating where t reaches the charge gaps, J/h as the exchange frequency, no division by zero where t
vanishes, and the Pauli readout at the gated t. Analytic limits and invariants only."""

import warnings

import numpy as np
import pytest

from phantom.charge.ci import ConstantInteractionParams
from phantom.charge.dot import Barrier
from phantom.charge.hubbard import HubbardParams
from phantom.charge.hybrid import InterdotCoupling
from phantom.constants import HBAR_EV_S
from phantom.contract import reference_device
from phantom.device import Device
from phantom.spin.exchange import exchange_j_ev, local_exchange_ev, symmetric_point_v

T_EV = 20e-6


@pytest.fixture(scope="module")
def device() -> Device:
    return reference_device()


def _pair(t_ev: float) -> tuple[ConstantInteractionParams, HubbardParams]:
    """The default CI model with a constant coupling t, and the Hubbard model built from it."""
    ci = ConstantInteractionParams()
    coupling = InterdotCoupling(v_ref_v=(0.0,) * 3, t_ref_ev=t_ev, t_max_ev=1e-2, log_slopes_per_v=(0.0,) * 3)
    coupled = ConstantInteractionParams(ci.gate_names, ci.c_dot_af, ci.c_gate_af, ci.temperature_k, ci.window, coupling)
    hub = HubbardParams.from_ci(ci, Barrier("B", t0_ev=float(coupling.tunnel_coupling_ev(np.zeros(3))), v_ref_v=0.0, v_scale_v=1e9))
    return coupled, hub


def _gaps(ci: ConstantInteractionParams, v: np.ndarray, n: np.ndarray) -> tuple[float, float]:
    """(Δ_R, Δ_L): (N₁−1, N₂+1) and (N₁+1, N₂−1) above (N₁, N₂), from the CI energies."""
    u = ci.state_energies_ev(v, np.array([n, n + [-1, 1], n + [1, -1]]))
    return float(u[1] - u[0]), float(u[2] - u[0])


def test_it_is_the_hubbard_exchange_at_the_symmetric_point():
    """With one electron on each dot, the valence sector is the Hubbard model's two-electron sector."""
    ci, hub = _pair(T_EV)
    v = symmetric_point_v(hub)
    assert np.array_equal(ci.ground_state_occupation(v[:3]), [1, 1])
    assert float(local_exchange_ev(ci, v[:3])) == pytest.approx(float(exchange_j_ev(hub, v)), rel=1e-12)


def test_far_from_the_anticrossings_it_is_second_order_in_t():
    t_small = 5e-6  # the default CI model's gaps are about a meV
    ci, hub = _pair(t_small)
    centre = symmetric_point_v(hub)[:3]
    for shift in (-0.002, 0.0, 0.002):
        v = centre + np.array([shift, -shift, 0.0])
        n = ci.ground_state_occupation(v)
        d_r, d_l = _gaps(ci, v, n)
        assert min(d_r, d_l) > 100 * t_small
        assert float(local_exchange_ev(ci, v)) == pytest.approx(2 * t_small**2 / d_r + 2 * t_small**2 / d_l, rel=1e-3)


def test_on_the_device_it_follows_b_m_as_t_squared_and_saturates(device: Device):
    """Along B_M in one charge cell of the reference device, J tracks 2t²(1/Δ_R + 1/Δ_L) while t ≪ Δ and
    rises with t; once t passes a twentieth of the gaps the second-order form overestimates it."""
    v0 = np.array(device.programmed_point_v)
    jb = device.gate_names.index("B_M")
    db = np.linspace(-0.06, 0.06, 25)
    v = np.broadcast_to(v0, db.shape + v0.shape).copy()
    v[:, jb] += db
    s = device.at(v)
    ci = device.runtime_ci()
    cell = np.all(s.n_ground == s.n_ground[len(db) // 2], axis=-1)
    second = np.array([2 * t**2 * (1 / dr + 1 / dl) for t, (dr, dl) in zip(s.tunnel_coupling_ev, (_gaps(ci, x, n) for x, n in zip(s.v_carried_v, s.n_ground, strict=True)), strict=True)])
    gaps = np.array([min(_gaps(ci, x, n)) for x, n in zip(s.v_carried_v, s.n_ground, strict=True)])
    weak = cell & (s.tunnel_coupling_ev < 1e-2 * gaps)
    strong = cell & (s.tunnel_coupling_ev > 0.05 * gaps)
    assert weak.sum() >= 5 and strong.any()
    np.testing.assert_allclose(s.exchange_ev[weak], second[weak], rtol=1e-3)
    assert np.all(s.exchange_ev[strong] < second[strong]), "saturating below the second-order form"
    assert np.all(np.diff(s.exchange_ev[cell]) > 0), "opening B_M raises J"
    np.testing.assert_array_equal(s.exchange_frequency_hz, s.exchange_ev / (2 * np.pi * HBAR_EV_S))


def test_no_division_by_zero_where_t_vanishes(device: Device):
    v0 = np.array(device.programmed_point_v)
    jb = device.gate_names.index("B_M")
    v = np.broadcast_to(v0, (41,) + v0.shape).copy()
    v[:, jb] += np.linspace(-0.4, 0.0, 41)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        s = device.at(v)
    assert np.all(np.isfinite(s.exchange_ev)) and np.all(s.exchange_ev >= 0)
    ci, _ = _pair(1e-300)
    v = symmetric_point_v(_pair(T_EV)[1])[:3] + np.array([0.02, -0.02, 0.0])  # without t, J is the singlet's charge gap
    n = ci.ground_state_occupation(v)
    d_r, d_l = _gaps(ci, v, n)
    assert float(local_exchange_ev(ci, v)) == pytest.approx(-min(0.0, d_r, d_l), abs=1e-15)


def test_the_pauli_readout_reads_at_the_gated_coupling(device: Device, monkeypatch):
    from phantom.device import session

    seen = []
    original = session.psb_occupation
    monkeypatch.setattr(session, "psb_occupation", lambda eps, e_v, t, kt: seen.append(t) or original(eps, e_v, t, kt))
    v = np.array(device.programmed_point_v)
    v[device.gate_names.index("B_M")] -= 0.05
    run = device.session(noise=False)
    run.hold(v, 1e-6)
    run.spin_readout(0.5, 50e-6, shots=10, t_integration_s=1e-6, t1_read_s=1e-3)
    expected = float(device.runtime_ci().tunnel_coupling_ev(device.effective_voltages(v)))
    assert seen == [pytest.approx(expected, rel=1e-15)]
    assert expected < 1e-2 * device.tunnel_coupling_ev, "closing B_M closes the readout's coupling"
