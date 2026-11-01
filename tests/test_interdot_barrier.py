"""P0 for the interdot barrier's one transmission (#70, ADR 0034): its top follows every gate through the
screened slopes and the shape grid, continued linearly past the grid's end; Kemble's transmission of it
is shared by t, the interdot rate and the conductance, bit for bit; a barrier gate opened far past the
merge keeps the dots one at every occupancy; the screening factors are compile parameters, and at 1 the
bare potential; the reference device's honeycomb-to-superdot transition spans tens of millivolts of B_M
and the plungers act several times more weakly. Analytic limits and invariants only."""

import dataclasses
import math

import numpy as np
import pytest

from phantom.api import CompileParams, compile_device, two_dot_layout
from phantom.charge.barrier import InterdotBarrier, kemble_transmission
from phantom.charge.shape import ShapeGrid, row_length
from phantom.contract import reference_device
from phantom.device import Device


@pytest.fixture(scope="module")
def device() -> Device:
    return reference_device()


def _along(device: Device, gate: str, dv: np.ndarray) -> np.ndarray:
    v0 = device.operating_point()
    return v0 + np.asarray(dv)[:, None] * np.eye(len(v0))[device.gate_names.index(gate)]


def _programmed(device: Device, v_carried: np.ndarray) -> np.ndarray:
    """Programmed voltages that carry `v_carried` (the wiring draw is affine: invert it by its Jacobian)."""
    v0 = np.array(device.programmed_point_v)
    base = device.effective_voltages(v0)
    jac = np.stack([device.effective_voltages(v0 + e) - base for e in np.eye(len(v0))], axis=1)
    return v0 + np.linalg.solve(jac, v_carried - base)


def _crossing(x: np.ndarray, y: np.ndarray, level: float) -> float:
    k = int(np.argmax(y >= level))
    assert k > 0, "the level is crossed inside the range"
    return float(np.interp(level, y[k - 1 : k + 1], x[k - 1 : k + 1]))


def _top_only_grid(n_gates: int, j: int, nodes: tuple[float, ...], tops: tuple[float, ...]) -> ShapeGrid:
    """A grid along gate `j` whose only non-zero column is the barrier top's."""
    rows = []
    for top in tops:
        row = [0.0] * row_length(n_gates)
        row[-1] = top
        rows.append(tuple(row))
    return ShapeGrid((j,), n_gates, (nodes,), (tuple(rows),), (False,), (True,))


def test_the_top_is_linear_in_the_gates_with_the_screened_slopes_and_the_transmission_is_kemble_s():
    b = InterdotBarrier(v_ref_v=(0.1, 0.2, 0.3), top_ev=4e-3, top_slopes_ev_per_v=(-0.15, -0.12, -0.02), hbar_omega_b_ev=4e-3, screening=(0.4, 0.12, 0.12))
    v = np.array([[0.1, 0.2, 0.3], [0.2, 0.2, 0.3], [0.1, 0.3, 0.3], [0.1, 0.2, 1.3]])
    expected = 4e-3 + np.array([0.0, -0.15 * 0.4 * 0.1, -0.12 * 0.12 * 0.1, -0.02 * 0.12 * 1.0])
    np.testing.assert_allclose(b.barrier_top_ev(v), expected, rtol=1e-14)
    np.testing.assert_allclose(b.transmission(v), kemble_transmission(expected, 4e-3), rtol=1e-12)
    at_top = np.array([0.1 + 4e-3 / (0.15 * 0.4), 0.2, 0.3])
    assert float(b.transmission(at_top)) == pytest.approx(0.5, rel=1e-12), "½ when the level reaches the top"
    np.testing.assert_allclose(b.log_odds_change(v), -2 * np.pi * (expected - 4e-3) / 4e-3, rtol=1e-12, atol=1e-15)
    bare = dataclasses.replace(b, screening=None)
    np.testing.assert_allclose(bare.barrier_top_ev(v), 4e-3 + np.array([0.0, -0.015, -0.012, -0.02]), rtol=1e-14), "no screening: the bare slopes"
    with pytest.raises(ValueError):
        InterdotBarrier(v_ref_v=(0.0,) * 3, top_ev=0.0, top_slopes_ev_per_v=(0.0,) * 2, hbar_omega_b_ev=4e-3)


def test_along_a_shape_gate_the_top_follows_the_grid_and_continues_linearly_past_its_ends():
    """Inside the grid the top is the grid's (PCHIP through the nodes, times the gate's screening);
    beyond an end it continues with the interpolant's end slope, where `change` holds."""
    grid = _top_only_grid(2, 0, (0.0, 0.1, 0.2, 0.3), (0.03, 0.02, 0.01, 0.0))  # a top falling 0.1 eV/V
    b = InterdotBarrier(v_ref_v=(0.3, 0.0), top_ev=5e-3, top_slopes_ev_per_v=(-0.5, -0.01), hbar_omega_b_ev=4e-3, screening=(0.5, 1.0), shape=grid)
    v = np.array([[0.3, 0.0], [0.15, 0.0], [0.0, 0.0], [-0.2, 0.0], [0.5, 0.0], [0.5, 0.1]])
    top = b.barrier_top_ev(v)
    assert top[0] == 5e-3, "at the reference point: the compiled top"
    assert top[1] == pytest.approx(5e-3 + 0.5 * 0.015, rel=1e-12), "between nodes: the grid, linear here"
    assert top[2] == pytest.approx(5e-3 + 0.5 * 0.03, rel=1e-12), "at the end node"
    assert top[3] == pytest.approx(5e-3 + 0.5 * 0.05, rel=1e-12), "beyond the low end: continued at −0.1 eV/V"
    assert top[4] == pytest.approx(5e-3 - 0.5 * 0.02, rel=1e-12), "beyond the merged end: the top keeps sinking"
    assert top[5] == pytest.approx(5e-3 - 0.5 * 0.02 - 0.01 * 0.1, rel=1e-12), "another gate adds through its slope"
    held = grid.quantity(grid.change(v), "barrier_top_ev")[..., 0]
    assert held[4] == pytest.approx(0.0, abs=1e-15) and held[3] == pytest.approx(0.03, rel=1e-12), "`change` holds the end nodes' values"
    np.testing.assert_allclose(grid.quantity_along(v, "barrier_top_ev", beyond="hold")[..., 0], held, rtol=1e-14)
    assert grid.merged(v).tolist() == [False, False, False, False, True, True]
    with pytest.raises(ValueError):
        grid.quantity_along(v, "barrier_top_ev", beyond="extrapolate")


def test_t_the_interdot_rate_and_the_conductance_share_the_transmission(device: Device):
    """At every gate point t = t_max·√T, Γ_m = f_a·T and g = T, with the compiled record at the operating
    point: `tunnel_coupling_ev`, the interdot barrier's `rate_hz` and slopes, and the e-fold scale of t."""
    b = device.runtime_interdot_barrier()
    assert b is not None and b.screening == device.interdot_barrier_screening
    rng = np.random.default_rng(0)
    v = device.operating_point() + rng.uniform(-0.2, 0.3, size=(50, len(device.gate_names)))
    s = device.at(np.stack([_programmed(device, x) for x in v]))
    transmission = b.transmission(s.v_carried_v)
    c = device.runtime_coupling()
    np.testing.assert_allclose(s.tunnel_coupling_ev, c.t_max_ev * np.sqrt(transmission), rtol=1e-12)
    inter = device.barriers[2]
    np.testing.assert_allclose(s.gamma_interdot_hz, inter.attempt_rate_hz * transmission, rtol=1e-9)
    np.testing.assert_array_equal(s.interdot_conductance, transmission)
    np.testing.assert_allclose(s.interdot_barrier_top_ev, b.barrier_top_ev(s.v_carried_v), rtol=1e-12)
    t0 = float(kemble_transmission(device.interdot_barrier_top_ev, device.interdot_barrier_hbar_omega_ev))
    assert device.tunnel_coupling_ev == pytest.approx(c.t_max_ev * math.sqrt(t0), rel=1e-14)
    assert inter.rate_hz == pytest.approx(inter.attempt_rate_hz * t0, rel=1e-14)
    expected_slopes = -2 * np.pi * np.array(device.interdot_barrier_screening) * np.array(device.interdot_barrier_top_slopes_ev_per_v) / device.interdot_barrier_hbar_omega_ev
    np.testing.assert_allclose(inter.log_rate_slopes_per_v, expected_slopes, rtol=1e-12)
    jb = device.gate_names.index("B_M")
    assert device.barrier_v_scale_v == pytest.approx(2.0 / expected_slopes[jb], rel=1e-12)
    t_down = device.runtime_ci().tunnel_coupling_ev(_along(device, "B_M", [-0.1 - device.barrier_v_scale_v, -0.1]))
    assert t_down[1] / t_down[0] == pytest.approx(math.e, rel=0.03), "the e-fold scale is t's while the barrier is closed"
    assert float(c.tunnel_coupling_ev(device.operating_point())) == pytest.approx(device.tunnel_coupling_ev, rel=1e-14)


def test_the_transition_spans_tens_of_millivolts_of_b_m_and_the_plungers_act_several_times_more_weakly(device: Device):
    """The owner's requirement (#70): along B_M the splitting f goes from 0.2 to 0.9 over 50–100 mV, and
    along the plunger diagonal it takes at least 1.8 times as much, each plunger's effective slope on
    the barrier top being under a third of B_M's."""
    v0 = np.array(device.programmed_point_v)
    jb, i1, i2 = (device.gate_names.index(g) for g in ("B_M", "P1", "P2"))
    db = np.linspace(-0.1, 0.3, 801)
    v = np.broadcast_to(v0, db.shape + v0.shape).copy()
    v[:, jb] += db
    f = device.at(v).interdot_splitting
    width_b = _crossing(db, f, 0.9) - _crossing(db, f, 0.2)
    assert 0.05 <= width_b <= 0.1, width_b
    dp = np.linspace(-0.15, 0.4, 1101)
    v = np.broadcast_to(v0, dp.shape + v0.shape).copy()
    v[:, i1] += dp
    v[:, i2] += dp
    f = device.at(v).interdot_splitting
    width_p = _crossing(dp, f, 0.9) - _crossing(dp, f, 0.2)
    assert width_p > 1.8 * width_b, (width_p, width_b)
    b = device.runtime_interdot_barrier()
    eff = np.array(b.screening) * np.array(b.top_slopes_ev_per_v)
    assert max(abs(eff[i1]), abs(eff[i2])) < abs(eff[jb]) / 3


def test_a_barrier_gate_opened_far_past_the_merge_keeps_the_dots_one_at_every_occupancy(device: Device):
    """With B_M at 5 V the top keeps sinking past the shape grid's end, so lowering both plungers by
    0.15 V (a few electrons each) leaves one merged dot: g stays 1 and f at its ceiling, while the shape
    is still held and flagged. Below the compiled plungers, opening B_M raises t, g and f (#69)."""
    v0 = np.array(device.programmed_point_v)
    jb, i1, i2 = (device.gate_names.index(g) for g in ("B_M", "P1", "P2"))
    v = np.broadcast_to(v0, (4,) + v0.shape).copy()
    v[:, jb] = 5.0
    v[:, i1] += [0.0, -0.075, -0.1, -0.15]
    v[:, i2] += [0.0, -0.075, -0.1, -0.15]
    s = device.at(v)
    assert np.all(s.dots_merged) and np.all(s.interdot_conductance > 1 - 1e-12) and np.all(s.interdot_splitting > 0.97)
    assert np.all(s.interdot_barrier_top_ev < -0.1), "the top is far below the level, not held"
    assert s.n_mean[0, 0] - s.n_mean[-1, 0] > 2, "the plungers still remove electrons"
    np.testing.assert_allclose(s.dot_positions_nm, np.broadcast_to(s.dot_positions_nm[0], (4, 2, 2)), atol=1e-9), "the shape is held"
    low = v0.copy()
    low[i1], low[i2] = 0.428, 0.449
    opening = np.broadcast_to(low, (4,) + v0.shape).copy()
    opening[:, jb] = [0.2, 0.25, 0.3, 0.4]
    s = device.at(opening)
    assert np.all(np.diff(s.tunnel_coupling_ev) > 0) and s.tunnel_coupling_ev[-1] > 100 * s.tunnel_coupling_ev[0]
    assert s.interdot_splitting[0] < 0.25 and s.interdot_splitting[-1] > 0.8


def test_the_screening_factors_are_compile_parameters_and_one_is_the_bare_potential(device: Device):
    """Compiled at screening 1 and 1, the interdot barrier's slopes are the bare top's slopes times
    −2π/ħω_b, the top itself, t and the shape grid are the same device's: the screening scales how the
    gates move the top, not where it is. In the bare potential a plunger moves the top nearly as much as
    B_M does (#48)."""
    layout = two_dot_layout(material="SiGe", pitch_nm=80.0)
    op = dict(zip(device.gate_names, device.programmed_point_v, strict=True))
    bare = compile_device(layout, CompileParams(operating_point_v=op, interdot_barrier_screening=1.0, interdot_cross_screening=1.0), seed=device.seed)
    assert bare.interdot_barrier_screening == (1.0,) * len(device.gate_names)
    assert bare.interdot_barrier_top_ev == device.interdot_barrier_top_ev and bare.interdot_barrier_hbar_omega_ev == device.interdot_barrier_hbar_omega_ev
    assert bare.interdot_barrier_top_slopes_ev_per_v == device.interdot_barrier_top_slopes_ev_per_v
    assert bare.shape.changes == device.shape.changes
    assert bare.tunnel_coupling_ev == device.tunnel_coupling_ev and bare.barriers[2].rate_hz == device.barriers[2].rate_hz
    np.testing.assert_allclose(bare.barriers[2].log_rate_slopes_per_v, -2 * np.pi * np.array(bare.interdot_barrier_top_slopes_ev_per_v) / bare.interdot_barrier_hbar_omega_ev, rtol=1e-12)
    jb, i1 = (device.gate_names.index(g) for g in ("B_M", "P1"))
    assert device.interdot_barrier_screening[jb] == pytest.approx(0.4) and device.interdot_barrier_screening[i1] == pytest.approx(0.12)
    b = bare.runtime_interdot_barrier()

    def drop(gate: str) -> float:
        return float(b.barrier_top_ev(_along(bare, gate, [0.0])[0]) - b.barrier_top_ev(_along(bare, gate, [0.1])[0]))

    assert drop("P1") > 0.5 * drop("B_M") > 0
