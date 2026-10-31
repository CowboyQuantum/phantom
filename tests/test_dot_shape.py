"""P0 for the double dot's shape along the barrier gates (#59, ADR 0030): exactly the compiled device at
the operating point, the compile tier at the grid's nodes and close to it between them, continuous
capacitances and induced charges, the dots drawn together and coupled more strongly as B_M opens until
they merge, each lead barrier moving its own dot, and t and the rates following the moved dots. Analytic
limits and invariants only."""

import numpy as np
import pytest

from phantom.charge.derived import capacitance_matrix_af, derive_double_dot
from phantom.compile.device import SHAPE_JUMP_NM
from phantom.constants import E_UNIT_V_PER_E_AF
from phantom.contract import reference_device
from phantom.device import Device
from phantom.electrostatics.dots import refined_dot


@pytest.fixture(scope="module")
def device() -> Device:
    return reference_device()


def _along(device: Device, gate: str, dv: np.ndarray) -> np.ndarray:
    v0 = device.operating_point()
    return v0 + np.asarray(dv)[:, None] * np.eye(len(v0))[device.gate_names.index(gate)]


def _compile_tier(device: Device, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(refined positions (2, 2), capacitance matrix) re-derived from the potential at carried `v`."""
    lay, basis = device.basis()
    derived = derive_double_dot(lay, v, basis=basis)
    u = basis.energy_ev(v)
    d1, d2 = (refined_dot(basis, u, d) for d in derived.dots)
    return np.array([[d1.x_nm, d1.y_nm], [d2.x_nm, d2.y_nm]]), capacitance_matrix_af(d1, d2, lay.stack.epsilon_well)


def test_at_the_operating_point_it_is_the_compiled_device(device: Device):
    v0 = device.operating_point()
    assert np.all(device.shape.change(v0) == 0.0)
    c, n_g = device.runtime_ci().local(v0)
    np.testing.assert_array_equal(c, np.array(device.c_dot_af))
    np.testing.assert_allclose(n_g, v0 @ np.array(device.c_gate_af).T / E_UNIT_V_PER_E_AF, rtol=1e-12)
    s = device.at(v0)
    np.testing.assert_array_equal(s.dot_positions_nm, [[d.x_nm, d.y_nm] for d in device.dots])
    assert not s.dots_merged
    assert float(s.tunnel_coupling_ev) == pytest.approx(device.tunnel_coupling_ev, rel=1e-12)


@pytest.mark.parametrize("gate", ["B_L", "B_M", "B_R"])
def test_at_the_nodes_and_between_them_it_follows_the_compile_tier(device: Device, gate: str):
    """At a node the positions and the capacitances' changes are the compile tier's; halfway between
    nodes the interpolant stays within 0.5 nm and 2 % of a direct derivation. The compile tier takes a
    dot's curvature at a grid point, so its capacitances move in small steps as the minimum changes
    pixel (about 1 %); the interpolant runs smoothly through the nodes."""
    g = device.shape
    k = g.gate_indices.index(device.gate_names.index(gate))
    nodes = np.array(g.nodes_v[k])
    v0 = device.operating_point()
    j = device.gate_names.index(gate)
    ref_pos, ref_c = _compile_tier(device, v0)
    for x in list(nodes[::3]) + list(0.5 * (nodes[1:] + nodes[:-1]))[::3]:
        v = v0.copy()
        v[j] = x
        pos, c = _compile_tier(device, v)
        s = device.at(device_programmed(device, v))
        runtime_c = device.runtime_ci().local(v)[0]
        on_node = np.any(np.isclose(nodes, x, atol=1e-12))
        np.testing.assert_allclose(s.dot_positions_nm - [[d.x_nm, d.y_nm] for d in device.dots], pos - ref_pos, atol=1e-9 if on_node else 0.5)
        np.testing.assert_allclose(runtime_c - np.array(device.c_dot_af), c - ref_c, atol=1e-12 if on_node else 0.02 * np.abs(ref_c).max())


def device_programmed(device: Device, v_carried: np.ndarray) -> np.ndarray:
    """Programmed voltages that carry `v_carried` (the wiring draw is affine: invert it by its Jacobian)."""
    v0 = np.array(device.programmed_point_v)
    base = device.effective_voltages(v0)
    jac = np.stack([device.effective_voltages(v0 + e) - base for e in np.eye(len(v0))], axis=1)
    return v0 + np.linalg.solve(jac, v_carried - base)


def test_the_capacitances_and_induced_charge_are_continuous(device: Device):
    """Along each barrier gate in 1 mV steps the capacitances and the induced charge have no jumps: their
    second differences stay far below their first, at the nodes and where the potential's grid minimum
    changes pixel alike."""
    ci = device.runtime_ci()
    for gate in ("B_L", "B_M", "B_R"):
        j = device.gate_names.index(gate)
        nodes = device.shape.nodes_v[device.shape.gate_indices.index(j)]
        v = _along(device, gate, np.arange(nodes[0], nodes[-1], 1e-3) - device.operating_point()[j])  # inside the grid: held beyond it
        c, n_g = ci.local(v)
        # n_g = Cφ/E counts the charge absolutely (about 23 e), so a dot growing by 10 % takes about two
        # more electrons at the same potential: a steep slope, not a jump
        for x, step in ((c, 0.05), (n_g, 0.2)):
            first, second = np.abs(np.diff(x, axis=0)).max(), np.abs(np.diff(x, n=2, axis=0)).max()
            assert first < step and second < 0.3 * first, gate  # a jump would make them equal


def test_opening_b_m_draws_the_dots_together_until_they_merge(device: Device):
    g = device.shape
    k = g.gate_indices.index(device.gate_names.index("B_M"))
    assert g.merged_above[k], "the reference device's dots merge as B_M opens"
    top = g.nodes_v[k][-1] - device.operating_point()[device.gate_names.index("B_M")]
    v = _along(device, "B_M", np.linspace(-0.15, top, 31))
    states = device.at(np.stack([device_programmed(device, x) for x in v]))
    gap = np.hypot(*(states.dot_positions_nm[:, 1] - states.dot_positions_nm[:, 0]).T)
    cm = -device.runtime_ci().local(v)[0][:, 0, 1]
    assert np.all(np.diff(gap) <= 1e-9) and gap[-1] < gap[0] - 5, "the dots approach"
    assert np.all(np.diff(cm) >= -1e-9) and cm[-1] > 1.2 * cm[0], "and couple more strongly"
    assert np.all(np.diff(states.tunnel_coupling_ev) >= 0) and states.tunnel_coupling_ev[-1] > 1e3 * states.tunnel_coupling_ev[0], "t rises"
    assert not states.dots_merged.any()
    beyond = device.at(device_programmed(device, _along(device, "B_M", [top + 0.02])[0]))
    assert bool(beyond.dots_merged)
    np.testing.assert_allclose(beyond.dot_positions_nm, states.dot_positions_nm[-1], atol=1e-9), "held at the last node"
    assert SHAPE_JUMP_NM > 0


@pytest.mark.parametrize(("gate", "dot"), [("B_L", 0), ("B_R", 1)])
def test_a_lead_barrier_moves_its_own_dot(device: Device, gate: str, dot: int):
    """Opening a lead barrier pulls its dot towards the lead by several nm; the other dot moves less
    than a third of that."""
    v = _along(device, gate, [-0.1, 0.1])
    s = device.at(np.stack([device_programmed(device, x) for x in v]))
    moved = np.abs(s.dot_positions_nm[1] - s.dot_positions_nm[0])[:, 0]
    assert moved[dot] > 4.0
    assert moved[1 - dot] < moved[dot] / 3


def test_the_rates_follow_the_moved_dots(device: Device):
    """Along B_M the interdot rate's exponent is the grid's change of the action at the fixed level and
    t's is half of it, so below saturation Γ_m/t² stays fixed while Γ_m moves by decades."""
    v = _along(device, "B_M", np.linspace(-0.15, 0.0, 7))
    rates = device.runtime_rates().rates_hz(v)
    t = device.runtime_coupling().tunnel_coupling_ev(v)
    ratio = np.log(rates[:, 2]) - 2 * np.log(t)
    np.testing.assert_allclose(ratio, ratio[-1], atol=0.1), "up to the logistic's saturation"
    assert np.all(np.diff(rates[:, 2]) > 0) and rates[-1, 2] > 1e3 * rates[0, 2]
