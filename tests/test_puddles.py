"""P0 for the dots growing with their occupancy (#72, ADR 0035): the puddles' scale √(n/N_ref), the
capacitances scaled bit-exactly at the operating point, the classical N ∝ φ^(3/2), the barrier seen from
the puddles' edges (zero shift at the anchor, higher for smaller puddles, the saddle for larger ones),
and on the reference device the addition energy falling as N^(−1/3), the geometric splitting growing,
the barrier top rising toward (1,1) at fixed gates and the superdot front moving from high occupancy
toward (1,1) as B_M opens. Analytic limits and invariants only."""

import dataclasses
import math

import numpy as np
import pytest

from phantom.api import CompileParams, compile_device, two_dot_layout
from phantom.charge.puddle import PuddleGrowth
from phantom.contract import reference_device
from phantom.device import Device


@pytest.fixture(scope="module")
def device() -> Device:
    return reference_device()


def _puddle(n_ref=(20.0, 20.0), turning=(15.0, 15.0)) -> PuddleGrowth:
    """A symmetric cut: the bare potential −10 meV at each dot rising to +5 meV at the middle, 60 nm apart."""
    s = np.linspace(0.0, 60.0, 121)
    u = -0.010 + 0.015 * np.sin(np.pi * s / 60.0)  # zero (the level) at 15.0 nm and 45.0 nm
    return PuddleGrowth(n_ref, turning, tuple(s), tuple(u))


def test_the_scale_is_the_square_root_of_the_induced_charge_and_the_capacitances_follow_it():
    p = _puddle()
    n = np.array([[20.0, 20.0], [5.0, 20.0], [80.0, 45.0], [0.1, 20.0]])
    s = p.scale(n)
    np.testing.assert_allclose(s, np.sqrt(np.array([[1.0, 1.0], [0.25, 1.0], [4.0, 2.25], [0.5 / 20.0, 1.0]])), rtol=1e-15), "√(n/N_ref), floored at half an electron"
    c = np.array([[85.0, -20.0], [-20.0, 90.0]])
    scaled = p.capacitances_af(np.broadcast_to(c, (4, 2, 2)), s)
    np.testing.assert_array_equal(scaled[0], c), "exactly the compiled matrix at the operating point"
    assert scaled[1, 0, 0] == pytest.approx(65.0 * 0.5 + 20.0 * 0.5, rel=1e-12), "dot 1's self-capacitance halved, the mutual one with it"
    assert scaled[1, 1, 1] == pytest.approx(70.0 + 20.0 * 0.5, rel=1e-12), "dot 2's self-capacitance unchanged, the mutual one halved"
    assert -scaled[2, 0, 1] == pytest.approx(20.0 * 2.0 * 1.5, rel=1e-12), "C_m ∝ s₁s₂"
    np.testing.assert_allclose(p.radii_nm(s), s * 15.0, rtol=1e-15)
    with pytest.raises(ValueError):
        PuddleGrowth((20.0, 20.0), (15.0, 15.0), (0.0, 1.0), (0.0,))


def test_the_barrier_is_seen_from_the_puddles_edges():
    """Zero shift at the anchor, where each puddle fills the bare allowed region; a smaller puddle sees
    the bare profile from further in, so the top above its edge is higher; a larger one reaches toward
    the saddle and lowers it, by at most the barrier's height; the two dots' shifts are averaged."""
    p = _puddle()

    def u(r_nm: float) -> float:
        return -0.010 + 0.015 * math.sin(math.pi * r_nm / 60.0)

    assert p.saddle_nm == pytest.approx(30.0)
    assert float(p.barrier_shift_ev(np.array([1.0, 1.0]))) == 0.0
    smaller = float(p.barrier_shift_ev(np.array([0.5, 1.0])))  # dot 1's edge moves from 15 nm in to 7.5 nm, where the profile is lower
    assert smaller == pytest.approx(0.5 * (u(15.0) - u(7.5)), rel=1e-6) and smaller > 0
    larger = float(p.barrier_shift_ev(np.array([1.0, 1.5])))  # dot 2's edge moves out to 22.5 nm from it, toward the saddle
    assert larger == pytest.approx(0.5 * (u(15.0) - u(22.5)), rel=1e-6) and larger < 0
    assert float(p.barrier_shift_ev(np.array([10.0, 10.0]))) == pytest.approx(u(15.0) - u(30.0), rel=1e-6), "both at the saddle: the barrier is gone"
    s = np.stack([np.linspace(0.3, 3.0, 50)] * 2, axis=-1)
    assert np.all(np.diff(p.barrier_shift_ev(s)) <= 1e-15), "monotone in the puddles' size"


def test_at_the_operating_point_nothing_changes_and_the_induced_charge_is_the_classical_n_of_phi(device: Device):
    """s = 1 exactly at the operating point: the capacitances, t, the rates and the splitting are the
    compiled device's. Along a plunger the induced charge with the scaled capacitances is n·√(n/N_ref),
    the classical N ∝ φ^(3/2) of a puddle whose capacitance grows with its radius."""
    ci = device.runtime_ci()
    assert ci.puddle is not None and device.puddle is not None
    v0 = device.operating_point()
    np.testing.assert_array_equal(ci.occupancy_scale(v0), [1.0, 1.0])
    c, n_g = dataclasses.replace(ci, splitting=None).local(v0)
    np.testing.assert_array_equal(c, np.array(device.c_dot_af))
    np.testing.assert_allclose(n_g, device.puddle.n_ref, rtol=1e-15)
    assert float(ci.tunnel_coupling_ev(v0)) == device.tunnel_coupling_ev
    np.testing.assert_allclose(device.runtime_rates().rates_hz(v0), [b.rate_hz for b in device.barriers], rtol=1e-12)
    s = device.at(np.array(device.programmed_point_v))
    np.testing.assert_allclose(s.puddle_radii_nm, device.puddle.turning_nm, rtol=1e-3), "the programmed point sits a fraction of an electron from the operating point"
    np.testing.assert_array_equal(device.puddle.radii_nm(ci.occupancy_scale(v0)), device.puddle.turning_nm)
    fixed = dataclasses.replace(ci, puddle=None)
    j = device.gate_names.index("P1")
    v = v0 + np.linspace(-0.3, 0.2, 11)[:, None] * np.eye(len(v0))[j]
    n_fixed = fixed.induced_charge(v)
    n = ci.induced_charge(v)
    expected = n_fixed * np.sqrt(np.maximum(n_fixed, 0.5) / np.array(device.puddle.n_ref))
    # the mutual term scales as s₁s₂ rather than s_i, so the relation holds to the mutual capacitance's share
    np.testing.assert_allclose(n, expected, rtol=0.02)
    assert np.all(np.diff(n[:, 0]) > np.diff(n_fixed[:, 0]) * 1.2), "more electrons per volt than the fixed dot, as N ∝ φ^(3/2)"


def test_the_dots_grow_as_n_to_the_third_and_the_cells_shrink_toward_high_occupancy(device: Device):
    """Walking both plungers from (1,1) to beyond the compiled occupancy at the programmed B_M: the
    puddle radius grows as N^(1/3), the addition energy falls as N^(−1/3), the geometric splitting f_C
    grows, and the barrier top at fixed gates falls, so t rises by orders of magnitude."""
    v0 = np.array(device.programmed_point_v)
    i1, i2 = (device.gate_names.index(g) for g in ("P1", "P2"))
    dp = np.linspace(-0.9, 0.1, 11)
    v = np.broadcast_to(v0, dp.shape + v0.shape).copy()
    v[:, i1] += dp
    v[:, i2] += dp
    s = device.at(v)
    n = s.n_mean
    assert n[0, 0] < 3 and n[-1, 0] > 20, "from the first electrons to beyond the compiled occupancy"
    ratio = s.puddle_radii_nm / s.puddle_radii_nm[-1]
    induced = device.runtime_ci().induced_charge(s.v_carried_v)  # the puddles' own N, continuous
    np.testing.assert_allclose(ratio, np.cbrt(np.maximum(induced, 0.5) / induced[-1]), rtol=0.03), "R ∝ N^(1/3), from the induced charge"
    e_add = s.addition_energies_ev
    assert np.all(np.diff(e_add, axis=0) < 0), "the addition energy falls with occupancy"
    assert e_add[0, 0] / e_add[-1, 0] > 2.0, "by more than a factor 2 from the first electrons to about 25"
    geometric = dataclasses.replace(device.runtime_ci(), splitting=None)
    f_c = geometric.interdot_splitting(s.v_carried_v)[1]
    assert np.all(np.diff(f_c) > 0) and f_c[0] < 0.5 * f_c[-1], "a grid near (1,1), hexagons at high occupancy"
    assert np.all(np.diff(s.interdot_barrier_top_ev) < 0) and s.interdot_barrier_top_ev[0] > s.interdot_barrier_top_ev[-1] + 0.02, "the barrier is higher by tens of meV at (1,1)"
    assert np.all(np.diff(s.tunnel_coupling_ev) > 0) and s.tunnel_coupling_ev[-1] > 1e6 * s.tunnel_coupling_ev[0]
    bare = dataclasses.replace(device.runtime_interdot_barrier(), puddle=None)
    gates_only = bare.barrier_top_ev(s.v_carried_v)
    assert np.all(s.interdot_barrier_top_ev[:6] > gates_only[:6] + 1e-3), "the puddles raise the top above the gates' own action at low occupancy"


def test_the_superdot_front_moves_from_high_occupancy_toward_the_first_electrons_as_b_m_opens(device: Device):
    """At each B_M the plunger diagonal has an occupancy above which the splitting exceeds 0.9; that
    occupancy falls monotonically as B_M rises, and no B_M below 0.9 V merges the (1,1) cell."""
    v0 = np.array(device.programmed_point_v)
    jb, i1, i2 = (device.gate_names.index(g) for g in ("B_M", "P1", "P2"))
    dp = np.linspace(-0.9, 0.3, 241)
    onsets = []
    for b_m in (0.3, 0.4, 0.5, 0.6, 0.8):
        v = np.broadcast_to(v0, dp.shape + v0.shape).copy()
        v[:, jb] = b_m
        v[:, i1] += dp
        v[:, i2] += dp
        s = device.at(v)
        k = np.flatnonzero(s.interdot_splitting > 0.9)
        assert len(k) and k[0] > 0, b_m
        onsets.append(float(s.n_mean[k[0]].sum()))
    assert np.all(np.diff(onsets) < 0), onsets
    assert onsets[0] > 50 and onsets[-1] < 25


def test_without_puddles_the_dots_keep_their_compiled_size(device: Device):
    layout = two_dot_layout(material="SiGe", pitch_nm=80.0)
    op = dict(zip(device.gate_names, device.programmed_point_v, strict=True))
    fixed = compile_device(layout, CompileParams(operating_point_v=op, puddles=False), seed=device.seed)
    assert fixed.puddle is None
    assert fixed.c_dot_af == device.c_dot_af and fixed.tunnel_coupling_ev == device.tunnel_coupling_ev
    v0 = np.array(device.programmed_point_v)
    i1 = device.gate_names.index("P1")
    v = np.stack([v0, v0 - 0.3 * np.eye(len(v0))[i1]])
    s = fixed.at(v)
    assert s.puddle_radii_nm is None
    np.testing.assert_allclose(s.addition_energies_ev[1], s.addition_energies_ev[0], rtol=1e-3), "fixed up to the splitting's effective capacitance"
