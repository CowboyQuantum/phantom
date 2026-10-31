"""P0 for the triple-point splitting set by the interdot conductance (#66, ADR 0032): Golden & Halperin's
weak and strong limits and a monotone curve between them, Kemble's transmission, the CI model's measured
splitting equal to its target, nothing changed without conductance, and the reference device going from
its geometric splitting with B_M closed to a superdot as B_M opens. Analytic limits and invariants only."""

import dataclasses
import math

import numpy as np
import pytest

from phantom.charge.ci import ConstantInteractionParams
from phantom.charge.splitting import (
    F_MAX,
    InterdotSplitting,
    golden_halperin_splitting,
    kemble_transmission,
    strong_splitting,
    weak_splitting,
)
from phantom.contract import reference_device
from phantom.device import Device


@pytest.fixture(scope="module")
def device() -> Device:
    return reference_device()


def test_the_golden_halperin_limits_and_a_monotone_curve_between():
    g = np.linspace(0.0, 1.0, 20001)
    f = golden_halperin_splitting(g)
    assert f[0] == 0.0 and f[-1] == pytest.approx(1.0, abs=1e-12)
    assert np.all(np.diff(f) >= 0), "monotone"
    small = np.array([1e-4, 1e-2, 0.1])
    np.testing.assert_allclose(golden_halperin_splitting(small), 2 * math.log(2) / math.pi**2 * 2 * small, rtol=1e-12), "weak: (2 ln 2/π²)·N_ch·g"
    knots = np.array([0.85, 0.99, 0.999])
    np.testing.assert_allclose(golden_halperin_splitting(knots), strong_splitting(knots), rtol=1e-12), "strong: 1 + (16e^γ/π³)(1−g)ln(1−g)"
    np.testing.assert_allclose(golden_halperin_splitting(0.9), strong_splitting(0.9), rtol=1e-5), "between the strong knots"
    assert weak_splitting(0.15) < golden_halperin_splitting(0.5) < strong_splitting(0.85)


def test_kemble_transmission():
    hw = 4e-3
    assert float(kemble_transmission(0.0, hw)) == 0.5
    du = np.linspace(-0.02, 0.02, 41)
    np.testing.assert_allclose(kemble_transmission(du, hw) + kemble_transmission(-du, hw), 1.0, rtol=1e-15)
    np.testing.assert_allclose(kemble_transmission(0.02, hw), math.exp(-2 * math.pi * 0.02 / hw), rtol=1e-12)


def _symmetric(top_ev: float) -> ConstantInteractionParams:
    ci = ConstantInteractionParams(c_dot_af=np.array([[85.0, -20.0], [-20.0, 85.0]]), c_gate_af=np.array([[12.0, 2.5, 0.3], [2.5, 12.0, 0.3]]))
    return dataclasses.replace(ci, splitting=InterdotSplitting(v_ref_v=(0.0,) * 3, top_ev=top_ev, top_slopes_ev_per_v=(0.0,) * 3, hbar_omega_b_ev=4e-3))


@pytest.mark.parametrize("top_ev", [0.05, 2e-3, 0.0, -2e-3, -0.01])
def test_the_ci_model_s_triple_points_split_by_the_target(top_ev: float):
    """Along the symmetric axis n_g = (x, x) the total charge steps 0 → 1 at x₁ and 1 → 2 at x₂; the
    fractional splitting 2(x₂ − x₁) of the CI energies with the effective capacitance is the target f."""
    ci = _symmetric(top_ev)
    v = np.zeros(3)
    c_inv = np.linalg.inv(ci.local(v)[0])
    x = np.linspace(0.0, 1.0, 200001)
    states = np.array([[0, 0], [1, 0], [0, 1], [1, 1], [2, 0], [0, 2]])
    d = states[None, :, :] - np.stack([x, x], axis=-1)[:, None, :]
    u = 0.5 * np.einsum("nki,ij,nkj->nk", d, c_inv, d)
    total = states.sum(axis=1)[np.argmin(u, axis=1)]
    x1, x2 = x[np.argmax(total >= 1)], x[np.argmax(total >= 2)]
    _, f = ci.interdot_splitting(v)
    assert 2 * (x2 - x1) == pytest.approx(float(f), abs=2e-5)


def test_without_conductance_the_charge_model_is_unchanged():
    ci = _symmetric(1.0)  # 1 eV above the level: g = e^{-2π/0.004}, nothing
    v = np.zeros(3)
    np.testing.assert_array_equal(ci.local(v)[0], ci.c_dot_af)
    g, f = ci.interdot_splitting(v)
    assert float(g) < 1e-300 and float(f) == pytest.approx(2 * 20.0 / (85.0 + 20.0), rel=1e-15)


def test_opening_b_m_takes_the_reference_device_from_its_geometric_splitting_to_a_superdot(device: Device):
    v0 = device.operating_point()
    jb = device.gate_names.index("B_M")
    v = v0 + np.linspace(-0.2, 0.2, 81)[:, None] * np.eye(len(v0))[jb]
    ci = device.runtime_ci()
    g, f = ci.interdot_splitting(v)
    _, f_c = dataclasses.replace(ci, splitting=None).interdot_splitting(v)
    assert np.all(np.diff(g) >= 0) and np.all(f >= f_c - 1e-12)
    rising = g > 1e-2  # with B_M closed f is the geometric splitting, which the shape moves a little
    assert np.all(np.diff(f[rising]) >= -1e-12)
    assert g[0] < 1e-15 and f[0] == pytest.approx(f_c[0], abs=1e-12), "closed: the geometric splitting, a near-grid"
    assert f[-1] == pytest.approx(F_MAX, abs=1e-6) and f[-1] > 0.95, "open: a superdot"
    merged = device.shape.merged(v)
    assert np.all(f[merged] > 0.95)
    s = device.at(device.programmed_point_v)
    assert 0 < float(s.interdot_conductance) < 1 and 0 < float(s.interdot_splitting) < 1
