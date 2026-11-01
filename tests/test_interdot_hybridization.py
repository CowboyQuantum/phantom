"""P0 for the interdot tunnel coupling in the charge states (#58, ADR 0029): the t → 0 limit is the CI
model, an interdot line follows the two-level polarization ½(1 − ε/Ω·tanh(Ω/2kT)), the rate equation over
the eigenstates keeps detailed balance and relaxes a pair at Γ_m at every detuning, t follows the gates
from its compiled value and saturates at ħω₀/π, the barrier gate widens the interdot line of the
reference device, and a session reads eigenstates. Analytic limits and invariants only."""

import math

import numpy as np
import pytest

from phantom.charge.ci import ConstantInteractionParams
from phantom.charge.hybrid import InterdotCoupling, boltzmann
from phantom.constants import K_B_EV_PER_K
from phantom.contract import reference_device
from phantom.device import Device
from phantom.readout.master_equation import RateEquation, project

T_EV = 50e-6
"""A coupling well above kT (13 µeV at 0.15 K), well below the charging energies (meV)."""


@pytest.fixture(scope="module")
def device() -> Device:
    return reference_device()


def _coupled(t_ev: float, ci: ConstantInteractionParams | None = None) -> ConstantInteractionParams:
    """The default CI model with a constant coupling t."""
    ci = ci or ConstantInteractionParams()
    c = InterdotCoupling(v_ref_v=(0.0,) * len(ci.gate_names), t_ref_ev=t_ev, t_max_ev=1e-2, log_slopes_per_v=(0.0,) * len(ci.gate_names))
    return ConstantInteractionParams(ci.gate_names, ci.c_dot_af, ci.c_gate_af, ci.temperature_k, ci.window, c)


def _across_interdot_line(ci: ConstantInteractionParams, half_width: float = 0.05, n: int = 401) -> np.ndarray:
    """Gate points (n, 3) along the detuning n_g + x(1, −1)/2 through the middle of the (1,0)–(0,1) line."""
    tp1, tp2 = ci.triple_points_ng()
    x = np.linspace(-half_width, half_width, n)
    ng = 0.5 * (tp1 + tp2) + x[:, None] * np.array([0.5, -0.5])
    v = np.array([ci.gate_voltages_for_ng(g) for g in ng])
    return np.concatenate([v, np.zeros((n, 1))], axis=1)


def _two_level_n1(ci: ConstantInteractionParams, v: np.ndarray, a: tuple[int, int], b: tuple[int, int], t: float | np.ndarray) -> np.ndarray:
    """⟨N₁⟩ of the pair a = (N₁+1, N₂), b = (N₁, N₂+1) alone: N₁ + ½(1 − ε/Ω tanh(Ω/2kT))."""
    kt = K_B_EV_PER_K * ci.temperature_k
    eps = np.array([float(np.diff(ci.state_energies_ev(p, np.array([b, a])))[0]) for p in v])  # U_a − U_b
    omega = np.sqrt(eps**2 + 4 * t**2)
    return b[0] + 0.5 * (1 - eps / omega * np.tanh(omega / (2 * kt)))


def test_the_t_to_zero_limit_is_the_ci_model():
    """Statics at any point; dynamics away from the resonance |ε| ≲ t, where the eigenstates share the
    electron even for a vanishing t (the projection assumes a dwell long against ħ/Ω, ADR 0029)."""
    ci = ConstantInteractionParams()
    tiny = _coupled(1e-15)
    rng = np.random.default_rng(3)
    v = np.concatenate([rng.uniform(-0.05, 0.25, (500, 2)), np.zeros((500, 1))], axis=1)  # the empty corner included
    np.testing.assert_allclose(tiny.mean_occupation(v), ci.mean_occupation(v), atol=1e-12)
    np.testing.assert_array_equal(tiny.ground_state_occupation(v), ci.ground_state_occupation(v))
    line = _across_interdot_line(ci, n=40)  # no point at ε = 0
    slow, fast = RateEquation(ci, 1e3, 1e3, 1e2), RateEquation(tiny, 1e3, 1e3, 1e2)
    np.testing.assert_allclose(fast.sweep(line, 1e-4)[0], slow.sweep(line, 1e-4)[0], atol=1e-9)


def test_without_a_coupling_the_eigenstates_are_the_charge_states_bit_for_bit():
    ci = ConstantInteractionParams()
    rates = RateEquation(ci, 1e3, 2e3, 5e2)
    v = _across_interdot_line(ci, n=3)[1]
    states = rates.window_states(v)
    u, w = rates.rate_matrix(v, states)
    e, c, w2 = rates.eigen_rate_matrix(v, states)
    np.testing.assert_array_equal(e, u)
    np.testing.assert_array_equal(c, np.eye(len(u)))
    np.testing.assert_array_equal(w2, w)
    p = np.random.default_rng(1).dirichlet(np.ones(len(u)))
    np.testing.assert_array_equal(project(states, c, p, states, c), p / p.sum())


def test_an_interdot_line_follows_the_two_level_polarization():
    """Across the (1,0)–(0,1) line ⟨N₁⟩ is the pair's thermal polarization with gap 2t, up to the
    admixture of the far states, (t/E_c)² ≈ 10⁻⁴; with t ≫ kT it is visibly wider than the CI step."""
    ci = _coupled(T_EV)
    v = _across_interdot_line(ci)
    expected = _two_level_n1(ci, v, (1, 0), (0, 1), T_EV)
    n1 = ci.mean_occupation(v)[:, 0]
    np.testing.assert_allclose(n1, expected, atol=5e-4)
    classical = ConstantInteractionParams().mean_occupation(v)[:, 0]
    assert np.max(np.abs(classical - expected)) > 0.1, "the CI step is far narrower"
    np.testing.assert_allclose(ci.mean_occupation(v).sum(axis=1), 1.0, atol=1e-9), "one electron, shared"


def test_the_rate_equation_over_eigenstates_keeps_detailed_balance():
    ci = _coupled(T_EV)
    rates = RateEquation(ci, 1e6, 3e5, 2e7)
    for v in _across_interdot_line(ci, n=7):
        states = rates.window_states(v)
        e, c, w = rates.eigen_rate_matrix(v, states)
        p = boltzmann(e, rates.kt_ev)
        assert np.abs(w @ p).max() <= 1e-9 * np.abs(w).max()
        np.testing.assert_allclose(c.T @ c, np.eye(len(e)), atol=1e-12)
        _, pop = rates.stationary(v)
        np.testing.assert_allclose(pop @ states, ci.mean_occupation(v), atol=1e-9), "the stationary state is the thermal one"


def test_a_pair_relaxes_at_the_interdot_rate_at_every_detuning():
    """For the two states of an interdot line alone, g = 1 for any mixing, so the pair's total rate
    is Γ_m: Γ_m·f(Δ/kT) + Γ_m·f(−Δ/kT)."""
    ci = _coupled(T_EV)
    g_m = 2e7
    rates = RateEquation(ci, 0.0, 0.0, g_m)
    for v in _across_interdot_line(ci, n=9):
        pair = np.array([[0, 1], [1, 0]])
        _, _, w = rates.eigen_rate_matrix(v, pair)
        assert w[0, 1] + w[1, 0] == pytest.approx(g_m, rel=1e-12)


def test_the_coupling_follows_the_gates_from_its_compiled_value(device: Device):
    c = device.runtime_coupling()
    v0 = device.operating_point()
    assert c.tunnel_coupling_ev(v0) == pytest.approx(device.tunnel_coupling_ev, rel=1e-12)
    inter = next(b for b in device.barriers if b.between == ("dot1", "dot2"))
    np.testing.assert_allclose(c.log_slopes_per_v, 0.5 * np.array(inter.log_rate_slopes_per_v), rtol=1e-15)
    j = device.gate_names.index(device.barrier_gate)
    along = v0 + np.linspace(-0.2, 1.0, 241)[:, None] * np.eye(len(v0))[j]
    t = c.tunnel_coupling_ev(along)
    assert np.all(np.diff(t) >= 0) and t[-1] > 1e4 * t[0], "opening the barrier gate raises t"
    # t = t_max·√T of the interdot barrier's one transmission (#70): along the barrier gate the top sinks
    # through the tunnelling level and far below it, so T → 1 and t reaches the ceiling ħω₀/π, which
    # bounds every gate
    assert np.all(t <= c.t_max_ev) and t[-1] == pytest.approx(c.t_max_ev, rel=1e-9), "to ħω₀/π as the barrier sinks"
    assert c.barrier is not None and np.all(t == c.t_max_ev * np.sqrt(c.barrier.transmission(along)))
    other = v0 + np.linspace(0.0, 3.0, 31)[:, None] * np.eye(len(v0))[device.gate_names.index("P1")]
    assert c.tunnel_coupling_ev(other)[-1] == pytest.approx(c.t_max_ev, rel=1e-6), "a plunger far up takes it to the ceiling too"
    mass = 0.19
    hbar_omega = math.sqrt(2 * 0.0380998212 * device.dots[0].curvature_ev_per_nm2 / mass)
    assert c.t_max_ev == pytest.approx(hbar_omega / math.pi, rel=1e-6)
    s = device.at(np.array(device.programmed_point_v))
    assert float(s.tunnel_coupling_ev) == pytest.approx(float(c.tunnel_coupling_ev(s.v_carried_v)), rel=1e-15)


def _interdot_point(device: Device, v0: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(the middle of the interdot segment of the honeycomb cell around `v0`, the detuning direction):
    the CI energies depend on N − n_g alone, so the segment's middle is floor(n_g) + ½(tp₁ + tp₂) in
    induced charge, reached with P1 and P2; along the direction (unit, programmed volts) n_g moves by
    (½, −½), one electron from dot 2 to dot 1 at a fixed total."""
    ci = device.runtime_ci()
    j = [device.gate_names.index(g) for g in ("P1", "P2")]

    def ng(v: np.ndarray) -> np.ndarray:
        return ci.induced_charge(device.effective_voltages(v))

    jac = np.stack([ng(v0 + np.eye(len(v0))[k]) - ng(v0) for k in j], axis=1)  # electrons per volt, linear
    tp1, tp2 = ci.triple_points_ng()
    v = v0.copy()
    v[j] += np.linalg.solve(jac, np.floor(ng(v0)) + 0.5 * (tp1 + tp2) - ng(v0))
    d = np.zeros(len(v0))
    d[j] = np.linalg.solve(jac, np.array([0.5, -0.5]))
    return v, d / np.linalg.norm(d)


def test_the_barrier_gate_widens_the_interdot_line(device: Device):
    """Opening B_M raises t, and the interdot line of the reference device widens from the thermal step
    (0.2 mV) towards the two-level polarization with gap 2t: its 10–90 % width grows with B_M, and ⟨N₁⟩
    follows the pair's formula up to the far states' admixture, (t/E_add)²."""
    jb = device.gate_names.index("B_M")
    e_add = float(np.min(device.runtime_ci().addition_energies_ev()))
    widths, couplings = [], []
    for db in (-0.05, 0.0, 0.016):
        v0 = np.array(device.programmed_point_v)
        v0[jb] += db
        centre, d = _interdot_point(device, v0)
        x = np.linspace(-0.004, 0.004, 801)
        s = device.at(centre + x[:, None] * d)
        n1 = s.n_mean[:, 0] - np.floor(s.n_mean[400, 0])
        widths.append(abs(x[np.argmin(np.abs(n1 - 0.9))] - x[np.argmin(np.abs(n1 - 0.1))]))
        t = float(s.tunnel_coupling_ev[400])
        couplings.append(t)
        base = np.floor(s.n_mean[400]).astype(int)  # the pair shares one electron above it
        a, b = (base[0] + 1, base[1]), (base[0], base[1] + 1)
        expected = _two_level_n1(device.runtime_ci(), s.v_carried_v, a, b, s.tunnel_coupling_ev)
        np.testing.assert_allclose(s.n_mean[:, 0], expected, atol=4 * (t / e_add) ** 2 + 1e-9)
    assert couplings[0] < K_B_EV_PER_K * device.temperature_k < couplings[2], "from below kT to above it"
    assert widths[0] < widths[1] < widths[2]
    assert widths[2] > 2.5 * widths[0]


def test_a_session_reads_the_eigenstates(device: Device):
    """At the middle of an interdot line with t > kT, a `sampled` session draws eigenstates: their total
    charge is an integer and their split is the bonding state's, near half an electron each way; a long
    `mean` dwell gives the thermal ⟨N⟩."""
    jb = device.gate_names.index("B_M")
    v0 = np.array(device.programmed_point_v)
    v0[jb] += 0.016
    centre, _ = _interdot_point(device, v0)
    s = device.at(centre)
    assert float(s.tunnel_coupling_ev) > 2 * K_B_EV_PER_K * device.temperature_k
    n = device.session(noise=False).sense(centre, 1e-3, samples=200).n_present
    np.testing.assert_allclose(n.sum(axis=1), np.rint(n.sum(axis=1)), atol=1e-9)
    frac = n[:, 0] - np.floor(n[:, 0])
    assert np.all(np.abs(frac - 0.5) < 0.2), "each sample reads an eigenstate's split, not 0 or 1"
    m = device.session(noise=False, charge="mean").sense(centre, 1.0).n_present
    np.testing.assert_allclose(m, s.n_mean, atol=1e-9)


def test_the_middle_of_an_interdot_segment_has_no_seam(device: Device):
    """The window is centred on rint(n_g), which changes at the middle of every interdot segment
    (n_g = N + ½). With t near a meV a chain truncated at ±2 states would leave a seam of 2·10⁻⁵ e there;
    the statics' ±4 window and the dynamics' window chosen by t keep it below 10⁻¹² (ADR 0014, 0029)."""
    jb = device.gate_names.index("B_M")
    v0 = np.array(device.programmed_point_v)
    v0[jb] += 0.08  # with the screened barrier (#70) t passes 0.5 meV about 80 mV up
    centre, d = _interdot_point(device, v0)
    vc = device.effective_voltages(centre)
    ci = device.runtime_ci()
    assert np.allclose(ci.induced_charge(vc) % 1.0, 0.5, atol=1e-9), "on the tie"
    assert float(ci.tunnel_coupling_ev(vc)) > 5e-4
    step = 1e-12 * d
    below, above = ci.mean_occupation(device.effective_voltages(centre - step)), ci.mean_occupation(device.effective_voltages(centre + step))
    np.testing.assert_allclose(above, below, atol=1e-10)
    m = device.session(noise=False, charge="mean").sense(centre, 1.0).n_present
    np.testing.assert_allclose(m, device.at(centre).n_mean, atol=1e-9), "the dynamics' window agrees with the statics'"
