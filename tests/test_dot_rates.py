"""P0 for the qubit dots' gate-dependent rates (#53, ADR 0027): the runtime rates against the compile tier
at and near the operating point, each barrier gate controlling its own rate, detailed balance far from the
anchor, and a session whose dynamics telegraph and latch once dot 1 is isolated and reduce to independent
Boltzmann draws when every rate is fast. With the interdot coupling (#58) dot 1 is isolated only when B_M
closes too: with B_L closed alone it loads through dot 2's lead at Γ_R (t/Δ)². Analytic limits and
invariants only."""

import math

import numpy as np
import pytest

from phantom.contract import reference_device
from phantom.device import Device
from phantom.electrostatics.dots import DerivedDot, barrier_action

CLOSED_B_L_V = -0.1
"""B_L below the programmed point that slows the left lead into the 10 kHz range on the reference device."""
CLOSED = {"B_L": -0.165, "B_M": -0.12}
"""B_L and B_M below the programmed point: dot 1 isolated, its lead in the 10 kHz range and t below
1 µeV, so nothing reaches it through dot 2 (#58). With the screened interdot barrier (#70) t falls by a
factor e every 23 mV of B_M, so B_M closes by 0.12 V; closing B_M also moves the left lead (#48), so
B_L closes further than `CLOSED_B_L_V`."""


@pytest.fixture(scope="module")
def device() -> Device:
    return reference_device()


def _derived(d) -> DerivedDot:
    return DerivedDot(d.name, d.x_nm, d.y_nm, d.energy_ev, d.curvature_ev_per_nm2, d.radius_nm, np.array(d.lever_arms))


def _compile_tier_rates_hz(device: Device, v: np.ndarray) -> np.ndarray:
    """(Γ_L, Γ_R, Γ_m) recomputed from the WKB action at each barrier's fixed tunnelling level."""
    lay, basis = device.basis()
    d1, d2 = (_derived(d) for d in device.dots)
    res = lay.reservoir_points_nm()
    ends = ((d1, res["R_L"]), (d2, res["R_R"]), (d1, (d2.x_nm, d2.y_nm)))
    out = []
    for b, (dot, point) in zip(device.barriers, ends, strict=True):
        action = barrier_action(basis, v, dot, point, 0.19, level_ev=b.tunnelling_level_ev)[2]
        out.append(b.attempt_rate_hz * math.exp(-2.0 * action))
    return np.array(out)


def _point(device: Device, **shift: float) -> np.ndarray:
    v = np.array(device.programmed_point_v)
    for gate, dv in shift.items():
        v[device.gate_names.index(gate)] += dv
    return v


def _dot1_degeneracy(device: Device, v0: np.ndarray, n1_from: int | None = None) -> np.ndarray:
    """The programmed point along P1 from `v0` where dot 1 is half-way through a transition to its lead;
    with `n1_from`, the transition from that many electrons, searched farther: the barrier gates move
    dot 1's lines along P1 as they reshape the dots (#59)."""
    j = device.gate_names.index("P1")
    half, points = (0.04, 16001) if n1_from is None else (0.1, 40001)  # 5 µV steps: the thermal step is about 0.2 mV wide
    grid = v0 + np.linspace(-half, half, points)[:, None] * np.eye(len(v0))[j]
    n = device.at(grid).n_mean
    frac = n[:, 0] - np.floor(n[:, 0])
    flat2 = np.abs(n[:, 1] - np.rint(n[:, 1])) < 1e-2  # dot 2 settled, up to its interdot admixture (#58)
    if n1_from is not None:
        flat2 &= np.floor(n[:, 0]) == n1_from
    k = int(np.argmin(np.where(flat2, np.abs(frac - 0.5), np.inf)))
    assert abs(frac[k] - 0.5) < 0.05
    return grid[k]


def _compile_tier_interdot_rate_hz(device: Device, v: np.ndarray) -> float:
    """Γ_m recomputed as f_a·T, T Kemble's transmission of the interdot barrier's top from the potential
    with the dots held, moved from the compiled top by each gate's screening factor (#70, ADR 0034)."""
    from phantom.charge.barrier import kemble_transmission
    from phantom.compile.device import interdot_barrier_top

    _, basis = device.basis()
    d1, d2 = (_derived(d) for d in device.dots)
    inter = device.barriers[2]
    top0 = device.interdot_barrier_top_ev
    bare = interdot_barrier_top(basis, v, d1, d2, inter.tunnelling_level_ev, 0.19)[0] - top0
    # one gate moves at a time here, so one factor applies: the one of the gate that moved
    moved = np.flatnonzero(np.asarray(v) != device.operating_point())
    sigma = device.interdot_barrier_screening[int(moved[0])] if len(moved) else 1.0
    # the puddles' edges (#72): the occupancy's own shift of the top, from the stored cut
    shift = 0.0 if device.puddle is None else float(device.puddle.barrier_shift_ev(device.runtime_ci().occupancy_scale(np.asarray(v, dtype=float))))
    return inter.attempt_rate_hz * float(kemble_transmission(top0 + sigma * bare + shift, device.interdot_barrier_hbar_omega_ev))


def test_the_runtime_rates_are_the_compiled_ones_at_the_anchor_and_follow_the_compile_tier_nearby(device: Device):
    rates = device.runtime_rates()
    v0 = device.operating_point()
    compiled = [b.rate_hz for b in device.barriers]
    np.testing.assert_allclose(rates.rates_hz(v0), compiled, rtol=1e-12)
    np.testing.assert_allclose(_compile_tier_rates_hz(device, v0)[:2], compiled[:2], rtol=1e-12), "the lead rates: the WKB action"
    assert _compile_tier_interdot_rate_hz(device, v0) == pytest.approx(compiled[2], rel=1e-12), "the interdot rate: Kemble's transmission (#70)"
    for gate in ("P1", "P2", "R_L", "R_R"):  # the barrier gates move the dots: the shape grid (#59) below
        for dv in (-0.01, -0.003, 0.003, 0.01):
            v = v0.copy()
            v[device.gate_names.index(gate)] += dv
            # up to a factor e every 4 mV on these gates: the linearised exponent is within 20 % of the
            # action's change, plus the grid's step in the action (about 0.1 in ln Γ on the 2 nm grid)
            tier = _compile_tier_rates_hz(device, v)[:2]
            change = np.log(tier / np.array(compiled[:2]))
            error = np.log(rates.rates_hz(v)[:2] / tier)
            assert np.all(np.abs(error) <= 0.1 + 0.2 * np.abs(change)), f"{gate} {dv:+}: {error} against {change}"
            # the interdot top is linear in the gates with the dots held: the screened slope reproduces it
            assert rates.rates_hz(v)[2] == pytest.approx(_compile_tier_interdot_rate_hz(device, v), rel=0.05), f"{gate} {dv:+}"
    # along a barrier gate the lead exponents are the shape grid's, the actions with the dots moved, at its
    # nodes, and the interdot exponent is −2πσ·Δ(top)/ħω_b of the grid's barrier top (#70)
    g = device.shape
    f_a = np.array([b.attempt_rate_hz for b in device.barriers])
    for gate in ("B_L", "B_M", "B_R"):
        j = device.gate_names.index(gate)
        k = g.gate_indices.index(j)
        for node, row in zip(g.nodes_v[k], g.changes[k], strict=True):
            v = v0.copy()
            v[j] = node
            x = np.log(compiled) - np.log(f_a - compiled) + np.array(row)[g.layout["log_rate"]]
            shift = 0.0 if device.puddle is None else float(device.puddle.barrier_shift_ev(device.runtime_ci().occupancy_scale(v)))  # the puddles' edges (#72)
            x[2] = np.log(compiled[2]) - np.log(f_a[2] - compiled[2]) - 2 * np.pi * (device.interdot_barrier_screening[j] * np.array(row)[g.layout["barrier_top_ev"]][0] + shift) / device.interdot_barrier_hbar_omega_ev
            np.testing.assert_allclose(rates.rates_hz(v), f_a * (1 / (1 + np.exp(-x))), rtol=1e-9)


def test_the_state_carries_the_rates_at_its_voltages(device: Device):
    v = np.stack([_point(device), _point(device, B_L=CLOSED_B_L_V), _point(device, B_M=-0.15)])
    s = device.at(v)
    expected = device.runtime_rates().rates_hz(device.effective_voltages(v))
    for k, name in enumerate(("gamma_left_hz", "gamma_right_hz", "gamma_interdot_hz")):
        assert getattr(s, name).shape == (3,)
        np.testing.assert_array_equal(getattr(s, name), expected[:, k])
    assert s.gamma_left_hz[1] < 1e-4 * s.gamma_left_hz[0], "closing B_L slows the left lead"
    assert s.gamma_interdot_hz[2] < 1e-3 * s.gamma_interdot_hz[0], "closing B_M slows interdot transfer"


def test_each_barrier_gate_controls_its_own_rate_and_a_rate_saturates_at_the_attempt_rate(device: Device):
    left, right, inter = device.barriers
    j = {g: device.gate_names.index(g) for g in ("B_L", "B_M", "B_R")}
    assert left.log_rate_slopes_per_v[j["B_L"]] > right.log_rate_slopes_per_v[j["B_L"]] > 0
    assert right.log_rate_slopes_per_v[j["B_R"]] > left.log_rate_slopes_per_v[j["B_R"]] > 0
    assert inter.log_rate_slopes_per_v[j["B_M"]] > max(left.log_rate_slopes_per_v[j["B_M"]], right.log_rate_slopes_per_v[j["B_M"]])
    # Γ_m ∝ T and t ∝ √T along the same cut: the interdot rate and t share the barrier gate's scale
    assert inter.log_rate_slopes_per_v[j["B_M"]] == pytest.approx(2.0 / device.barrier_v_scale_v, rel=1e-6)
    rates = device.runtime_rates()
    v0 = device.operating_point()
    shut = v0.copy()
    shut[j["B_L"]] -= 0.05
    ratio = rates.rates_hz(shut)[0] / rates.rates_hz(v0)[0]
    assert ratio == pytest.approx(math.exp(-0.05 * left.log_rate_slopes_per_v[j["B_L"]]), rel=0.3), "exponential while closed, near the compiled slope"
    # one barrier gate opened until its action vanishes (T = 1) gives f_a·σ(ln(T₀/(1 − T₀)) − ln T₀) =
    # f_a/(2 − T₀); the barrier gates' changes add (#59), so all three open approach f_a, never past it
    f_a = np.array([b.attempt_rate_hz for b in device.barriers])
    t0 = compiled_t0 = np.array([b.rate_hz for b in device.barriers]) / f_a
    one = v0.copy()
    one[j["B_L"]] += 1.0
    assert rates.rates_hz(one)[0] == pytest.approx(f_a[0] / (2 - t0[0]), rel=0.1), "the grid's action is the refined dot's"
    wide = v0.copy()
    wide[[j["B_L"], j["B_M"], j["B_R"]]] += 1.0
    assert np.all(rates.rates_hz(wide) <= f_a) and np.all(rates.rates_hz(wide) > 0.5 * f_a), compiled_t0


def test_detailed_balance_holds_far_from_the_anchor(device: Device):
    """Only the prefactors move: W[j, i]·p_i = W[i, j]·p_j with the Boltzmann p, at a closed barrier
    where the rates span many decades."""
    rates = device.runtime_rates()
    v = device.effective_voltages(_dot1_degeneracy(device, _point(device, B_L=CLOSED_B_L_V, B_M=-0.12)))
    states = rates.window_states(v)
    _, w = rates.rate_matrix(v, states)
    u = device.runtime_ci().state_energies_ev(v, states)
    boltzmann = np.exp(-(u - u.min()) / rates.kt_ev)
    flux = w * boltzmann[None, :]
    off = ~np.eye(len(u), dtype=bool)
    np.testing.assert_allclose(flux[off], flux.T[off], rtol=1e-9, atol=1e-12 * np.abs(w).max())
    assert np.abs(w @ boltzmann).max() <= 1e-9 * np.abs(w).max() * boltzmann.max()


def _isolated_centre(device: Device) -> np.ndarray:
    """Dot 1's transition with `CLOSED`, on the line the B_L-alone transition is on."""
    alone = _dot1_degeneracy(device, _point(device, B_L=CLOSED_B_L_V))
    return _dot1_degeneracy(device, _point(device, **CLOSED), n1_from=int(np.floor(device.at(alone).n_mean[0])))


def _pair_switching_hz(device: Device, v: np.ndarray) -> float:
    """The total rate between the two lowest eigenstates at programmed `v`, Hz."""
    rates = device.runtime_rates()
    vc = device.effective_voltages(v)
    e, _, w = rates.eigen_rate_matrix(vc, rates.window_states(vc))
    a, b = np.argsort(e)[:2]
    return float(w[a, b] + w[b, a])


@pytest.mark.parametrize("gamma_dwell", [0.02, 0.5])
def test_an_isolated_dot_gives_the_two_state_switching_statistics(device: Device, gamma_dwell: float):
    """B_L and B_M closed put Γ_L in the 10 kHz range at dot 1's transition and leave no path through
    dot 2: the charge is a two-state Markov process with total rate Γ_L, changing between samples with
    probability 2p(1−p)(1 − e^{−Γ_L·dwell})."""
    v = _isolated_centre(device)
    s = device.at(v)
    gamma = float(s.gamma_left_hz)
    assert 1e3 < gamma < 1e5
    assert _pair_switching_hz(device, v) == pytest.approx(gamma, rel=0.02), "nothing but the left lead"
    n = np.rint(device.session(noise=False).sense(v, gamma_dwell / gamma, samples=6000).n_present[:, 0])
    p = float(np.mean(n > n.min()))
    changes = float(np.mean(n[1:] != n[:-1]))
    expected = 2 * p * (1 - p) * (1 - math.exp(-gamma_dwell))
    assert set(np.unique(n)) <= {n.min(), n.min() + 1}
    assert abs(changes - expected) < 6 * math.sqrt(expected / 6000) + 0.002


def test_with_b_l_closed_alone_dot_1_loads_through_dot_2(device: Device):
    """With B_L closed and B_M open, dot 2's lead loads dot 1 in second order in t, as elastic
    cotunnelling: the electron-like path through (N₁, N₂+1), Δ above (N₁+1, N₂), and the hole-like one
    through (N₁+1, N₂−1), Δ' above (N₁, N₂), add in amplitude. The pair switches at
    Γ_L + Γ_R t²(1/Δ + 1/Δ')², orders of magnitude faster than Γ_L alone (#58)."""
    v = _dot1_degeneracy(device, _point(device, B_L=CLOSED_B_L_V))
    s = device.at(v)
    vc = s.v_carried_v
    n1, n2 = (int(x) for x in np.floor(s.n_mean))
    u = device.runtime_ci().state_energies_ev(vc, np.array([[n1, n2], [n1 + 1, n2], [n1, n2 + 1], [n1 + 1, n2 - 1]]))
    t = float(s.tunnel_coupling_ev)
    delta, delta_h = float(u[2] - u[1]), float(u[3] - u[0])
    assert min(delta, delta_h) > 50 * t
    expected = float(s.gamma_left_hz) + float(s.gamma_right_hz) * t**2 * (1 / delta + 1 / delta_h) ** 2
    assert expected > 30 * float(s.gamma_left_hz)
    assert _pair_switching_hz(device, v) == pytest.approx(expected, rel=0.02)


def test_with_every_rate_fast_the_session_reduces_to_independent_boltzmann_draws(device: Device):
    v = _dot1_degeneracy(device, _point(device))
    s = device.at(v)
    slowest = min(float(s.gamma_left_hz), float(s.gamma_right_hz), float(s.gamma_interdot_hz))
    p = float(s.n_mean[0] - np.floor(s.n_mean[0]))
    for dynamics in (True, False):
        n = device.session(noise=False, dynamics=dynamics).sense(v, 50.0 / slowest, samples=6000).n_present[:, 0]
        upper = np.rint(n) > np.floor(s.n_mean[0])
        assert abs(upper.mean() - p) < 5 * math.sqrt(p * (1 - p) / 6000)
        changes = float(np.mean(upper[1:] != upper[:-1]))
        assert abs(changes - 2 * p * (1 - p)) < 6 * math.sqrt(2 * p * (1 - p) / 6000), "uncorrelated samples"


def _crossing_v(p1: np.ndarray, n1: np.ndarray, level: float) -> float:
    k = int(np.flatnonzero(np.diff(np.sign(n1 - level)) != 0)[0])
    return float(p1[k] + (level - n1[k]) * (p1[k + 1] - p1[k]) / (n1[k + 1] - n1[k]))


@pytest.mark.parametrize("closed", [False, True], ids=["open", "closed"])
def test_a_closed_barrier_latches_the_charge_along_a_sweep(device: Device, closed: bool):
    """A P1 sweep across dot 1's transition at 3 µs per point, up then down, in `mean` mode: with dot 1
    isolated (B_L and B_M closed, Γ_L about 10⁴ Hz there) the charge lags the gate and the transition
    opens into a hysteresis loop a few mV wide; with the barriers open (Γ_L about 10¹¹ Hz) the two
    directions agree."""
    centre = _isolated_centre(device) if closed else _dot1_degeneracy(device, _point(device))
    j = device.gate_names.index("P1")
    p1 = centre[j] + np.linspace(-0.01, 0.01, 201)
    up = np.broadcast_to(centre, p1.shape + centre.shape).copy()
    up[:, j] = p1
    run = device.session(noise=False, charge="mean")
    run.hold(up[0], 1.0)
    n_up = run.sense(up, 3e-6).n_present[:, 0]
    n_down = run.sense(up[::-1], 3e-6).n_present[::-1, 0]
    half = math.floor(float(device.at(centre).n_mean[0])) + 0.5
    width = _crossing_v(p1, n_up, half) - _crossing_v(p1, n_down, half)
    if closed:
        assert width > 2e-3
    else:
        assert abs(width) < 1e-5
