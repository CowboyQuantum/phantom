"""Step 6 protocols: exchange, the fingerprint and house cells, spin resonance, Rabi and Ramsey, and
Pauli spin blockade readout with the valley draw."""

from __future__ import annotations

import numpy as np

from gallery.runner import Context, Protocol
from phantom.charge import Barrier, ConstantInteractionParams
from phantom.charge.hubbard import HubbardParams
from phantom.constants import HBAR_EV_S
from phantom.disorder import NoiseParams, draw_noise
from phantom.spin import (
    H_EV_S,
    MU_B_EV_PER_T,
    Lindblad,
    SpinParams,
    detuning_noise_sigma_ev,
    detuning_voltages,
    draw_spin,
    exchange_j_ev,
    exchange_second_order_ev,
    exchange_sensitivity_ev_per_ev,
    probability,
    psb_contrast,
    psb_occupation,
    quasi_static_average,
    quasi_static_sigma_ev,
    rabi_probability,
    ramsey_probability,
    rotating_frame_hamiltonian_ev,
    singlet_return_probability,
    symmetric_point_v,
    t2_star_s,
)
from phantom.valley import ValleyParams, draw_valley_splitting_uev

GHZ_PER_EV = 1e-9 / H_EV_S


def _axis(spec: dict) -> np.ndarray:
    return np.linspace(float(spec["from"]), float(spec["to"]), int(spec["points"]))


def _pair(p: Protocol) -> tuple[HubbardParams, np.ndarray, float]:
    """The default double dot as a Hubbard model with the barrier gate `B`, its symmetric point, and the
    rms gate noise of the device's noise draw (the whole 1/f band, quasi-static for a spin experiment)."""
    inst = p.instance
    barrier = Barrier("B", t0_ev=float(inst.get("t0_uev", 20.0)) * 1e-6, v_scale_v=float(inst.get("v_scale_v", 0.05)))
    hub = HubbardParams.from_ci(ConstantInteractionParams(), barrier)
    v0 = symmetric_point_v(hub, {"B": float(inst.get("v_barrier_v", 0.0))})
    seed = int(inst.get("seed", 0))
    noise = draw_noise(NoiseParams.for_material(inst.get("material", "SiGe")), np.random.default_rng(seed), seed)
    return hub, v0, float(noise.one_over_f_amplitude_v)


def _anticrossing_ev(hub: HubbardParams) -> float:
    return float(hub.onsite_u_ev[1] - hub.mutual_v_ev)


# --- exchange ------------------------------------------------------------------------------------


def exchange_vs_detuning(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """J(ε) from the Hubbard two-electron sector at three barrier settings, with the second-order formula
    at the middle one and dJ/dε, the charge-noise sensitivity."""
    with ctx.compile():
        hub, v0, _ = _pair(p)
        eac = _anticrossing_ev(hub)
        vbs = [float(x) for x in p.instance.get("v_barrier_values_v", [-0.05, 0.0, 0.05])]
    eps = _axis(p.sweep[0]) * eac
    out: dict[str, np.ndarray] = {"eps_mev": eps * 1e3, "eps_over_anticrossing": eps / eac}
    with ctx.run(len(eps)):
        for k, vb in enumerate(vbs):
            v = detuning_voltages(hub, v0, eps)
            v[..., hub.gate_names.index("B")] = vb
            out[f"j_vb{k}_ghz"] = exchange_j_ev(hub, v) * GHZ_PER_EV
            if k == 1:
                out["j_second_order_ghz"] = exchange_second_order_ev(hub, v) * GHZ_PER_EV
                out["dj_deps_vb1"] = exchange_sensitivity_ev_per_ev(hub, v0, eps)
        ctx.count(model_evals=len(eps) * (len(vbs) + 2))
    return out


def exchange_vs_barrier(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """J at the symmetric point against the interdot barrier voltage: exponential with slope 2/v_scale."""
    with ctx.compile():
        hub, v0, _ = _pair(p)
    vb = _axis(p.sweep[0])
    v = np.broadcast_to(v0, (len(vb), len(hub.gate_names))).copy()
    v[:, hub.gate_names.index("B")] = vb
    with ctx.run(len(vb)):
        j = exchange_j_ev(hub, v)
        ctx.count(model_evals=len(vb))
    j_law = j[0] * np.exp(2 * (vb - vb[0]) / hub.barrier.v_scale_v)
    return {"v_barrier_v": vb, "j_ghz": j * GHZ_PER_EV, "exponential_law_ghz": j_law * GHZ_PER_EV, "t_uev": hub.barrier.tunnel_coupling_ev(vb) * 1e6}


def exchange_fingerprint(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Singlet return against detuning and evolution time, averaged over the quasi-static detuning noise
    of the device's draw: the oscillation frequency J(ε)/h rises toward the anticrossing and the
    visibility decays fastest where dJ/dε is largest."""
    with ctx.compile():
        hub, v0, sigma_v = _pair(p)
        eac = _anticrossing_ev(hub)
        sigma_eps = detuning_noise_sigma_ev(hub, sigma_v)
        nodes = int(p.instance.get("quadrature_nodes", 31))
    eps = _axis(p.sweep[0]) * eac
    tau = _axis(p.sweep[1]) * 1e-9
    with ctx.run(len(eps) * len(tau)):

        def f(delta: np.ndarray) -> np.ndarray:
            j = exchange_j_ev(hub, detuning_voltages(hub, v0, eps[:, None] + delta[None, :]))  # (eps, nodes)
            return singlet_return_probability(tau[None, :, None], j[:, None, :])  # (eps, tau, nodes)

        p_s = quasi_static_average(f, sigma_eps, nodes)
        ctx.count(model_evals=len(eps) * len(tau))  # the Hubbard evaluations per column (nodes) are the per-point cost
    j0 = exchange_j_ev(hub, detuning_voltages(hub, v0, eps))
    return {"eps_mev": eps * 1e3, "tau_ns": tau * 1e9, "singlet_probability": p_s, "j_ghz": j0 * GHZ_PER_EV, "sigma_eps_uev": np.array(sigma_eps * 1e6)}


def exchange_house_cells(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Singlet return against the two plungers at fixed evolution time with the interdot barrier open
    (J high): inside (1,1) the return probability is cos(J τ/ħ) averaged over the detuning noise, so
    cells of constant return are bounded by contours of constant J and by the charge transitions where
    the pair leaves (1,1) (a singlet in (0,2) or (2,0) returns as a singlet; other charge states carry
    no spin information and read 0.5)."""
    with ctx.compile():
        hub, v0, sigma_v = _pair(p)
        sigma_eps = detuning_noise_sigma_ev(hub, sigma_v)
        tau = float(p.instance.get("tau_ns", 2.0)) * 1e-9
        nodes = int(p.instance.get("quadrature_nodes", 7))
        half = float(p.sweep[0]["half_span_v"])
        n = int(p.sweep[0]["points"])
        i1, i2 = hub.gate_names.index("P1"), hub.gate_names.index("P2")
    p1 = np.linspace(v0[i1] - half, v0[i1] + half, n)
    p2 = np.linspace(v0[i2] - half, v0[i2] + half, n)
    v = np.broadcast_to(v0, (n, n, len(hub.gate_names))).copy()
    v[..., i1] = p1[:, None]
    v[..., i2] = p2[None, :]
    with ctx.run(n * n):
        occ = hub.ground_state_occupation(v)
        n1, n2 = np.rint(occ[..., 0]), np.rint(occ[..., 1])
        in_11 = (n1 == 1) & (n2 == 1)
        doubly = ((n1 == 0) & (n2 == 2)) | ((n1 == 2) & (n2 == 0))

        def f(delta: np.ndarray) -> np.ndarray:
            j = exchange_j_ev(hub, detuning_voltages(hub, v[..., None, :], delta[None, None, :]))  # (n, n, nodes)
            return singlet_return_probability(tau, j)

        p_s = quasi_static_average(f, sigma_eps, nodes)
        p_s = np.where(in_11, p_s, np.where(doubly, 1.0, 0.5))
        ctx.count(model_evals=n * n * (nodes + 1))
    j_map = exchange_j_ev(hub, v) * GHZ_PER_EV
    return {"v_p1_v": p1, "v_p2_v": p2, "singlet_probability": p_s, "j_ghz": np.where(in_11, j_map, np.nan), "n_total": n1 + n2}


# --- single-spin resonance, Rabi, Ramsey -------------------------------------------------------------


def _spins(p: Protocol) -> tuple[SpinParams, tuple[float, ...], float]:
    inst = p.instance
    params = SpinParams(**inst.get("spin", {}))
    seed = int(inst.get("seed", 0))
    draw = draw_spin(params, np.random.default_rng(seed), seed)
    noise = draw_noise(NoiseParams.for_material(inst.get("material", "SiGe")), np.random.default_rng(seed), seed)
    return params, draw.g_factors, float(noise.one_over_f_amplitude_v)


def esr_line(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Spin-flip probability against field and drive frequency after a π burst at the Rabi strength,
    the frequency axis measured from the g = 2 line (2 μ_B B / h): one line per dot whose offset grows
    linearly with field at (g_i − 2) μ_B / h, so the two lines separate by the g-factor difference of
    the draw; their width is the Rabi power broadening."""
    with ctx.compile():
        _, g, _ = _spins(p)
        omega = H_EV_S * float(p.instance.get("f_rabi_hz", 5e6))
        tau_pi = np.pi * HBAR_EV_S / omega
    b = _axis(p.sweep[0])
    df = _axis(p.sweep[1])
    with ctx.run(len(b) * len(df)):
        flip = np.zeros((len(b), len(df)))
        for gi in g:
            delta = (gi - 2.0) * MU_B_EV_PER_T * b[:, None] - H_EV_S * df[None, :]
            flip += rabi_probability(tau_pi, delta, omega)
        ctx.count(model_evals=len(b) * len(df) * len(g))
    return {"b_t": b, "f_offset_mhz": df * 1e-6, "flip_probability": np.clip(flip, 0, 1), "g_factors": np.array(g), "g2_line_ghz": 2.0 * MU_B_EV_PER_T * b / H_EV_S * 1e-9}


def rabi_chevron(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Spin-flip probability of dot 1 against drive detuning and burst time from the Lindblad propagator
    with the Markovian T₁ and T_φ of the parameters: the chevron, narrowing and slowing off resonance."""
    with ctx.compile():
        params, _, _ = _spins(p)
        omega = H_EV_S * float(p.instance.get("f_rabi_hz", 1e7))
    df = _axis(p.sweep[0])
    tau = _axis(p.sweep[1]) * 1e-9
    dt = tau[1] - tau[0]
    up, down = np.array([1, 0], dtype=complex), np.array([0, 1], dtype=complex)
    with ctx.run(len(df) * len(tau)):
        flip = np.zeros((len(df), len(tau)))
        for k, d in enumerate(df):
            h = rotating_frame_hamiltonian_ev(np.array([H_EV_S * d]), np.array([omega]))
            rho = Lindblad.for_spins(h, t1_s=params.t1_s, t_phi_s=params.t_phi_markov_s).evolve(np.outer(up, up), dt, len(tau) - 1)
            flip[k] = probability(rho, down)
        ctx.count(model_evals=len(df) * len(tau))
    return {"detuning_mhz": df * 1e-6, "tau_ns": tau * 1e9, "flip_probability": flip}


def ramsey_fringes(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Ramsey fringes of each dot at a fixed drive detuning: the fringe frequency is the detuning, the
    Gaussian envelope is set by the quasi-static noise (hyperfine plus g-noise through dg/dV from the
    device's charge-noise draw) and the exponential one by the Markovian dephasing."""
    with ctx.compile():
        params, g, sigma_v = _spins(p)
        b = float(p.instance.get("b_t", 0.5))
        delta = H_EV_S * float(p.instance.get("f_detuning_hz", 2e6))
        sig = [quasi_static_sigma_ev(params, b, sigma_v) for _ in g]
    tau = _axis(p.sweep[0]) * 1e-6
    out: dict[str, np.ndarray] = {"tau_us": tau * 1e6}
    with ctx.run(len(tau) * len(g)):
        for k, s in enumerate(sig):
            out[f"p_down_dot{k}"] = ramsey_probability(tau, delta, s, params.t_phi_markov_s)
            out[f"t2_star_dot{k}_us"] = np.array(t2_star_s(s) * 1e6)
        out["envelope"] = 0.5 + 0.5 * np.exp(-((tau / t2_star_s(sig[0])) ** 2)) * np.exp(-tau / params.t_phi_markov_s)
        ctx.count(model_evals=len(tau) * len(g))
    return out


# --- Pauli spin blockade ---------------------------------------------------------------------------


def psb_readout_window(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """(0,2) occupation against detuning for a singlet and for a triplet, for two valley draws of the
    receiving dot: the singlet steps at the anticrossing, the triplet only once the detuning exceeds the
    valley splitting, and the contrast between them is the readout window."""
    inst = p.instance
    with ctx.compile():
        t = float(inst.get("t0_uev", 20.0)) * 1e-6
        kt = 8.617333262e-5 * float(inst.get("temperature_k", 0.15))
        radius = float(inst.get("radius_nm", 10.0))
        seeds = [int(s) for s in inst.get("valley_seeds", [1, 3])]
        params = ValleyParams.for_material(inst.get("material", "SiGe"))
        e_v = [float(draw_valley_splitting_uev(params, np.array([radius]), np.random.default_rng(s))[0]) * 1e-6 for s in seeds]
    eps = _axis(p.sweep[0]) * 1e-6
    out: dict[str, np.ndarray] = {"eps_uev": eps * 1e6}
    with ctx.run(len(eps) * len(e_v)):
        for k, ev in enumerate(e_v):
            s, tr = psb_occupation(eps, ev, t, kt)
            out["occupation_singlet"] = s
            out[f"occupation_triplet_draw{k}"] = tr
            out[f"contrast_draw{k}"] = psb_contrast(eps, ev, t, kt)
            out[f"valley_draw{k}_uev"] = np.array(ev * 1e6)
        ctx.count(model_evals=len(eps) * len(e_v))
    return out
