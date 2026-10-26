"""Step 6 protocols, second set: bias triangles with Pauli blockade, the spin funnel, Elzerman readout
traces, charge-noise spectroscopy by CPMG and by Ramsey against J, and the exchange-gate fidelity."""

from __future__ import annotations

import numpy as np

from gallery.protocols.step6 import GHZ_PER_EV, _axis, _pair
from gallery.runner import Context, Protocol
from phantom.charge import ConstantInteractionParams
from phantom.constants import HBAR_EV_S
from phantom.disorder import NoiseParams, draw_noise
from phantom.spin import (
    MU_B_EV_PER_T,
    BlockadeRates,
    SpinParams,
    bias_triangle_current_a,
    coherence,
    detuning_noise_sigma_ev,
    detuning_voltages,
    draw_spin,
    elzerman_shots,
    empty_probability,
    exchange_j_ev,
    exchange_sensitivity_ev_per_ev,
    filter_frequency_hz,
    funnel_singlet_probability,
    log_grid_hz,
    omega_spectrum,
    one_over_f_amplitude_sq,
    one_over_f_psd,
    quasi_static_average,
    singlet_return_probability,
    swap_gate_fidelity,
    t2_star_s,
)
from phantom.valley import ValleyParams, draw_valley_splitting_uev


def psb_bias_triangles(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Current against both plungers around the (1,1)/(0,2) vertex at finite bias: forward bias in a
    field (blocked base of the triangles up to the valley splitting), forward bias at zero field (the
    blockade lifted by mixing), and reverse bias (no blockade)."""
    inst = p.instance
    with ctx.compile():
        ci = ConstantInteractionParams()
        rates = BlockadeRates(**inst.get("rates", {}))
        v_sd = float(inst.get("v_sd_v", 1e-3))
        b = float(inst.get("b_t", 0.1))
        params = ValleyParams.for_material(inst.get("material", "SiGe"))
        e_v = float(draw_valley_splitting_uev(params, np.array([float(inst.get("radius_nm", 10.0))]), np.random.default_rng(int(inst.get("valley_seed", 3))))[0]) * 1e-6
        tp1, tp2 = ci.triple_points_ng()
        c = ci.gate_voltages_for_ng((tp1 + tp2) / 2 + np.array([0.0, 1.0]))
        half, n = float(p.sweep[0]["half_span_v"]), int(p.sweep[0]["points"])
    p1 = np.linspace(c[0] - half, c[0] + half, n)
    p2 = np.linspace(c[1] - half, c[1] + half, n)
    v = np.zeros((n, n, 3))
    v[..., 0], v[..., 1] = p1[:, None], p2[None, :]
    with ctx.run(n * n):
        i_fwd = bias_triangle_current_a(ci, v, v_sd, e_v, b, rates)
        i_zero = bias_triangle_current_a(ci, v, v_sd, e_v, 0.0, rates)
        i_rev = bias_triangle_current_a(ci, v, -v_sd, e_v, b, rates)
        ctx.count(model_evals=3 * n * n)
    return {"v_p1_v": p1, "v_p2_v": p2, "current_forward_pa": i_fwd * 1e12, "current_forward_zero_field_pa": i_zero * 1e12, "current_reverse_pa": i_rev * 1e12, "valley_uev": np.array(e_v * 1e6)}


def spin_funnel(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Singlet probability after separating a singlet into (1,1) against detuning and field, with a
    jittered wait: the S–T₋ funnel along J(ε) = ḡ μ_B B and the S–T₀ loss where the Zeeman difference
    of the g draw is comparable to J."""
    inst = p.instance
    with ctx.compile():
        hub, v0, _ = _pair(p)
        eac = float(hub.onsite_u_ev[1] - hub.mutual_v_ev)
        params = SpinParams(**inst.get("spin", {}))
        seed = int(inst.get("seed", 0))
        g = draw_spin(params, np.random.default_rng(seed), seed).g_factors
        delta_mix = float(inst.get("delta_mix_nev", 20.0)) * 1e-9
    eps = _axis(p.sweep[0]) * eac
    b = _axis(p.sweep[1])
    with ctx.run(len(eps) * len(b)):
        j = exchange_j_ev(hub, detuning_voltages(hub, v0, eps))
        ez = 0.5 * (g[0] + g[1]) * MU_B_EV_PER_T * b
        dz = (g[0] - g[1]) * MU_B_EV_PER_T * b
        p_s = funnel_singlet_probability(j[:, None], ez[None, :], delta_mix, dz[None, :])
        ctx.count(model_evals=len(eps) * len(b))
    return {"eps_mev": eps * 1e3, "b_mt": b * 1e3, "singlet_probability": p_s, "j_ghz": j * GHZ_PER_EV, "funnel_b_mt": j / (0.5 * (g[0] + g[1]) * MU_B_EV_PER_T) * 1e3}


def elzerman_trace(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Charge traces after pulsing to the read level for eight shots (offset for display), and the
    ensemble probability that the dot is empty for an upper-state electron, sampled over many shots and
    analytic: the read window a consumer thresholds."""
    inst = p.instance
    with ctx.compile():
        g_out, g_in, t1 = float(inst.get("gamma_out_hz", 2e4)), float(inst.get("gamma_in_hz", 1e4)), float(inst.get("t1_s", 1e-2))
        dt, n_t = float(p.sweep[0]["to"]) * 1e-6 / (int(p.sweep[0]["points"]) - 1), int(p.sweep[0]["points"])
        n_shots, n_show = int(inst.get("n_shots", 2000)), int(inst.get("n_show", 8))
        seed = int(inst.get("seed", 0))
    t = np.arange(n_t) * dt
    with ctx.run(n_shots * n_t):
        traces, spins = elzerman_shots(np.random.default_rng(seed), n_shots, float(inst.get("p_upper", 1.0)), g_out, g_in, t1, dt, n_t)
        ctx.count(model_evals=n_shots * n_t)
    out: dict[str, np.ndarray] = {"t_us": t * 1e6, "empty_sampled": 1 - traces.mean(axis=0), "empty_analytic": empty_probability(t, g_out, g_in, t1)}
    for k in range(n_show):
        out[f"shot{k}"] = traces[k].astype(float) + 1.2 * k
    out["spins_shown"] = spins[:n_show].astype(float)
    return out


def _qubit_spectrum(p: Protocol, sensitivity_ev_per_v: float):
    inst = p.instance
    seed = int(inst.get("seed", 0))
    noise = draw_noise(NoiseParams.for_material(inst.get("material", "SiGe")), np.random.default_rng(seed), seed)
    a2 = one_over_f_amplitude_sq(noise.one_over_f_amplitude_v, noise.f_min_hz, noise.f_max_hz)
    f_high = float(inst.get("f_high_hz", 1e9))
    f = log_grid_hz(noise.f_min_hz, f_high)
    s_v = one_over_f_psd(f, a2, noise.f_min_hz, f_high)
    return f, s_v, omega_spectrum(s_v, sensitivity_ev_per_v), noise


def charge_noise_spectroscopy_cpmg(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Coherence of the exchange qubit at an operating detuning after N CPMG π pulses in a fixed total
    time, on the device's 1/f charge-noise spectrum (the same amplitude as the sensor's `noise-psd`,
    extended to the qubit's band), plotted against the filter frequency N/2τ; and the detuning and
    gate-referred spectra themselves."""
    inst = p.instance
    with ctx.compile():
        hub, v0, _ = _pair(p)
        eac = float(hub.onsite_u_ev[1] - hub.mutual_v_ev)
        eps_op = float(inst.get("eps_over_anticrossing", 0.6)) * eac
        dj = float(exchange_sensitivity_ev_per_ev(hub, v0, np.array([eps_op]))[0])
        sens_v = dj * detuning_noise_sigma_ev(hub, 1.0)  # eV per volt of gate noise, through dJ/dε
        f, s_v, s_omega, noise = _qubit_spectrum(p, sens_v)
        tau = float(inst.get("tau_us", 10.0)) * 1e-6
        n_list = [int(n) for n in inst.get("n_pulses", [0, 1, 2, 4, 8, 16, 32, 64, 128])]
    with ctx.run(len(n_list)):
        w = np.array([coherence(tau, n, f, s_omega) for n in n_list])
        ctx.count(model_evals=len(n_list))
    f_n = np.array([filter_frequency_hz(tau, n) if n > 0 else np.nan for n in n_list])
    s_eps = s_v * detuning_noise_sigma_ev(hub, 1.0) ** 2
    return {
        "n_pulses": np.array(n_list, dtype=float),
        "coherence": w,
        "filter_frequency_hz": f_n,
        "f_hz": f,
        "s_v_v2_per_hz": s_v,
        "s_eps_uev2_per_hz": s_eps * 1e12,
        "dj_deps": np.array(dj),
        "sigma_eps_uev": np.array(detuning_noise_sigma_ev(hub, noise.one_over_f_amplitude_v) * 1e6),
    }


def charge_noise_spectroscopy_j_vs_ramsey(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Ramsey T₂* of the exchange qubit against detuning: analytic √2ħ/(σ_ε dJ/dε) and the 1/e time of
    the quadrature-averaged singlet return, with 1/(dJ/dε) overlaid; the direct check that the same
    charge-noise draw governs the qubit and the sensor."""
    inst = p.instance
    with ctx.compile():
        hub, v0, sigma_v = _pair(p)
        eac = float(hub.onsite_u_ev[1] - hub.mutual_v_ev)
        sigma_eps = detuning_noise_sigma_ev(hub, sigma_v)
        nodes = int(inst.get("quadrature_nodes", 31))
    eps = _axis(p.sweep[0]) * eac
    with ctx.run(len(eps)):
        dj = exchange_sensitivity_ev_per_ev(hub, v0, eps)
        t2_analytic = np.array([t2_star_s(abs(d) * sigma_eps) for d in dj])
        j0 = exchange_j_ev(hub, detuning_voltages(hub, v0, eps))
        t2_num = np.empty(len(eps))
        for k, e in enumerate(eps):
            t_max = 3 * min(t2_analytic[k], 1e-3)
            period = 2 * np.pi * HBAR_EV_S / max(j0[k], 1e-12)
            tau = np.linspace(0, t_max, int(np.clip(20 * t_max / period, 400, 20000)))

            def f(delta, e=e, tau=tau):
                j = exchange_j_ev(hub, detuning_voltages(hub, v0, e + delta))
                return singlet_return_probability(tau[:, None], j[None, :])

            p_s = quasi_static_average(f, sigma_eps, nodes)
            # the envelope of the oscillation: the running maximum of |2 P_S − 1| sampled per period
            env = np.abs(2 * p_s - 1)
            win = max(int(period / (tau[1] - tau[0])), 1)
            env_max = np.array([env[i : i + win].max() for i in range(0, len(env), win)])
            t_env = tau[::win][: len(env_max)]
            below = np.flatnonzero(env_max < np.exp(-1))
            t2_num[k] = t_env[below[0]] if len(below) else np.nan
        ctx.count(model_evals=len(eps))
    return {"eps_mev": eps * 1e3, "t2_star_analytic_us": t2_analytic * 1e6, "t2_star_from_decay_us": t2_num * 1e6, "inverse_sensitivity_scaled_us": np.sqrt(2) * HBAR_EV_S / (np.abs(dj) * sigma_eps) * 1e6, "dj_deps": dj, "sigma_eps_uev": np.array(sigma_eps * 1e6)}


def exchange_gate_fidelity(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Average fidelity of the SWAP gate (J on for πħ/J) against operating detuning for the noiseless
    device and for two noise draws, with the Markovian T₁ and T_φ of the parameters: the sweet spot at
    the symmetric point and the loss toward the anticrossing where dJ/dε grows; the gate time alongside."""
    inst = p.instance
    with ctx.compile():
        hub, v0, _ = _pair(p)
        eac = float(hub.onsite_u_ev[1] - hub.mutual_v_ev)
        params = SpinParams(**inst.get("spin", {}))
        seeds = [int(s) for s in inst.get("noise_seeds", [0, 1])]
        sigmas = [detuning_noise_sigma_ev(hub, float(draw_noise(NoiseParams.for_material(inst.get("material", "SiGe")), np.random.default_rng(s), s).one_over_f_amplitude_v)) for s in seeds]
        nodes = int(inst.get("quadrature_nodes", 15))
    eps = _axis(p.sweep[0]) * eac
    out: dict[str, np.ndarray] = {"eps_mev": eps * 1e3}
    with ctx.run(len(eps)):
        clean = np.empty(len(eps))
        t_gate = np.empty(len(eps))
        for k, e in enumerate(eps):
            clean[k], t_gate[k] = swap_gate_fidelity(hub, v0, float(e), 0.0, params, n_nodes=nodes)
        out["fidelity_noiseless"] = clean
        for s, sig in enumerate(sigmas):
            out[f"fidelity_draw{s}"] = np.array([swap_gate_fidelity(hub, v0, float(e), sig, params, n_nodes=nodes)[0] for e in eps])
            out[f"sigma_eps_draw{s}_uev"] = np.array(sig * 1e6)
        ctx.count(model_evals=len(eps))
    out["infidelity_draw0"] = 1 - out["fidelity_draw0"]
    out["infidelity_draw1"] = 1 - out["fidelity_draw1"]
    out["infidelity_noiseless"] = 1 - clean
    out["gate_time_ns"] = t_gate * 1e9
    return out
