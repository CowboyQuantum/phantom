"""Step 7 protocols, second set: the optimum's coupling target, read point and sweet spot checked
against the compile tier and against scans, and a population of compiled devices across seeds."""

from __future__ import annotations

import numpy as np

from gallery.protocols.step2 import _layout
from gallery.protocols.step7 import _compiled
from gallery.runner import Context, Protocol
from phantom.charge import derive_double_dot
from phantom.compile import CompileParams, compile_device, interdot_tunnel_coupling_ev
from phantom.constants import K_B_EV_PER_K
from phantom.layout import MaterialClass
from phantom.optimum import coherence_optimum, coupling_optimum, readout_optimum, sensor_optimum
from phantom.spin import SpinParams, empty_probability, psb_contrast, swap_gate_fidelity
from phantom.spin.exchange import detuning_noise_sigma_ev, detuning_voltages, exchange_j_ev

GHZ_PER_EV = 1.0 / 4.135667696e-15 * 1e-9


def _marker(x: np.ndarray, x0: float, height: float) -> np.ndarray:
    """A spike of `height` at the grid point nearest `x0`, zero elsewhere."""
    out = np.zeros_like(x, dtype=float)
    out[int(np.argmin(np.abs(x - x0)))] = height
    return out


def optimum_coupling_target(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """The interdot tunnel coupling against the barrier gate: recomputed at the compile tier (WKB action
    along the barrier of the landscape at every barrier voltage) and from the compiled exponential the
    runtime uses. For three target couplings, the barrier voltage the optimum chooses (from the
    exponential) beside the voltage where the compile-tier curve crosses the target."""
    inst = p.instance
    with ctx.compile():
        dev = _compiled(p)
        lay, basis = dev.basis()
        v0 = dev.operating_point()
        derived = derive_double_dot(lay, v0, temperature_k=dev.temperature_k, basis=basis)
        d1, d2 = derived.dots
        mass = MaterialClass.by_name(dev.material).effective_mass_me
        jb = dev.gate_names.index(dev.barrier_gate)
        vb = v0[jb] + np.linspace(float(p.sweep[0]["from"]), float(p.sweep[0]["to"]), int(p.sweep[0]["points"]))
        t_compile = np.empty(len(vb))
        for k, x in enumerate(vb):
            v = v0.copy()
            v[jb] = x
            t_compile[k] = interdot_tunnel_coupling_ev(basis, v, d1, d2, mass)[0]
        targets = np.array([float(t) for t in inst.get("targets_uev", [10.0, 30.0, 100.0])]) * 1e-6
    with ctx.run(len(vb) + len(targets)):
        barrier = dev.runtime_barrier()
        t_runtime = barrier.tunnel_coupling_ev(vb)
        chosen = np.array([coupling_optimum(dev, float(t)).v_barrier_v for t in targets])
        ctx.count(model_evals=len(vb) + len(targets))
    # where the compile-tier curve crosses each target, by log-linear interpolation on its rising part
    # (it flattens where the barrier disappears and the dots merge)
    log_t = np.log(np.maximum(t_compile, 1e-300))
    rising = np.concatenate([[True], np.cumprod(np.diff(log_t) > 1e-9).astype(bool)])
    brute = np.interp(np.log(targets), log_t[rising], vb[rising], left=np.nan, right=np.nan)
    return {
        "v_barrier_v": vb,
        "t_compile_uev": t_compile * 1e6,
        "t_runtime_uev": t_runtime * 1e6,
        "target_uev": targets * 1e6,
        "v_chosen_by_optimum_v": chosen,
        "v_from_compile_tier_v": brute,
        "operating_barrier_v": np.array(v0[jb]),
    }


def optimum_readout_point(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Pauli-blockade contrast against detuning for the compiled device (the receiving dot's valley
    draw, the compiled tunnel coupling, the electron temperature) with the read point the optimum
    chooses marked; and the Elzerman empty probability against time for the compiled lead rate with the
    optimum's read time marked."""
    with ctx.compile():
        dev = _compiled(p)
        kt = K_B_EV_PER_K * dev.temperature_k
        e_v = dev.dots[1].valley_splitting_uev * 1e-6
        eps = np.linspace(float(p.sweep[0]["from"]), float(p.sweep[0]["to"]), int(p.sweep[0]["points"])) * 1e-6
        t1 = float(p.instance.get("t1_read_s", 1e-3))
        g = dev.barriers[0].rate_hz
        t = np.linspace(0.0, 6.0 / g, int(p.sweep[1]["points"]))
    with ctx.run(len(eps) + len(t)):
        r = readout_optimum(dev, t1)
        c = psb_contrast(eps, e_v, dev.tunnel_coupling_ev, kt)
        pe = empty_probability(t, g, g, t1)
        ctx.count(model_evals=len(eps) + len(t))
    return {
        "eps_uev": eps * 1e6,
        "psb_contrast": c,
        "psb_chosen_marker": _marker(eps, r.psb_detuning_ev, r.psb_contrast),
        "t_ns": t * 1e9,
        "p_empty": pe,
        "elzerman_chosen_marker": _marker(t, r.elzerman_peak_time_s, r.elzerman_peak_probability),
        "valley_receiving_uev": np.array(e_v * 1e6),
        "tunnel_coupling_uev": np.array(dev.tunnel_coupling_ev * 1e6),
        "lead_rate_hz": np.array(g),
    }


def optimum_sweet_spot(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """SWAP-gate infidelity and the exchange energy against detuning through the sweet spot the optimum
    chooses, for the compiled device with its noise draw: the infidelity is lowest where dJ/dε
    vanishes, at the chosen point."""
    with ctx.compile():
        dev = _compiled(p)
        hub = dev.runtime_hubbard()
        spin = SpinParams(**dev.params["spin"])
        c = coherence_optimum(dev)
        v_sweet = np.array(c.sweet_spot_v)
        sigma = detuning_noise_sigma_ev(hub, dev.noise_sigma_gate_v())
        eac = float(hub.onsite_u_ev[1] - hub.mutual_v_ev)
        frac = np.linspace(float(p.sweep[0]["from"]), float(p.sweep[0]["to"]), int(p.sweep[0]["points"]))
        eps = frac * eac
    with ctx.run(len(eps)):
        f = np.array([swap_gate_fidelity(hub, v_sweet, float(e), sigma, spin)[0] for e in eps])
        j = exchange_j_ev(hub, detuning_voltages(hub, v_sweet, eps))
        ctx.count(model_evals=len(eps))
    infid = np.maximum(1.0 - f, 1e-12)
    return {
        "eps_mev": eps * 1e3,
        "swap_infidelity": infid,
        "chosen_marker": _marker(eps, 0.0, float(infid.max())),
        "j_ghz": j * GHZ_PER_EV,
        "fidelity_at_chosen": np.array(c.swap_fidelity_at_sweet_spot),
        "sigma_eps_uev": np.array(sigma * 1e6),
    }


def compiled_device_draws(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Twelve devices compiled from the same layout and operating point with seeds 0–11: the interdot
    coupling and the lead rates, the valley splittings, and what the optimum reaches for each (sensor
    slope, blockade contrast, SWAP infidelity at the sweet spot of the compiled operating point). A
    device whose disorder stops a dot or the sensor from forming is recorded as not formed (NaN) — a
    yield, not an error."""
    inst = p.instance
    seeds = [int(s) for s in inst.get("seeds", list(range(12)))]
    with ctx.compile():
        lay = _layout(p)
        params = CompileParams(operating_point_v={k: float(v) for k, v in inst["operating_point_v"].items()}, temperature_k=float(inst.get("temperature_k", 0.15)))
        devices = []
        for s in seeds:
            try:
                devices.append(compile_device(lay, params, s))
            except ValueError:
                devices.append(None)
    n = len(seeds)
    nan = np.full(n, np.nan)
    t_uev, g_l, g_r, ev1, ev2, dg, slope, contrast, infid = (nan.copy() for _ in range(9))
    with ctx.run(n):
        for k, dev in enumerate(devices):
            if dev is None:
                continue
            t_uev[k] = dev.tunnel_coupling_ev * 1e6
            g_l[k], g_r[k] = dev.barriers[0].rate_hz, dev.barriers[1].rate_hz
            ev1[k], ev2[k] = dev.valley_splittings_uev()[:2]
            dg[k] = abs(dev.g_factors[0] - dev.g_factors[1])
            slope[k] = sensor_optimum(dev).slope_a_per_v
            contrast[k] = readout_optimum(dev).psb_contrast
            infid[k] = max(1.0 - coherence_optimum(dev).swap_fidelity_at_sweet_spot, 1e-12)
        ctx.count(model_evals=n)
    return {
        "seed": np.array(seeds, dtype=float),
        "formed": np.array([d is not None for d in devices], dtype=float),
        "log10_t_uev": np.log10(t_uev),
        "log10_gamma_left_hz": np.log10(g_l),
        "log10_gamma_right_hz": np.log10(g_r),
        "valley_dot1_uev": ev1,
        "valley_dot2_uev": ev2,
        "g_difference": dg,
        "sensor_slope_na_per_v": slope * 1e9,
        "psb_contrast": contrast,
        "log10_swap_infidelity": np.log10(infid),
    }
