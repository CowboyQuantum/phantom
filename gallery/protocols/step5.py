"""Step 5 protocols: the valley-splitting distribution and its three signatures."""

from __future__ import annotations

import numpy as np

from gallery.runner import Context, Protocol
from phantom.valley import (
    ValleyParams,
    draw_valley_splitting_uev,
    kink_field_t,
    pulsed_spectroscopy_occupation,
    relaxation_rate_hz,
    second_electron_shift_uev,
)


def valley_splitting_draws(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Histograms of the valley splitting over many dots for Si/SiGe and Si MOS at one dot radius, the
    fraction of dots below a usable threshold, and the mean against dot radius for SiGe."""
    inst = p.instance
    with ctx.compile():
        n = int(inst.get("n_dots", 4000))
        radius = float(inst.get("radius_nm", 10.0))
        threshold = float(inst.get("threshold_uev", 50.0))
        radii = np.linspace(float(inst.get("radius_min_nm", 6.0)), float(inst.get("radius_max_nm", 30.0)), 25)
        sige, mos = ValleyParams.for_material("SiGe"), ValleyParams.for_material("MOS")
    with ctx.run(2 * n + 25 * n):
        d_sige = draw_valley_splitting_uev(sige, np.full(n, radius), np.random.default_rng(int(inst.get("seed", 0))))
        d_mos = draw_valley_splitting_uev(mos, np.full(n, radius), np.random.default_rng(int(inst.get("seed", 0)) + 1))
        means = np.array([draw_valley_splitting_uev(sige, np.full(n, r), np.random.default_rng(100 + k)).mean() for k, r in enumerate(radii)])
        ctx.count(model_evals=2 * n + 25 * n)
    edges = np.linspace(0, max(d_sige.max(), d_mos.max()), 61)
    centres = 0.5 * (edges[1:] + edges[:-1])
    h_sige, _ = np.histogram(d_sige, edges)
    h_mos, _ = np.histogram(d_mos, edges)
    return {
        "splitting_uev": centres,
        "count_sige": h_sige.astype(float),
        "count_mos": h_mos.astype(float),
        "fraction_below_threshold_sige": np.array((d_sige < threshold).mean()),
        "fraction_below_threshold_mos": np.array((d_mos < threshold).mean()),
        "radius_nm": radii,
        "mean_sige_uev": means,
        "inverse_radius_law_uev": means[0] * radii[0] / radii,
    }


def magnetospectroscopy_transition(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """The second-electron transition position against in-plane field for three valley draws: a kink
    at the field where the Zeeman energy equals the valley splitting."""
    inst = p.instance
    with ctx.compile():
        params = ValleyParams.for_material(inst.get("material", "SiGe"))
        radius = float(inst.get("radius_nm", 10.0))
        lever = float(inst.get("lever_arm_ev_per_v", 0.2))
        seeds = [int(s) for s in inst.get("seeds", [1, 2, 3])]
        draws = [float(draw_valley_splitting_uev(params, np.array([radius]), np.random.default_rng(s))[0]) for s in seeds]
    b = np.linspace(float(p.sweep[0]["from"]), float(p.sweep[0]["to"]), int(p.sweep[0]["points"]))
    out: dict[str, np.ndarray] = {"b_t": b}
    with ctx.run(len(b) * len(draws)):
        for k, e_v in enumerate(draws):
            out[f"shift_draw{k}_uv"] = second_electron_shift_uev(e_v, b) / lever  # µeV / (eV/V) = µV
            out[f"valley_draw{k}_uev"] = np.array(e_v)
            out[f"kink_draw{k}_t"] = np.array(kink_field_t(e_v))
        ctx.count(model_evals=len(b) * len(draws))
    return out


def t1_hot_spot(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Spin relaxation time against field for three valley draws: the hot spot dips where the Zeeman
    energy crosses each draw's valley splitting."""
    inst = p.instance
    with ctx.compile():
        params = ValleyParams.for_material(inst.get("material", "SiGe"))
        radius = float(inst.get("radius_nm", 10.0))
        seeds = [int(s) for s in inst.get("seeds", [1, 2, 3])]
        draws = [float(draw_valley_splitting_uev(params, np.array([radius]), np.random.default_rng(s))[0]) for s in seeds]
        hs = inst.get("hot_spot", {})
    b = np.linspace(float(p.sweep[0]["from"]), float(p.sweep[0]["to"]), int(p.sweep[0]["points"]))
    out: dict[str, np.ndarray] = {"b_t": b}
    with ctx.run(len(b) * len(draws)):
        for k, e_v in enumerate(draws):
            rate = relaxation_rate_hz(b, e_v, float(hs.get("background_hz", 1.0)), float(hs.get("peak_hz", 1e5)), float(hs.get("width_uev", 5.0)))
            out[f"t1_draw{k}_s"] = 1.0 / rate
            out[f"valley_draw{k}_uev"] = np.array(e_v)
        ctx.count(model_evals=len(b) * len(draws))
    return out


def pulsed_gate_spectroscopy(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Time-averaged occupation of the second electron against plunger and pulse amplitude for one
    valley draw: the transition splits by the amplitude, and the excited-state line runs E_V / α away."""
    inst = p.instance
    with ctx.compile():
        params = ValleyParams.for_material(inst.get("material", "SiGe"))
        radius = float(inst.get("radius_nm", 10.0))
        e_v = float(draw_valley_splitting_uev(params, np.array([radius]), np.random.default_rng(int(inst.get("seed", 1))))[0])
        lever = float(inst.get("lever_arm_ev_per_v", 0.2))
    v = np.linspace(float(p.sweep[0]["from"]), float(p.sweep[0]["to"]), int(p.sweep[0]["points"]))
    a = np.linspace(float(p.sweep[1]["from"]), float(p.sweep[1]["to"]), int(p.sweep[1]["points"]))
    with ctx.run(len(v) * len(a)):
        n = pulsed_spectroscopy_occupation(v, a, 0.0, e_v, lever, float(inst.get("kt_uev", 13.0)), float(inst.get("excited_weight", 0.4)))
        ctx.count(model_evals=len(v) * len(a))
    return {"v_plunger_mv": v * 1e3, "pulse_amplitude_mv": a * 1e3, "occupation": n.T, "valley_uev": np.array(e_v), "excited_line_offset_mv": np.array(e_v * 1e-6 / lever * 1e3)}
