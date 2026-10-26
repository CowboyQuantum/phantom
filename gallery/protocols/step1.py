"""Step 1 protocols: single dot and sensor, double-dot stability, noise. Each function reports its
model evaluations into `ctx` and returns the arrays named by the protocol's `render` block."""

from __future__ import annotations

import numpy as np

from gallery.runner import Context, Protocol
from phantom.charge import Barrier
from phantom.readout import (
    LeadRateLandscape,
    LockinParams,
    SensedDoubleDot,
    SensorSETParams,
    lockin_sweep_dynamic_a_per_v,
    raster_times_s,
)
from phantom.readout.transport import ChannelTurnOn, pinch_off_current_a


def _device(p: Protocol) -> SensedDoubleDot:
    inst = p.instance
    if inst.get("noise", "default") == "none":
        return SensedDoubleDot()
    return SensedDoubleDot.default_noisy(seed=int(inst.get("seed", 0)))


def _axis(spec: dict) -> np.ndarray:
    return np.linspace(float(spec["from"]), float(spec["to"]), int(spec["points"]))


def _grid3(p1: np.ndarray, p2: np.ndarray, sg: float = 0.0) -> np.ndarray:
    v = np.zeros((len(p1), len(p2), 3))
    v[..., 0] = p1[:, None]
    v[..., 1] = p2[None, :]
    v[..., 2] = sg
    return v


# --- A/B: electrical health shapes and the single dot -------------------------------------------


def turn_on_single_gate(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    with ctx.compile():
        channel = ChannelTurnOn(**p.instance.get("channel", {}))
    v = _axis(p.sweep[0])
    with ctx.run(len(v)):
        i = channel.current_a(v)
        ctx.count(model_evals=len(v))
    return {"v_gate_v": v, "current_a": i}


def pinch_off_barrier(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    with ctx.compile():
        barrier = Barrier("B", **p.instance.get("barrier", {}))
    v = _axis(p.sweep[0])
    with ctx.run(len(v)):
        i = pinch_off_current_a(barrier, v, float(p.instance.get("i_sat_a", 2e-9)))
        ctx.count(model_evals=len(v))
    return {"v_barrier_v": v, "current_a": i}


def coulomb_oscillations_single_dot(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    with ctx.compile():
        d = _device(p)
    sg = _axis(p.sweep[0])
    v = np.zeros((len(sg), 3))
    v[:, 2] = sg
    with ctx.run(len(sg)):
        i = d.sweep_a(v, raster_times_s(v.shape[:-1], 1e-3) if p.instance.get("noise") != "none" else None)
        ctx.count(model_evals=len(sg))
    return {"v_sg_v": sg, "current_a": i}


def diamonds_single_dot(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    with ctx.compile():
        d = _device(p)
    sg, vsd = _axis(p.sweep[0]), _axis(p.sweep[1])
    v = np.zeros((len(sg), 3))
    v[:, 2] = sg
    with ctx.run(len(sg) * len(vsd)):
        n = d.charge.mean_occupation(v)
        phi = d.sensor.potential(v, n)
        i = d.sensor.current_a(phi[:, None], vsd[None, :])
        ctx.count(model_evals=len(sg) * len(vsd))
    return {"v_sg_v": sg, "v_sd_v": vsd, "current_a": i}


def sensor_coulomb_peak(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    with ctx.compile():
        d = _device(p)
        centre = d.sensor.phase_v
    half = float(p.sweep[0].get("half_span_v", 0.15))
    sg = np.linspace(centre - half, centre + half, int(p.sweep[0]["points"]))
    v = np.zeros((len(sg), 3))
    v[:, 2] = sg
    with ctx.run(len(sg)):
        i = d.current_a(v)
        ctx.count(model_evals=len(sg))
    slope = np.gradient(i, sg)
    return {"v_sg_v": sg, "current_a": i, "slope_a_per_v": slope}


# --- C: double-dot stability ---------------------------------------------------------------------


def stability_coarse_2dot(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    with ctx.compile():
        d = _device(p)
    p1, p2 = _axis(p.sweep[0]), _axis(p.sweep[1])
    v = _grid3(p1, p2, float(p.instance.get("sg_v", 0.0)))
    with ctx.run(v.shape[0] * v.shape[1]):
        times = None if p.instance.get("noise") == "none" else raster_times_s(v.shape[:-1], float(p.readout.get("dt_s", 1e-3)))
        i = d.sweep_a(v, times)
        ctx.count(model_evals=v.shape[0] * v.shape[1])
    return {"v_p1_v": p1, "v_p2_v": p2, "current_a": i}


def _fine_window(d: SensedDoubleDot, half_v: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    tp1, tp2 = d.charge.triple_points_ng()
    v1 = d.charge.gate_voltages_for_ng(tp1)
    v2 = d.charge.gate_voltages_for_ng(tp2)
    c = (v1 + v2) / 2
    return np.linspace(c[0] - half_v, c[0] + half_v, n), np.linspace(c[1] - half_v, c[1] + half_v, n)


def stability_fine_2dot(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    with ctx.compile():
        d = _device(p)
        p1, p2 = _fine_window(d, float(p.sweep[0]["half_span_v"]), int(p.sweep[0]["points"]))
    v = _grid3(p1, p2)
    with ctx.run(v.shape[0] * v.shape[1]):
        times = None if p.instance.get("noise") == "none" else raster_times_s(v.shape[:-1], float(p.readout.get("dt_s", 1e-3)))
        i = d.sweep_a(v, times)
        ctx.count(model_evals=v.shape[0] * v.shape[1])
    return {"v_p1_v": p1, "v_p2_v": p2, "current_a": i}


def stability_lockin(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Chopped-gate stability map with tunnel-rate dynamics. The chop direction (P1, P2, common mode
    or detuning) comes from the protocol's readout.lockin block; the lead-rate landscape from its
    instance block. Scan order: P2 slow, P1 fast."""
    inst = p.instance
    with ctx.compile():
        d = _device(p)
        lk = LockinParams(**p.readout["lockin"])
        land = inst.get("lead_rates", {})
        rates = LeadRateLandscape(
            d.charge.gate_names,
            rate0_hz=np.array(land.get("rate0_hz", [1e3, 1e3]), dtype=float),
            weights_per_v=np.array(land.get("weights_per_v", [[125.0, 0.0, 0.0], [0.0, 125.0, 0.0]]), dtype=float),
        )
        rng = np.random.default_rng(int(inst.get("seed", 0)) + 101)
    p1, p2 = _axis(p.sweep[0]), _axis(p.sweep[1])
    v = np.zeros((len(p2), len(p1), 3))
    v[..., 0] = p1[None, :]
    v[..., 1] = p2[:, None]
    with ctx.run(v.shape[0] * v.shape[1]):
        t0 = None if inst.get("noise") == "none" else 0.0
        out = lockin_sweep_dynamic_a_per_v(d, v, lk, rates, rng, t0_s=t0)
        ctx.count(model_evals=v.shape[0] * v.shape[1] * lk.evaluations_per_point)
    return {
        "v_p1_v": p1,
        "v_p2_v": p2,
        "lockin_a_per_v": out["signal_a_per_v"].T,
        "log10_rate_hz": out["log10_rate_hz"].T,
        "n1_present": out["n_present_mean"][..., 0].T,
    }


def plunger_vs_barrier_map(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """P1 (fast) against the left lead barrier B_L (slow) on the layout: for every barrier setting the
    lead rate is the WKB rate of the landscape and the charge follows the master equation along the
    row at the dwell time, carried across the flyback. Dot 1 is the only dot whose lead is swept;
    where the rate falls below the inverse dwell the lines fade, latch and carry across rows."""
    from gallery.protocols.step2 import _layout, _operating_point
    from gallery.protocols.step4 import derive_dynamics
    from phantom.electrostatics import GateBasis

    inst = p.instance
    with ctx.compile():
        lay = _layout(p)
        v0 = _operating_point(lay, p)
        basis = GateBasis.from_layout(lay)
        j_p1, j_bl = lay.gate_names.index("P1"), lay.gate_names.index("B_L")
        temperature = float(inst.get("temperature_k", 0.15))
        dwell = float(p.readout.get("dt_s", 1e-3))
        p1 = v0[j_p1] + np.linspace(float(p.sweep[0]["from"]), float(p.sweep[0]["to"]), int(p.sweep[0]["points"]))
        vbl = np.linspace(float(p.sweep[1]["from"]), float(p.sweep[1]["to"]), int(p.sweep[1]["points"]))
        rows = []
        for vb in vbl:
            v = v0.copy()
            v[j_bl] = vb
            d, rate_eq, _ = derive_dynamics(lay, v, basis, temperature)
            rows.append((v, d, rate_eq))
    n1, n2 = len(p1), len(vbl)
    current = np.empty((n2, n1))
    n1_present = np.empty((n2, n1))
    n1_eq = np.empty((n2, n1))
    gamma = np.empty(n2)
    with ctx.run(n1 * n2):
        p_prev = states_prev = None
        for r, (v, d, rate_eq) in enumerate(rows):
            vv = np.repeat(v[None, :], n1, axis=0)
            vv[:, j_p1] = p1
            mean, p_prev, states_prev = rate_eq.sweep(vv, dwell, p0=p_prev, states0=states_prev)
            sensor = SensorSETParams(alpha_v_per_v=d.sensor_alpha_v_per_v)
            current[r] = sensor.current_a(sensor.potential(vv, mean))
            n1_present[r] = mean[:, 0]
            n1_eq[r] = d.ci.mean_occupation(vv)[:, 0]
            gamma[r] = rate_eq.gamma_left_hz
        ctx.count(model_evals=n1 * n2)
    return {
        "v_p1_v": p1,
        "v_bl_v": vbl,
        "current_a": current.T,
        "n1_measured": n1_present.T,
        "n1_equilibrium": n1_eq.T,
        "log10_gamma_left_hz": np.log10(np.maximum(gamma, 1e-3)),
    }


# --- D: noise --------------------------------------------------------------------------------------


def _flank_point(d: SensedDoubleDot) -> np.ndarray:
    """A gate point on the steep flank of a sensor peak, plungers at 0 V."""
    return np.array([0.0, 0.0, d.sensor.phase_v + 0.9 * d.sensor.width_v])


def _trace(ctx: Context, p: Protocol) -> tuple[np.ndarray, np.ndarray]:
    with ctx.compile():
        d = _device(p)
        point = _flank_point(d)
    dt, duration = float(p.readout["dt_s"]), float(p.readout["duration_s"])
    t = np.arange(0.0, duration, dt)
    with ctx.run(len(t)):
        i = d.trace_a(point, t)
        ctx.count(model_evals=len(t))
    return t, i


def noise_time_trace(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    t, i = _trace(ctx, p)
    zoom = t < float(p.render.get("zoom_s", 2.0))
    return {"t_s": t, "current_a": i, "t_zoom_s": t[zoom], "current_zoom_a": i[zoom]}


def noise_psd(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    _, i = _trace(ctx, p)
    dt = float(p.readout["dt_s"])
    seg = int(p.readout.get("segment", 4096))
    n_seg = len(i) // seg
    x = (i[: n_seg * seg] - i.mean()).reshape(n_seg, seg)
    w = np.hanning(seg)
    scale = 1.0 / (np.sum(w**2) / dt)
    psd = (np.abs(np.fft.rfft(x * w, axis=1)) ** 2).mean(axis=0) * scale * 2
    f = np.fft.rfftfreq(seg, dt)
    keep = f > 0
    return {"f_hz": f[keep], "psd_a2_per_hz": psd[keep]}
