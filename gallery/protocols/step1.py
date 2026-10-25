"""Step 1 protocols: single dot and sensor, double-dot stability, noise. Each function reports its
model evaluations into `ctx` and returns the arrays named by the protocol's `render` block."""

from __future__ import annotations

import numpy as np

from gallery.runner import Context, Protocol
from phantom.charge import Barrier, ConstantInteractionParams
from phantom.readout import (
    LeadBarrier,
    LockinParams,
    SensedDoubleDot,
    SensorSETParams,
    latched_occupation,
    lockin_sweep_a_per_v,
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


def stability_lockin_2dot(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    with ctx.compile():
        d = _device(p)
        lk = LockinParams(**p.readout["lockin"])
    p1, p2 = _axis(p.sweep[0]), _axis(p.sweep[1])
    v = _grid3(p1, p2)
    with ctx.run(v.shape[0] * v.shape[1]):
        t0 = None if p.instance.get("noise") == "none" else 0.0
        g = lockin_sweep_a_per_v(d, v, lk, t0_s=t0)
        ctx.count(model_evals=v.shape[0] * v.shape[1] * lk.evaluations_per_point)
    return {"v_p1_v": p1, "v_p2_v": p2, "lockin_a_per_v": g}


def plunger_vs_barrier_map(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Dot 1 loads from its reservoir through the lead barrier BL; dot 2 has no reservoir of its own.
    P1 is the fast axis (one row per barrier value, flyback between rows); BL is the slow axis."""
    inst = p.instance
    with ctx.compile():
        base = SensedDoubleDot()
        bl = inst.get("lead_barrier_lever_arms_af", [3.0, 0.5])
        c_gate = np.concatenate([base.charge.c_gate_af, np.array(bl, dtype=float)[:, None]], axis=1)
        charge = ConstantInteractionParams(
            gate_names=base.charge.gate_names + ("BL",),
            c_gate_af=c_gate,
            temperature_k=base.charge.temperature_k,
        )
        alpha = np.concatenate([base.sensor.alpha_v_per_v, [float(inst.get("sensor_alpha_bl", 0.005))]])
        sensor = SensorSETParams(**{**base.sensor.__dict__, "alpha_v_per_v": alpha})
        lead = LeadBarrier("BL", **inst.get("lead_barrier", {}))
        p2_v = float(inst.get("p2_ng", 0.9)) * charge.transition_period_v(1, "P2")
        rng = np.random.default_rng(int(inst.get("seed", 0)))
    p1 = _axis(p.sweep[0])
    vbl = _axis(p.sweep[1])
    dwell = float(p.readout.get("dt_s", 1e-3))
    # scan order: for each barrier value (slow), sweep P1 (fast)
    v = np.zeros((len(vbl), len(p1), 4))
    v[..., 0] = p1[None, :]
    v[..., 1] = p2_v
    v[..., 3] = vbl[:, None]
    flat = v.reshape(-1, 4)
    with ctx.run(len(p1) * len(vbl)):
        n_eq = charge.mean_occupation(flat)
        n_meas = latched_occupation(n_eq, lead.rate_hz(flat[:, 3]), dwell, rng)
        i = sensor.current_a(sensor.potential(flat, n_meas))
        ctx.count(model_evals=len(p1) * len(vbl))
    shape = (len(vbl), len(p1))
    return {
        "v_p1_v": p1,
        "v_bl_v": vbl,
        "current_a": i.reshape(shape).T,
        "n1_measured": n_meas[:, 0].reshape(shape).T,
        "n1_equilibrium": n_eq[:, 0].reshape(shape).T,
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
