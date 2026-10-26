"""Step 4 protocols on a per-device noise draw: switcher dwell-time statistics and slow drift."""

from __future__ import annotations

import numpy as np

from gallery.runner import Context, Protocol
from phantom.disorder import NoiseParams, draw_noise
from phantom.readout import SensedDoubleDot, TelegraphNoise


def _drawn_device(p: Protocol, sample_dt_s: float):
    inst = p.instance
    seed = int(inst.get("seed", 0))
    draw = draw_noise(NoiseParams.for_material(inst.get("material", "SiGe")), np.random.default_rng(seed), seed, float(inst.get("severity", 1.0)))
    return draw, SensedDoubleDot.from_noise_draw(draw, sample_dt_s)


def _dwell_times(trace: np.ndarray, dt_s: float) -> tuple[np.ndarray, np.ndarray]:
    """Dwell times in the upper and lower level of a two-level trace."""
    hi = trace > trace.mean()
    edges = np.flatnonzero(np.diff(hi.astype(int)) != 0) + 1
    bounds = np.concatenate([[0], edges, [len(trace)]])
    lengths = np.diff(bounds) * dt_s
    levels = hi[bounds[:-1]]
    return lengths[levels], lengths[~levels]


def switcher_histogram(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """A long sensor trace on a drawn device and the dwell-time statistics of its dominant two-level
    fluctuator: the dwell times in each level are exponentially distributed with the drawn up and down
    rates, which a tune-up reads off exactly this way."""
    dt = float(p.readout["dt_s"])
    with ctx.compile():
        draw, dev = _drawn_device(p, dt)
        point = np.array([0.0, 0.0, dev.sensor.phase_v + 0.9 * dev.sensor.width_v])
        duration = float(p.readout["duration_s"])
        if not draw.switcher_amplitudes_v:
            raise ValueError("this seed drew no switcher; pick another seed")
        k = int(np.argmax(draw.switcher_amplitudes_v))
    t = np.arange(0.0, duration, dt)
    with ctx.run(len(t)):
        i = dev.trace_a(point, t)
        ctx.count(model_evals=len(t))
    sw = TelegraphNoise(draw.switcher_amplitudes_v[k], draw.switcher_rates_up_hz[k], draw.switcher_rates_down_hz[k], seed=draw.seed + 10 + k).trace(t)
    up_dwell, down_dwell = _dwell_times(sw, dt)
    edges = np.geomspace(dt, max(up_dwell.max(), down_dwell.max(), 10 * dt), 25)
    centres = np.sqrt(edges[1:] * edges[:-1])
    h_up, _ = np.histogram(up_dwell, edges)
    h_down, _ = np.histogram(down_dwell, edges)
    zoom = t < float(p.render.get("zoom_s", 20.0))
    return {
        "t_s": t[zoom],
        "current_a": i[zoom],
        "dwell_s": centres,
        "counts_up": h_up.astype(float) + 0.5,
        "counts_down": h_down.astype(float) + 0.5,
        "rate_down_hz": np.array(draw.switcher_rates_down_hz[k]),
        "rate_up_hz": np.array(draw.switcher_rates_up_hz[k]),
        "mean_dwell_up_s": np.array(up_dwell.mean()),
        "mean_dwell_down_s": np.array(down_dwell.mean()),
        "n_switchers": np.array(len(draw.switcher_amplitudes_v)),
        "one_over_f_amplitude_v": np.array(draw.one_over_f_amplitude_v),
    }


def drift_repeated_cuts(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """One P1 cut across a charge transition repeated every minute for two hours on a drawn device,
    with the noise clock running between cuts: the transition drifts with the 1/f ensemble's slowest
    fluctuators and jumps when a switcher flips."""
    inst = p.instance
    dt = float(p.readout["dt_s"])
    with ctx.compile():
        _, dev = _drawn_device(p, dt)
        n_pts = int(p.sweep[0]["points"])
        repeats = int(inst.get("repeats", 120))
        interval = float(inst.get("interval_s", 60.0))
        period = dev.charge.transition_period_v(0, "P1")
        p1 = 0.5 * period + np.linspace(-0.3 * period, 0.3 * period, n_pts)
        sg = dev.sensor.phase_v + 0.9 * dev.sensor.width_v
    v = np.zeros((repeats, n_pts, 3))
    v[..., 0] = p1[None, :]
    v[..., 2] = sg
    times = np.arange(repeats)[:, None] * interval + np.arange(n_pts)[None, :] * dt
    with ctx.run(repeats * n_pts):
        i = dev.sweep_a(v, times)
        ctx.count(model_evals=repeats * n_pts)
    # transition position per cut: where the smoothed cut crosses the midpoint between its two plateaus
    kernel = np.ones(9) / 9.0
    position = np.empty(repeats)
    for r in range(repeats):
        smooth = np.convolve(i[r], kernel, mode="same")
        lo, hi = np.median(smooth[: n_pts // 3]), np.median(smooth[-n_pts // 3 :])
        mid = 0.5 * (lo + hi)
        crossing = np.flatnonzero(np.diff(np.sign(smooth - mid)) != 0)
        position[r] = p1[crossing[0]] if len(crossing) else np.nan
    return {
        "v_p1_v": p1,
        "elapsed_min": np.arange(repeats) * interval / 60.0,
        "current_a": i.T,
        "transition_position_v": position,
        "transition_shift_uv": (position - position[0]) * 1e6,
    }
