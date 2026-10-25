"""Per-model micro-benchmarks: cost of one evaluation of each runtime model (ADR 0009, 0011).

Each entry times a representative call on a fixed workload, best of three, and reports microseconds per
point or sample. Numbers are machine-dependent; `docs/gallery/measured.json` records them with the
commit, platform and date, and `docs/MODELS.md` renders them next to the protocol budgets.
"""

from __future__ import annotations

import time
from collections.abc import Callable

import numpy as np

from phantom.charge import Barrier, ConstantInteractionParams, HubbardParams
from phantom.readout import (
    LeadRateLandscape,
    LockinParams,
    OneOverFNoise,
    SensedDoubleDot,
    TelegraphNoise,
    latched_occupation,
    lockin_sweep_a_per_v,
    lockin_sweep_dynamic_a_per_v,
    raster_times_s,
)


def _grid(n: int, n_g: int = 3) -> np.ndarray:
    side = int(np.sqrt(n))
    v = np.zeros((side, side, n_g))
    v[..., 0] = np.linspace(-0.05, 0.15, side)[:, None]
    v[..., 1] = np.linspace(-0.05, 0.15, side)[None, :]
    return v


def workloads() -> list[tuple[str, str, int, Callable[[], object]]]:
    """(model, what is timed, points, callable)."""
    ci = ConstantInteractionParams()
    hb = HubbardParams.from_ci(ci, Barrier("B", t0_ev=50e-6))
    d = SensedDoubleDot.default_noisy(seed=0)
    v3 = _grid(10_000)
    v4 = _grid(10_000, 4)
    times = raster_times_s(v3.shape[:-1], 1e-3)
    lk = LockinParams(gate="P1", cycles=4)
    v_lk = _grid(2_500)
    rates = LeadRateLandscape(ci.gate_names)
    phi = np.linspace(-0.6, 0.6, 100_000)
    n_eq = np.stack([np.repeat(np.arange(100.0), 1000), np.zeros(100_000)], axis=1)
    return [
        ("charge.ci.ConstantInteractionParams.mean_occupation", "thermal occupation per gate point", 10_000, lambda: ci.mean_occupation(v3)),
        ("charge.hubbard.HubbardParams.mean_occupation", "16x16 diagonalisation per gate point", 10_000, lambda: hb.mean_occupation(v4)),
        ("readout.sensor.SensorSETParams.current_a", "sensor current per point", 100_000, lambda: d.sensor.current_a(phi)),
        ("readout.noise.TelegraphNoise.trace", "one fluctuator per sample", 100_000, lambda: TelegraphNoise(1.0, 10.0, 10.0, seed=0).trace(np.arange(100_000) * 1e-3)),
        ("readout.noise.OneOverFNoise.trace", "24 fluctuators per sample", 10_000, lambda: OneOverFNoise(1.0, seed=0).trace(np.arange(10_000) * 1e-3)),
        ("readout.measure.SensedDoubleDot.sweep_a (noisy)", "occupation + sensor + three noise traces per point", 10_000, lambda: d.sweep_a(v3, times)),
        ("readout.lockin.lockin_sweep_a_per_v", "8 evaluations per point, noisy", 2_500, lambda: lockin_sweep_a_per_v(d, v_lk, lk, t0_s=0.0)),
        ("readout.lockin.lockin_sweep_dynamic_a_per_v", "8 evaluations + rate dynamics per point, noisy", 2_500, lambda: lockin_sweep_dynamic_a_per_v(d, v_lk, lk, rates, np.random.default_rng(0), t0_s=0.0)),
        ("readout.loading.latched_occupation", "per scan step", 100_000, lambda: latched_occupation(n_eq, np.full(100_000, 1e3), 1e-3, np.random.default_rng(0))),
    ]


def run(repeats: int = 3) -> list[dict[str, object]]:
    rows = []
    for model, what, points, fn in workloads():
        best = min(_time(fn) for _ in range(repeats))
        rows.append({"model": model, "what": what, "points": points, "us_per_point": round(1e6 * best / points, 3)})
    return rows


def _time(fn: Callable[[], object]) -> float:
    t0 = time.perf_counter()
    fn()
    return time.perf_counter() - t0
