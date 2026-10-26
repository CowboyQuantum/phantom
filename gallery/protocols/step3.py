"""Step 3 protocols: the clean device beside disorder draws, and the statistics of unintentional dots."""

from __future__ import annotations

import numpy as np

from gallery.protocols.step2 import _layout, _operating_point, gate_outlines_nm
from gallery.runner import Context, Protocol
from phantom.charge import derive_double_dot
from phantom.disorder import ScenarioParams, apply_scenario, draw_scenario
from phantom.electrostatics import GateBasis
from phantom.electrostatics.dots import _region_mask, local_minima
from phantom.layout import MaterialClass
from phantom.readout import SensedDoubleDot, SensorSETParams, raster_times_s


def _draws(p: Protocol):
    lay = _layout(p)
    params = ScenarioParams.for_material(lay.material)
    seeds = [int(s) for s in p.instance.get("seeds", [1, 2, 3])]
    severity = float(p.instance.get("severity", 1.0))
    return lay, params, seeds, severity


def potential_map_disorder(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """The electron energy landscape of the clean layout and of three disorder draws at the same
    operating point: interface charges as local dips and bumps, the whole gate layer shifted."""
    with ctx.compile():
        lay, params, seeds, severity = _draws(p)
        v = _operating_point(lay, p)
        clean = GateBasis.from_layout(lay)
        realisations = [apply_scenario(lay, draw_scenario(lay, params, s, severity)) for s in seeds]
    n = clean.basis.shape[1] * clean.basis.shape[2]
    out = {"x_nm": clean.x_nm, "y_nm": clean.y_nm, "u_clean_mev": (clean.energy_ev(v) * 1e3).T, "gate_outlines_nm": gate_outlines_nm(lay)}
    with ctx.run(n * (1 + len(seeds))):
        for s, (perturbed, basis) in zip(seeds, realisations, strict=True):
            out[f"u_seed{s}_mev"] = (basis.energy_ev(v) * 1e3).T
            out[f"gate_outlines_seed{s}_nm"] = gate_outlines_nm(perturbed)
        ctx.count(model_evals=n * (1 + len(seeds)))
    return out


def _sensed(lay, v, derived, seed: int):
    sensor = SensorSETParams(alpha_v_per_v=derived.sensor_alpha_v_per_v)
    return SensedDoubleDot.default_noisy(seed=seed, charge=derived.ci, sensor=sensor)


def _flank(v, derived, sensor_gate_index: int, widths: float = 0.9):
    sensor = SensorSETParams(alpha_v_per_v=derived.sensor_alpha_v_per_v)
    n0 = derived.ci.mean_occupation(v)
    phi0 = float(sensor.potential(v, n0))
    k = np.rint((phi0 - sensor.phase_v) / sensor.period_v)
    target = sensor.phase_v + k * sensor.period_v + widths * sensor.width_v
    v = v.copy()
    v[sensor_gate_index] += (target - phi0) / sensor.alpha_v_per_v[sensor_gate_index]
    return v


def stability_coarse_2dot_disorder(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """The layout-derived coarse stability map for the clean layout and three disorder draws, each
    re-derived (dots, lever arms, capacitances, sensor coupling) from its own perturbed landscape and
    with its sensor set on its flank. Lines tilt and shift, periods change, and an extra minimum shows
    as an additional set of lines."""
    with ctx.compile():
        lay, params, seeds, severity = _draws(p)
        v0 = _operating_point(lay, p)
        j_sg = lay.gate_names.index(p.instance.get("sensor_gate", "SG"))
        temperature = float(p.instance.get("temperature_k", 0.15))
        cases = [("clean", lay, GateBasis.from_layout(lay))]
        for s in seeds:
            perturbed, basis = apply_scenario(lay, draw_scenario(lay, params, s, severity))
            cases.append((f"seed{s}", perturbed, basis))
        devices = []
        for name, layout_c, basis_c in cases:
            derived = derive_double_dot(layout_c, v0, temperature_k=temperature, basis=basis_c)
            v = _flank(v0, derived, j_sg)
            derived = derive_double_dot(layout_c, v, temperature_k=temperature, basis=basis_c)
            devices.append((name, v, _sensed(layout_c, v, derived, int(p.instance.get("seed", 0)))))
    ax1, ax2 = p.sweep
    i1, i2 = lay.gate_names.index(ax1["gate"]), lay.gate_names.index(ax2["gate"])
    n1, n2 = int(ax1["points"]), int(ax2["points"])
    out: dict[str, np.ndarray] = {}
    with ctx.run(n1 * n2 * len(devices)):
        for name, v, dev in devices:
            p1 = v[i1] + np.linspace(float(ax1["from"]), float(ax1["to"]), n1)
            p2 = v[i2] + np.linspace(float(ax2["from"]), float(ax2["to"]), n2)
            vv = np.broadcast_to(v, (n1, n2, len(v))).copy()
            vv[..., i1] = p1[:, None]
            vv[..., i2] = p2[None, :]
            times = None if p.instance.get("noise") == "none" else raster_times_s(vv.shape[:-1], float(p.readout.get("dt_s", 1e-3)))
            out[f"v_p1_{name}_v"] = p1
            out[f"v_p2_{name}_v"] = p2
            out[f"current_{name}_a"] = dev.sweep_a(vv, times)
        ctx.count(model_evals=n1 * n2 * len(devices))
    return out


def unintentional_dots(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Number of local minima of the electron energy inside the channel region, beyond the two intended
    dots, averaged over many draws, against the charge-density severity. The P1 row of Step 3 compares
    this curve with a published MOS versus SiGe comparison once its reference entry exists."""
    with ctx.compile():
        lay, params, _, _ = _draws(p)
        v = _operating_point(lay, p)
        mat = MaterialClass.by_name(lay.material)
        severities = np.array([float(s) for s in p.instance.get("severities", [0.0, 1.0, 3.0, 10.0, 30.0])])
        n_draws = int(p.instance.get("n_draws", 40))
        del mat
    counts = np.zeros((len(severities), n_draws))
    with ctx.run(len(severities) * n_draws):
        for i, sev in enumerate(severities):
            for k in range(n_draws):
                _, basis = apply_scenario(lay, draw_scenario(lay, params, 1000 + k, float(sev)))
                u = basis.energy_ev(v)
                minima = local_minima(u, _region_mask(basis, lay.channel_nm))
                counts[i, k] = len(minima)
        ctx.count(model_evals=len(severities) * n_draws)
    extra = np.maximum(counts - 2, 0)
    return {
        "severity": severities,
        "density_cm2": severities * params.fixed_charges.density_cm2,
        "mean_extra_minima": extra.mean(axis=1),
        "p90_extra_minima": np.percentile(extra, 90, axis=1),
        "fraction_with_extra": (extra > 0).mean(axis=1),
    }
