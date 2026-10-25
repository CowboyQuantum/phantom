"""Step 2 protocols: the potential landscape of a layout, and measurements derived from it alone."""

from __future__ import annotations

import numpy as np

from gallery.runner import Context, Protocol
from phantom.charge import derive_double_dot
from phantom.electrostatics import GateBasis, line_cut
from phantom.layout import Layout, two_dot_layout
from phantom.readout import SensedDoubleDot, SensorSETParams, raster_times_s


def _layout(p: Protocol) -> Layout:
    inst = p.instance
    if inst.get("layout_file"):
        return Layout.from_json(inst["layout_file"])
    return two_dot_layout(inst.get("material", "SiGe"), float(inst.get("pitch_nm", 80.0)))


def _operating_point(layout: Layout, p: Protocol) -> np.ndarray:
    op = p.instance["operating_point_v"]
    return np.array([float(op[name]) for name in layout.gate_names])


def gate_outlines_nm(layout: Layout) -> np.ndarray:
    """Closed outlines of every gate as one NaN-separated polyline (m, 2) for overlays."""
    rows = []
    for g in layout.gates:
        v = g.vertices
        rows.append(np.concatenate([v, v[:1], [[np.nan, np.nan]]], axis=0))
    return np.concatenate(rows, axis=0)


def potential_map_layout(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    with ctx.compile():
        lay = _layout(p)
        basis = GateBasis.from_layout(lay, float(p.instance.get("grid_pitch_nm", 2.0)))
        v = _operating_point(lay, p)
    with ctx.run(basis.basis.shape[1] * basis.basis.shape[2]):
        u = basis.energy_ev(v) * 1e3
        ctx.count(model_evals=basis.basis.shape[1] * basis.basis.shape[2])
    return {"x_nm": basis.x_nm, "y_nm": basis.y_nm, "u_mev": u.T, "gate_outlines_nm": gate_outlines_nm(lay)}


def potential_cut_along_channel(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    with ctx.compile():
        lay = _layout(p)
        basis = GateBasis.from_layout(lay, float(p.instance.get("grid_pitch_nm", 2.0)))
        v0 = _operating_point(lay, p)
        gate = p.instance.get("vary_gate", "B_M")
        values = [float(x) for x in p.instance.get("vary_values_v", [0.0, 0.2, 0.4])]
        x0, x1 = (float(x) for x in p.instance.get("x_range_nm", [-200.0, 200.0]))
    samples = int(p.sweep[0]["points"])
    out: dict[str, np.ndarray] = {}
    with ctx.run(samples * len(values)):
        j = lay.gate_names.index(gate)
        for val in values:
            v = v0.copy()
            v[j] = val
            s, u = line_cut(basis, v, x0, 0.0, x1, 0.0, samples)
            out[f"u_{gate}_{val:g}V_mev".replace(".", "p").replace("-", "m")] = u * 1e3
        out["x_nm"] = x0 + s
        ctx.count(model_evals=samples * len(values))
    return out


def _sensed_from_layout(p: Protocol):
    """Derive the double dot from the layout, then set the sensor gate so the sensor sits on the steep
    flank of its nearest Coulomb peak (the first thing a tune-up does with a sensor), and re-derive at
    that operating point. The sensor voltage actually used is returned with the arrays."""
    lay = _layout(p)
    v = _operating_point(lay, p)
    temperature = float(p.instance.get("temperature_k", 0.15))
    derived = derive_double_dot(lay, v, temperature_k=temperature)
    sensor = SensorSETParams(alpha_v_per_v=derived.sensor_alpha_v_per_v)
    flank = float(p.instance.get("sensor_flank_widths", 0.9))
    j_sg = lay.gate_names.index(p.instance.get("sensor_gate", "SG"))
    if sensor.alpha_v_per_v[j_sg] > 0:
        n0 = derived.ci.mean_occupation(v)
        phi0 = float(sensor.potential(v, n0))
        k = np.rint((phi0 - sensor.phase_v) / sensor.period_v)
        target = sensor.phase_v + k * sensor.period_v + flank * sensor.width_v
        v = v.copy()
        v[j_sg] += (target - phi0) / sensor.alpha_v_per_v[j_sg]
        derived = derive_double_dot(lay, v, temperature_k=temperature)
        sensor = SensorSETParams(alpha_v_per_v=derived.sensor_alpha_v_per_v)
    seed = int(p.instance.get("seed", 0))
    if p.instance.get("noise") == "none":
        dev = SensedDoubleDot(charge=derived.ci, sensor=sensor)
    else:
        dev = SensedDoubleDot.default_noisy(seed=seed, charge=derived.ci, sensor=sensor)
    return lay, v, derived, dev


def stability_coarse_2dot_layout(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """The coarse stability map with nothing typed but the layout, the operating point and the sensor's
    per-electron coupling: lever arms, capacitances and the sensor's gate coupling come from the basis."""
    with ctx.compile():
        lay, v0, derived, dev = _sensed_from_layout(p)
    ax1, ax2 = p.sweep
    i1, i2 = lay.gate_names.index(ax1["gate"]), lay.gate_names.index(ax2["gate"])
    p1 = v0[i1] + np.linspace(float(ax1["from"]), float(ax1["to"]), int(ax1["points"]))
    p2 = v0[i2] + np.linspace(float(ax2["from"]), float(ax2["to"]), int(ax2["points"]))
    v = np.broadcast_to(v0, (len(p1), len(p2), len(v0))).copy()
    v[..., i1] = p1[:, None]
    v[..., i2] = p2[None, :]
    with ctx.run(len(p1) * len(p2)):
        times = None if p.instance.get("noise") == "none" else raster_times_s(v.shape[:-1], float(p.readout.get("dt_s", 1e-3)))
        i = dev.sweep_a(v, times)
        ctx.count(model_evals=len(p1) * len(p2))
    return {
        "v_p1_v": p1,
        "v_p2_v": p2,
        "current_a": i,
        "lever_arms": derived.ci.lever_arms_ev_per_v(),
        "c_dot_af": derived.ci.c_dot_af,
        "c_gate_af": derived.ci.c_gate_af,
        "operating_point_v": v0,
    }


def pinch_off_crosstalk(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Current through the channel under one barrier against that barrier and a neighbouring gate.
    The channel is open where the electron energy under the barrier lies below a threshold; because
    the energy is linear in every gate voltage, the pinch-off voltage shifts linearly with the
    neighbour, with slope minus the ratio of the two lever arms at the barrier. The current shape is a
    soft step of width `width_ev` (phenomenological, calibration none)."""
    with ctx.compile():
        lay = _layout(p)
        basis = GateBasis.from_layout(lay, float(p.instance.get("grid_pitch_nm", 2.0)))
        v0 = _operating_point(lay, p)
        barrier, neighbour = p.instance["barrier"], p.instance["neighbour"]
        x_b, _ = lay.gate(barrier).centroid_nm
        lever = basis.at(x_b, 0.0)  # at the channel under the barrier's centre
        u_th = float(p.instance.get("threshold_ev", -0.05))
        width = float(p.instance.get("width_ev", 0.004))
        i_sat = float(p.instance.get("i_sat_a", 2e-9))
    ib, inb = lay.gate_names.index(barrier), lay.gate_names.index(neighbour)
    vb = v0[ib] + np.linspace(float(p.sweep[0]["from"]), float(p.sweep[0]["to"]), int(p.sweep[0]["points"]))
    vn = v0[inb] + np.linspace(float(p.sweep[1]["from"]), float(p.sweep[1]["to"]), int(p.sweep[1]["points"]))
    with ctx.run(len(vb) * len(vn)):
        v = np.broadcast_to(v0, (len(vb), len(vn), len(v0))).copy()
        v[..., ib] = vb[:, None]
        v[..., inb] = vn[None, :]
        u = -(v @ lever)  # eV under the barrier
        current = i_sat / (1.0 + np.exp((u - u_th) / width))
        ctx.count(model_evals=len(vb) * len(vn))
    others = sum(lever[j] * v0[j] for j in range(len(v0)) if j not in (ib, inb))
    v_po = (-u_th - others - lever[inb] * vn) / lever[ib]
    return {
        "v_barrier_v": vb,
        "v_neighbour_v": vn,
        "current_a": current,
        "v_pinch_off_v": v_po,
        "slope_v_per_v": np.array(-lever[inb] / lever[ib]),
    }
