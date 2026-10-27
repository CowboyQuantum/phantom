"""Step 7 protocols: the compiled device round trip, the virtual-gate lattice, and the sensor optimum."""

from __future__ import annotations

import numpy as np

from gallery.protocols.step2 import _layout
from gallery.runner import Context, Protocol
from phantom.compile import CompileParams, Device, compile_device
from phantom.optimum import sensor_optimum, virtual_gate_matrix


def _compiled(p: Protocol) -> Device:
    inst = p.instance
    lay = _layout(p)
    params = CompileParams(
        operating_point_v={k: float(v) for k, v in inst["operating_point_v"].items()},
        temperature_k=float(inst.get("temperature_k", 0.15)),
        severity=float(inst.get("severity", 1.0)),
        wiring_faults=bool(inst.get("wiring_faults", True)),
    )
    return compile_device(lay, params, int(inst.get("seed", 2)))


def _plunger_grid(dev: Device, v0: np.ndarray, half: float, n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    i1, i2 = dev.gate_names.index("P1"), dev.gate_names.index("P2")
    p1 = v0[i1] + np.linspace(-half, half, n)
    p2 = v0[i2] + np.linspace(-half, half, n)
    vv = np.broadcast_to(v0, (n, n, len(v0))).copy()
    vv[..., i1] = p1[:, None]
    vv[..., i2] = p2[None, :]
    return p1, p2, vv


def compiled_device_roundtrip(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Compile a device from the layout and a seed, serialise it to JSON text, reload it, and draw the
    coarse stability map from both: the two maps are bit-identical (their difference is all zeros).
    The compiled parameters that set the picture are recorded (not the file's byte count, whose last
    digits differ between platforms' float printing)."""
    with ctx.compile():
        dev = _compiled(p)
        text = dev.dumps()
        again = Device.loads(text)
        v0 = dev.operating_point()
        half, n = float(p.sweep[0]["half_span_v"]), int(p.sweep[0]["points"])
        p1, p2, vv = _plunger_grid(dev, v0, half, n)
    with ctx.run(2 * n * n):
        i_a = dev.runtime_sensed(noise=False).sweep_a(vv, None)
        i_b = again.runtime_sensed(noise=False).sweep_a(vv, None)
        ctx.count(model_evals=2 * n * n)
    return {
        "v_p1_v": p1,
        "v_p2_v": p2,
        "current_a": i_a,
        "difference_a": i_b - i_a,
        "tunnel_coupling_uev": np.array(dev.tunnel_coupling_ev * 1e6),
        "lead_rates_hz": np.array([b.rate_hz for b in dev.barriers]),
        "valley_splittings_uev": dev.valley_splittings_uev(),
        "g_factors": np.array(dev.g_factors),
        "operating_point_v": v0,
    }


def virtual_gates_map(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """The compiled device's coarse map in the raw plungers and in the virtual plungers of
    `optimum.virtual_gate_matrix`: the honeycomb's lines, tilted by the cross lever arms in raw
    voltages, come out axis-aligned in virtual coordinates because each virtual plunger shifts one
    dot's potential alone."""
    with ctx.compile():
        dev = _compiled(p)
        m = virtual_gate_matrix(dev)
        v0 = dev.operating_point()
        half, n = float(p.sweep[0]["half_span_v"]), int(p.sweep[0]["points"])
        p1, p2, vv = _plunger_grid(dev, v0, half, n)
        i1, i2 = dev.gate_names.index("P1"), dev.gate_names.index("P2")
        u = np.zeros((n, n, len(v0)))
        u[..., i1] = p1[:, None] - v0[i1]
        u[..., i2] = p2[None, :] - v0[i2]
        vv_virtual = v0 + u @ m.T
    with ctx.run(2 * n * n):
        sensed = dev.runtime_sensed(noise=False)
        i_raw = sensed.sweep_a(vv, None)
        i_virtual = sensed.sweep_a(vv_virtual, None)
        ctx.count(model_evals=2 * n * n)
    return {"v_p1_v": p1, "v_p2_v": p2, "current_raw_a": i_raw, "u_p1_v": p1 - v0[i1], "u_p2_v": p2 - v0[i2], "current_virtual_a": i_virtual, "virtual_gate_matrix": m}


def optimum_sensor_flank(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """The compiled device's sensor current and its slope against the sensor gate over one peak period,
    and the point `optimum.sensor_optimum` chooses (the steepest flank), marked by a spike in the
    `optimum_marker` trace."""
    with ctx.compile():
        dev = _compiled(p)
        opt = sensor_optimum(dev)
        j = dev.gate_names.index("SG")
        v0 = dev.operating_point()
        sg = v0[j] + np.linspace(float(p.sweep[0]["from"]), float(p.sweep[0]["to"]), int(p.sweep[0]["points"]))
        vv = np.broadcast_to(v0, (len(sg), len(v0))).copy()
        vv[:, j] = sg
    with ctx.run(len(sg)):
        i = dev.runtime_sensed(noise=False).sweep_a(vv, None)
        ctx.count(model_evals=len(sg))
    slope = np.abs(np.gradient(i, sg))
    marker = np.where(np.abs(sg - opt.v_sensor_gate_v) <= 0.5 * (sg[1] - sg[0]), slope.max(), 0.0)
    return {"v_sg_v": sg, "current_a": i, "slope_a_per_v": slope, "optimum_marker": marker, "optimum_v_sg_v": np.array(opt.v_sensor_gate_v), "optimum_slope_a_per_v": np.array(opt.slope_a_per_v)}
