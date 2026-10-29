"""Step 7 protocols, third set: the sensor dot's barriers and bath on the compiled device (ADR 0026, #47)."""

from __future__ import annotations

import numpy as np

from gallery.protocols.step2 import _axis
from gallery.protocols.step7 import _compiled
from gallery.runner import Context, Protocol
from phantom.constants import E_CHARGE_C, HBAR_EV_S
from phantom.electrostatics.dots import DerivedDot, barrier_action

CONDUCTANCE_QUANTUM_S = E_CHARGE_C / (2 * np.pi * HBAR_EV_S)
"""e²/h in siemens (e over h in eV·s)."""
CURRENT_FLOOR_A = 1e-15
"""A floor for logarithmic current panels, below any measurement (1 fA)."""


def _gate_grid(v0: np.ndarray, j: int, values: np.ndarray) -> np.ndarray:
    v = np.broadcast_to(v0, values.shape + v0.shape).copy()
    v[..., j] = values
    return v


def sensor_bath_turn_on(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """Each of the sensor's bath gates swept from depletion to accumulation with every other gate at the
    compiled operating point: the bath's occupation (the soft step at the accumulation threshold) and
    the sensor current. Below the turn-on the bath holds no electrons and the sensor carries none; above
    it the bath gate keeps lowering its side of the barrier and moving the sensor along its peaks."""
    with ctx.compile():
        dev = _compiled(p)
        leads = dev.sensor_leads
        v0 = dev.operating_point()
        x = _axis(p.sweep[0])
        sensed = dev.runtime_sensed(noise=False)
    out: dict[str, np.ndarray] = {"v_bath_v": x}
    turn_on = []
    with ctx.run(2 * len(x)):
        for k, (side, gate) in enumerate(zip(("left", "right"), leads.reservoir_gates, strict=True)):
            j = dev.gate_names.index(gate)
            v = _gate_grid(v0, j, x)
            out[f"current_{side}_bath_swept_a"] = sensed.current_a(v)
            out[f"log10_current_{side}_bath_swept_a"] = np.log10(np.maximum(out[f"current_{side}_bath_swept_a"], CURRENT_FLOOR_A))
            out[f"occupation_{side}_bath_swept"] = leads.bath_occupation(v)[:, k]
            lever = np.asarray(leads.bath_lever_arms[k])
            others = float(lever @ v0 - lever[j] * v0[j])
            turn_on.append((leads.bath_offset_ev[k] - leads.bath_threshold_ev - others) / lever[j])
        ctx.count(model_evals=2 * len(x))
    out["v_turn_on_v"] = np.array(turn_on)
    return out


def sensor_barrier_barrier_map(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """The sensor current against its two barrier gates at the compiled sensor-plunger voltage: no
    current where either barrier is closed (the pinch-off corner), Coulomb peaks running diagonally
    where both are open (each barrier also moves the dot's level), and the peaks broadening as the
    barriers open. The effective rate Γ_SΓ_D/(Γ_S + Γ_D) is mapped beside it."""
    with ctx.compile():
        dev = _compiled(p)
        v0 = dev.operating_point()
        ax1, ax2 = p.sweep
        i1, i2 = dev.gate_names.index(ax1["gate"]), dev.gate_names.index(ax2["gate"])
        b1 = v0[i1] + _axis(ax1)
        b2 = v0[i2] + _axis(ax2)
        v = np.broadcast_to(v0, (len(b1), len(b2), len(v0))).copy()
        v[..., i1] = b1[:, None]
        v[..., i2] = b2[None, :]
        sensed = dev.runtime_sensed(noise=False)
    with ctx.run(len(b1) * len(b2)):
        current = sensed.current_a(v)
        rates = sensed.sensor.rates_hz(v)
        gamma_eff = rates[..., 0] * rates[..., 1] / np.maximum(rates.sum(axis=-1), 1e-300)
        ctx.count(model_evals=len(b1) * len(b2))
    return {
        "v_sb_l_v": b1,
        "v_sb_r_v": b2,
        "current_a": current,
        "log10_gamma_eff_hz": np.log10(np.maximum(gamma_eff, 1e-300)),
        "operating_point_v": v0,
    }


def sensor_barrier_rates(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """The sensor's two tunnel rates against one sensor barrier gate: recomputed at the compile tier
    from the WKB action along each barrier of the landscape at the leads' fixed tunnelling level, beside
    the compiled exponential with its attempt-rate ceiling that the runtime uses. They agree at the
    anchor and drift apart as the barrier's shape changes away from it."""
    with ctx.compile():
        dev = _compiled(p)
        leads = dev.sensor_leads
        lay, basis = dev.basis()
        s = dev.sensor_dot
        dot = DerivedDot(s.name, s.x_nm, s.y_nm, s.energy_ev, s.curvature_ev_per_nm2, s.radius_nm, np.array(s.lever_arms))
        points = [lay.gate(n).centroid_nm for n in leads.reservoir_gates]
        v_ref = np.array(leads.v_ref_v)
        j = dev.gate_names.index(p.sweep[0]["gate"])
        vb = v_ref[j] + _axis(p.sweep[0])
        v = _gate_grid(v_ref, j, vb)
        compile_tier = np.empty((len(vb), 2))
        for i, vi in enumerate(v):
            for k, (point, level) in enumerate(zip(points, leads.tunnelling_level_ev, strict=True)):
                compile_tier[i, k] = leads.attempt_rate_hz * np.exp(-2.0 * barrier_action(basis, vi, dot, point, 0.19, level_ev=level)[2])
    with ctx.run(len(vb)):
        runtime = leads.barrier_rates_hz(v)
        ctx.count(model_evals=len(vb))
    return {
        "v_barrier_v": vb,
        "gamma_source_compile_hz": compile_tier[:, 0],
        "gamma_drain_compile_hz": compile_tier[:, 1],
        "gamma_source_runtime_hz": runtime[:, 0],
        "gamma_drain_runtime_hz": runtime[:, 1],
        "attempt_rate_hz": np.full(len(vb), leads.attempt_rate_hz),
    }


def sensor_peak_broadening(ctx: Context, p: Protocol) -> dict[str, np.ndarray]:
    """One sensor Coulomb peak against the sensor plunger, centred on the peak, with both sensor barriers
    moved together by each offset: closed barriers give a low, thermally broadened peak; opening them
    raises it until the lifetime width ħ(Γ_S + Γ_D) exceeds kT, after which it broadens and its height
    saturates at the conductance quantum (for symmetric rates). Per offset, the linear-response peak
    conductance in units of e²/h and the ratio of lifetime to thermal width."""
    inst = p.instance
    with ctx.compile():
        dev = _compiled(p)
        sensor = dev.runtime_sensor()
        v0 = dev.operating_point()
        jg = dev.gate_names.index("SG")
        jb = [dev.gate_names.index(g) for g in dev.sensor_leads.barrier_gates]
        offsets = np.array([float(x) for x in inst.get("barrier_offsets_v", [-0.02, -0.01, 0.0, 0.01, 0.02])])
        rel = _axis(p.sweep[0])
        n_dots = dev.runtime_ci().mean_occupation(v0)
        a = float(sensor.alpha_v_per_v[jg])
        settings = []
        for off in offsets:
            v = v0.copy()
            v[jb] += off
            phi = float(sensor.potential(v, n_dots))
            k = np.rint((phi - sensor.phase_v) / sensor.period_v)
            v[jg] += (sensor.phase_v + k * sensor.period_v - phi) / a  # onto the nearest peak
            settings.append(v)
    out: dict[str, np.ndarray] = {"sg_from_peak_v": rel, "barrier_offset_v": offsets}
    conductance, ratio = [], []
    with ctx.run(len(offsets) * len(rel)):
        for i, v in enumerate(settings):
            vv = _gate_grid(v, jg, v[jg] + rel)
            out[f"current_offset{i}_a"] = sensor.current_a(sensor.potential(vv, n_dots), None, vv)
            conductance.append(sensor.linear_response_peak_a(1e-7, v) / 1e-7 / CONDUCTANCE_QUANTUM_S)
            ratio.append(float(sensor.lifetime_width_v(v)) / sensor.width_v)
        ctx.count(model_evals=len(offsets) * len(rel))
    out["peak_conductance_quanta"] = np.array(conductance)
    out["lifetime_over_thermal_width"] = np.array(ratio)
    return out
