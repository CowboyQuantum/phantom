"""`optimum`: the best a compiled device can do, layer by layer (ADR 0023).

Everything here is cheap and reads the `Device`'s runtime objects; nothing calls a solver. The
answers are what a perfect tuner would find:

- **sensor** — the sensor-gate voltage where the sensor current is steepest (the flank of its Coulomb
  peak), found by a coarse scan over one peak period, a fine scan around its best point and a
  parabolic interpolation of |dI/dV|;
- **lattice** — the virtual-gate matrix: the plungers are replaced by linear combinations that each
  move one dot's potential alone, at the raw plunger's scale (inv(L_pp)·diag(L_pp)), that compensate
  the other gates' cross-talk, and that hold the sensor potential fixed, so the honeycomb's lines come
  out axis-aligned on a flat sensor background;
- **coupling** — the barrier-gate voltage that gives a target interdot tunnel coupling, from the
  compiled exponential t(V);
- **readout** — the Pauli-blockade read point with the largest singlet–triplet contrast for the
  receiving dot's valley draw, and the time at which the Elzerman empty-probability peaks for the
  compiled lead rates;
- **coherence** — the exchange sweet spot (the symmetric point, where dJ/dε = 0), the single-spin T₂*
  of the draw, and the SWAP fidelity at the sweet spot with the device's noise.

`optimum` is attained by brute force on a small device (P0). Consumers that hide truth do not
forward it.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np

from phantom.constants import K_B_EV_PER_K
from phantom.device import Device
from phantom.spin.dynamics import SpinParams, quasi_static_sigma_ev, t2_star_s
from phantom.spin.elzerman import empty_probability
from phantom.spin.exchange import detuning_noise_sigma_ev, symmetric_point_v
from phantom.spin.gates import swap_gate_fidelity
from phantom.spin.psb import best_read_point

MODEL_VERSION = 1


@dataclass(frozen=True)
class SensorOptimum:
    v_sensor_gate_v: float
    slope_a_per_v: float
    """|dI/dV| at the chosen point."""


@dataclass(frozen=True)
class LatticeOptimum:
    virtual_gate_matrix: tuple[tuple[float, ...], ...]
    """V = M u + v₀: programmed gate voltages from virtual coordinates u (one per gate)."""
    plungers: tuple[str, str]


@dataclass(frozen=True)
class CouplingOptimum:
    t_target_ev: float
    v_barrier_v: float


@dataclass(frozen=True)
class ReadoutOptimum:
    psb_detuning_ev: float
    psb_contrast: float
    elzerman_peak_time_s: float
    elzerman_peak_probability: float


@dataclass(frozen=True)
class CoherenceOptimum:
    sweet_spot_v: tuple[float, ...]
    """Gate voltages of the symmetric point of the (1,1) cell."""
    t2_star_single_spin_s: float
    swap_fidelity_at_sweet_spot: float
    swap_gate_time_s: float


@dataclass(frozen=True)
class Optimum:
    sensor: SensorOptimum
    lattice: LatticeOptimum
    coupling: CouplingOptimum
    readout: ReadoutOptimum
    coherence: CoherenceOptimum

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# --- sensor ---------------------------------------------------------------------------------------


def _sensor_slope(device: Device, v_sg: np.ndarray, dv: float = 1e-5) -> np.ndarray:
    sensed = device.runtime_sensed(noise=False)
    j = device.gate_names.index(device.params["sensor_gate"])
    v0 = device.operating_point()
    vv = np.broadcast_to(v0, v_sg.shape + (len(v0),)).copy()
    vv[..., j] = v_sg
    vp, vm = vv.copy(), vv.copy()
    vp[..., j] += dv
    vm[..., j] -= dv
    return np.abs(sensed.sweep_a(vp, None) - sensed.sweep_a(vm, None)) / (2 * dv)


def sensor_optimum(device: Device, coarse_points: int = 201, fine_points: int = 201) -> SensorOptimum:
    """The steepest point of the sensor current over one peak period around the compiled set-point: a
    coarse scan, a fine scan of one coarse step around its best point, and a parabolic interpolation
    through the fine scan's best three points. (A golden-section search converged to a point that
    differed between platforms by 1e-7 relative, because it compares slopes equal to the last bit on
    a flat maximum; the grid-and-parabola answer is deterministic to ADR 0014's tolerance.)"""
    sensor = device.runtime_sensor()
    j = device.gate_names.index(device.params["sensor_gate"])
    v0 = device.operating_point()[j]
    alpha = sensor.alpha_v_per_v[j]
    period_v = sensor.period_v / max(alpha, 1e-9)
    grid = v0 + np.linspace(-period_v / 2, period_v / 2, coarse_points)
    slope = _sensor_slope(device, grid)
    k = int(np.argmax(slope))
    step = grid[1] - grid[0]
    fine = grid[k] + np.linspace(-step, step, fine_points)
    fs = _sensor_slope(device, fine)
    m = int(np.clip(np.argmax(fs), 1, fine_points - 2))
    y0, y1, y2 = fs[m - 1], fs[m], fs[m + 1]
    denom = y0 - 2 * y1 + y2
    offset = 0.5 * (y0 - y2) / denom if denom < 0 else 0.0
    x = fine[m] + offset * (fine[1] - fine[0])
    return SensorOptimum(float(x), float(_sensor_slope(device, np.array([x]))[0]))


# --- lattice --------------------------------------------------------------------------------------


def virtual_gate_matrix(device: Device, plungers: tuple[str, str] = ("P1", "P2"), compensate_sensor: bool = True) -> np.ndarray:
    """M (gates × gates) with V = v₀ + M u: virtual plunger i shifts dot i alone, with the same lever
    arm on it as the raw plunger has (a virtual volt is worth a raw volt on its own dot and nothing on
    the other); every other gate's column compensates its cross-talk on both dots through the
    plungers; with `compensate_sensor` the sensor gate moves along so that the sensor potential stays
    fixed, and the plungers absorb the sensor gate's own lever arms on the dots. The constraints (dot
    energies, sensor potential) are solved jointly for the plungers and the sensor gate; the diagonal
    of the remaining gates is one."""
    lever = device.runtime_ci().lever_arms_ev_per_v()  # [dot, gate]
    names = device.gate_names
    ip = [names.index(g) for g in plungers]
    sensor_gate = device.params.get("sensor_gate", "SG")
    alpha = np.array(device.sensor_alpha_v_per_v)
    use_sensor = compensate_sensor and sensor_gate in names and alpha[names.index(sensor_gate)] > 0
    moving = ip + ([names.index(sensor_gate)] if use_sensor else [])
    # rows: the two dot energies, plus the sensor potential when it is held; columns: the moving gates
    a = np.concatenate([lever[:, moving], alpha[moving][None, :]], axis=0) if use_sensor else lever[:, moving]
    a_inv = np.linalg.inv(a)
    n = len(names)
    m = np.eye(n)
    target = np.diag(np.diag(lever[:, ip]))
    rhs = np.concatenate([target, np.zeros((1, 2))], axis=0) if use_sensor else target
    m[np.ix_(moving, ip)] = a_inv @ rhs
    for j in range(n):
        if j in moving:
            continue
        rhs_j = np.concatenate([-lever[:, j], [-alpha[j]]]) if use_sensor else -lever[:, j]
        m[moving, j] = a_inv @ rhs_j
    return m


def lattice_optimum(device: Device, plungers: tuple[str, str] = ("P1", "P2")) -> LatticeOptimum:
    m = virtual_gate_matrix(device, plungers)
    return LatticeOptimum(tuple(tuple(float(x) for x in row) for row in m), plungers)


# --- coupling -------------------------------------------------------------------------------------


def coupling_optimum(device: Device, t_target_ev: float) -> CouplingOptimum:
    b = device.runtime_barrier()
    v = b.v_ref_v + b.v_scale_v * np.log(t_target_ev / b.t0_ev)
    return CouplingOptimum(float(t_target_ev), float(v))


# --- readout --------------------------------------------------------------------------------------


def readout_optimum(device: Device, t1_read_s: float = 1e-3) -> ReadoutOptimum:
    kt = K_B_EV_PER_K * device.temperature_k
    e_v = device.dots[1].valley_splitting_uev * 1e-6
    eps = np.linspace(-200e-6, e_v + 200e-6, 4001)
    e_best, contrast = best_read_point(eps, e_v, device.tunnel_coupling_ev, kt)
    rates = {b.between: b.rate_hz for b in device.barriers}
    g_out = rates[("dot1", "lead_left")]
    g_in = g_out
    lam = g_out + 1.0 / t1_read_s
    t_peak = np.log(lam / g_in) / (lam - g_in) if abs(lam - g_in) > 1e-12 * lam else 1.0 / lam
    p_peak = float(empty_probability(np.array([t_peak]), g_out, g_in, t1_read_s)[0])
    return ReadoutOptimum(float(e_best), float(contrast), float(t_peak), p_peak)


# --- coherence ------------------------------------------------------------------------------------


def coherence_optimum(device: Device, b_t: float = 0.5, spin: SpinParams | None = None) -> CoherenceOptimum:
    spin = spin or SpinParams(**{k: v for k, v in device.params["spin"].items()})
    hub = device.runtime_hubbard()
    fixed = {g: device.operating_point_v[i] for i, g in enumerate(device.gate_names) if g not in ("P1", "P2")}
    v_sweet = symmetric_point_v(hub, fixed)
    sigma_v = device.noise_sigma_gate_v()
    sigma_eps = detuning_noise_sigma_ev(hub, sigma_v)
    f, t_gate = swap_gate_fidelity(hub, v_sweet, 0.0, sigma_eps, spin)
    t2 = t2_star_s(quasi_static_sigma_ev(spin, b_t, sigma_v))
    return CoherenceOptimum(tuple(float(x) for x in v_sweet), float(t2), float(f), float(t_gate))


def optimum(device: Device, t_target_ev: float = 20e-6) -> Optimum:
    return Optimum(sensor_optimum(device), lattice_optimum(device), coupling_optimum(device, t_target_ev), readout_optimum(device), coherence_optimum(device))
