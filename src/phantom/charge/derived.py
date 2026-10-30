"""Charge-model parameters derived from a layout instead of typed (Step 2, ADR 0013).

From a `Layout` and an operating gate vector: the gate basis gives lever arms at each dot minimum, the
curvature gives a dot radius and so a self-capacitance (`Dot`), the dot separation gives an estimate of
the mutual capacitance, and the capacitance picture of the CI model is completed by

    C_gate = C_dot · α        (α the dimensionless lever arms, [dot, gate]),

so that `ConstantInteractionParams.lever_arms_ev_per_v()` returns exactly the electrostatic lever arms.
The sensor's own lever arms per gate come from the basis at the sensor dot. The sensor is a third dot
of the same kind (ADR 0026): its total capacitance by the same rule, its mutual capacitance to each
qubit dot by the sphere formula, and from the inverse of the three-dot capacitance matrix

    E_c,s = e² (C⁻¹)_ss,      β_i = −(e / 1 V) (C⁻¹)_si   [V of the sensor potential per electron on dot i],

its charging energy and the shift of its potential per electron on each dot (an electron raises the
sensor's level, so β < 0). The qubit dots' 2×2 block is the double dot's own matrix, unchanged: the
sensor's charge does not act back on the dots (it changes during every sensor tunnelling event, faster
than any dot measurement resolves).

Approximations, stated: the dot is a point for lever arms (no wavefunction averaging); the mutual
capacitance of two dots at separation s is that of two spheres, 4πε r₁ r₂ / s (valid for s ≫ r); the
dot's own capacitance is the disk value from its harmonic radius, and its total capacitance is that
divided by (1 − Σ_j α_j) so that the gate capacitances C_dot·α are part of the total. Leads are not
counted. Step 7's solver replaces all of this.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phantom.charge.ci import ConstantInteractionParams
from phantom.charge.dot import EPSILON_0_F_PER_M, Dot
from phantom.constants import E_UNIT_V_PER_E_AF
from phantom.electrostatics.closed_form import GateBasis
from phantom.electrostatics.dots import DerivedBarrier, DerivedDot, barrier_between, find_dots
from phantom.layout.geometry import Layout, MaterialClass

MODEL_VERSION = 2
"""2: the sensor's charging energy and per-electron coupling from the three-dot capacitance matrix
(ADR 0026); every earlier output unchanged."""


@dataclass(frozen=True)
class DerivedDoubleDot:
    """What the layout gives: the basis, the two dots, the sensor dot, the barrier between the dots,
    the CI parameters built from them, and the sensor's charging energy and coupling to the dots."""

    basis: GateBasis
    dots: tuple[DerivedDot, DerivedDot]
    sensor_dot: DerivedDot | None
    interdot_barrier: DerivedBarrier
    ci: ConstantInteractionParams
    operating_point_v: np.ndarray
    sensor_charging_energy_ev: float | None = None
    """E_c,s = e²(C⁻¹)_ss of the sensor dot, eV; None without a sensor dot."""
    sensor_beta_v_per_e: np.ndarray | None = None
    """(2,) shift of the sensor potential per electron on (dot 1, dot 2), V (negative); None without a sensor dot."""

    @property
    def sensor_alpha_v_per_v(self) -> np.ndarray:
        """Sensor potential per volt of every gate: the basis at the sensor dot, or zeros."""
        if self.sensor_dot is None:
            return np.zeros(len(self.basis.gate_names))
        return self.sensor_dot.lever_arms


def mutual_capacitance_af(r1_nm: float, r2_nm: float, separation_nm: float, epsilon_r: float) -> float:
    """Two spheres of radii r1, r2 at centre separation s in a dielectric: C_m ≈ 4π ε r1 r2 / s (aF)."""
    s = max(separation_nm, r1_nm + r2_nm)
    return float(4 * np.pi * EPSILON_0_F_PER_M * epsilon_r * r1_nm * r2_nm / s * 1e-9 / 1e-18)


def total_capacitance_af(dot: DerivedDot, epsilon_r: float) -> float:
    """A dot's total capacitance: its disk self-capacitance over (1 − Σ_j α_j), so that the gate
    capacitances C·α are part of it (aF)."""
    return float(Dot(dot.name, dot.radius_nm, epsilon_r).self_capacitance_af / max(1.0 - float(np.sum(dot.lever_arms)), 0.05))


def sensor_coupling(c_dot_af: np.ndarray, dots: tuple[DerivedDot, DerivedDot], sensor: DerivedDot, epsilon_r: float) -> tuple[float, np.ndarray]:
    """(E_c,s in eV, β (2,) in V per electron) of the sensor dot from the three-dot capacitance matrix:
    the double dot's own 2×2 block, the sensor's total capacitance on the diagonal, and the sphere
    mutual capacitance between the sensor and each dot off it."""
    cm = np.array([mutual_capacitance_af(d.radius_nm, sensor.radius_nm, float(np.hypot(d.x_nm - sensor.x_nm, d.y_nm - sensor.y_nm)), epsilon_r) for d in dots])
    c = np.zeros((3, 3))
    c[:2, :2] = c_dot_af
    c[2, 2] = total_capacitance_af(sensor, epsilon_r)
    c[2, :2] = c[:2, 2] = -cm
    inv = np.linalg.inv(c)
    return float(E_UNIT_V_PER_E_AF * inv[2, 2]), -E_UNIT_V_PER_E_AF * inv[2, :2]


def derive_double_dot(
    layout: Layout,
    operating_point_v: np.ndarray,
    temperature_k: float = 0.15,
    pitch_nm: float = 2.0,
    basis: GateBasis | None = None,
) -> DerivedDoubleDot:
    """Find the two dots the layout forms at `operating_point_v` (one voltage per gate, in layout order)
    and build the CI model whose lever arms are the electrostatic ones. Compile-tier."""
    v = np.asarray(operating_point_v, dtype=float)
    if v.shape != (len(layout.gates),):
        raise ValueError(f"operating_point_v needs {len(layout.gates)} voltages, in layout gate order")
    mat = MaterialClass.by_name(layout.material)
    basis = basis or GateBasis.from_layout(layout, pitch_nm)
    d1, d2 = find_dots(basis, v, layout.channel_nm, 2, mat.effective_mass_me, names=("dot1", "dot2"))
    sensor = None
    if layout.sensor_nm is not None:
        sensor = find_dots(basis, v, layout.sensor_nm, 1, mat.effective_mass_me, names=("sensor",))[0]
    eps = layout.stack.epsilon_well
    alpha = np.stack([d1.lever_arms, d2.lever_arms])  # [dot, gate]
    # The CI diagonal is the dot's total capacitance: its own (disk) capacitance plus the gate
    # capacitances, and C_gate = C_dot·α, so C_ii = C_self / (1 − Σ_j α_ij).
    c_self = np.array([Dot("dot1", d1.radius_nm, eps).self_capacitance_af, Dot("dot2", d2.radius_nm, eps).self_capacitance_af])
    c_total = c_self / np.clip(1.0 - alpha.sum(axis=1), 0.05, None)
    cm = mutual_capacitance_af(d1.radius_nm, d2.radius_nm, np.hypot(d1.x_nm - d2.x_nm, d1.y_nm - d2.y_nm), eps)
    c_dot = np.array([[c_total[0], -cm], [-cm, c_total[1]]])
    ci = ConstantInteractionParams(
        gate_names=layout.gate_names,
        c_dot_af=c_dot,
        c_gate_af=c_dot @ alpha,
        temperature_k=temperature_k,
    )
    e_c, beta = sensor_coupling(c_dot, (d1, d2), sensor, eps) if sensor is not None else (None, None)
    return DerivedDoubleDot(basis, (d1, d2), sensor, barrier_between(basis, v, d1, d2), ci, v, e_c, beta)
