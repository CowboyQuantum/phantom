"""Dots and barriers as the named parts of a charge model.

A `Dot` carries the geometry that sets its charging energy; a `Barrier` carries the exponential
dependence of the tunnel coupling on its gate voltage. Both are parameters, not solvers: Step 2 derives
them from a layout, Step 1 types them in.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phantom.charge.ci import E_CHARGE_C

EPSILON_0_F_PER_M = 8.8541878128e-12
"""Vacuum permittivity in F/m (CODATA 2018)."""


@dataclass(frozen=True)
class Dot:
    """A quantum dot as a conducting disk of radius `radius_nm` in a dielectric `epsilon_r`.

    The self-capacitance of a thin disk in a uniform dielectric is C = 8 ε r (textbook electrostatics),
    which gives charging energies of a few meV for 20–40 nm dots in silicon (ε_r = 11.7). Range: radius
    10–60 nm. A reference entry replaces this closed form when Step 2 derives dots from a layout.
    """

    name: str
    radius_nm: float = 25.0
    epsilon_r: float = 11.7

    def __post_init__(self) -> None:
        if self.radius_nm <= 0 or self.epsilon_r <= 0:
            raise ValueError("radius_nm and epsilon_r must be positive")

    @property
    def self_capacitance_af(self) -> float:
        """C = 8 ε₀ ε_r r, in aF."""
        return 8 * EPSILON_0_F_PER_M * self.epsilon_r * self.radius_nm * 1e-9 / 1e-18

    @property
    def charging_energy_ev(self) -> float:
        """E_c = e² / C, in eV."""
        return E_CHARGE_C / (self.self_capacitance_af * 1e-18)


@dataclass(frozen=True)
class Barrier:
    """A tunnel barrier whose coupling falls exponentially as its gate closes.

    t(V_B) = t0 · exp((V_B − v_ref_v) / v_scale_v): `t0_ev` is the coupling at `v_ref_v`, and `v_scale_v`
    is the gate voltage that changes it by a factor e. Range: t0 1–200 µeV, v_scale 10–100 mV. The
    exponential form is the WKB limit for a smooth barrier; Step 2 fits it to the potential.
    """

    name: str
    t0_ev: float = 20e-6
    v_ref_v: float = 0.0
    v_scale_v: float = 0.05

    def __post_init__(self) -> None:
        if self.t0_ev < 0 or self.v_scale_v <= 0:
            raise ValueError("t0_ev must be >= 0 and v_scale_v > 0")

    def tunnel_coupling_ev(self, v_b: float | np.ndarray) -> np.ndarray:
        """t(V_B) in eV for barrier gate voltage `v_b` (V)."""
        return self.t0_ev * np.exp((np.asarray(v_b, dtype=float) - self.v_ref_v) / self.v_scale_v)
