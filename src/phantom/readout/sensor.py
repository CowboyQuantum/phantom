"""A charge sensor: a single-electron transistor (SET) read in sequential tunnelling.

Seeded from qdt-sim's `DoubleDotModel` sensor (Step 1). The sensor is itself a single quantum dot
with charging energy `E_c = charging_energy_ev` and lever arm `α_s = E_c / period_v` (eV per V of its
potential φ_s), biased by a source–drain voltage `V_sd`. Its levels are `μ_k = E_c·k − α_s·(φ_s − phase_v)`;
sequential tunnelling at the sensor's thermal energy `kT_s = α_s·width_v / 2` gives

    I = I_bg + I_0 · Σ_k [ f((μ_k − eV_sd/2)/kT_s) − f((μ_k + eV_sd/2)/kT_s) ],   f(x) = 1/(1+eˣ),

a periodic train of thermally broadened Coulomb peaks at small bias and Coulomb diamonds of height
`E_c/e` against `V_sd`. `I_0 = 4 kT_s · peak_a / (e · bias_v)` makes `peak_a` the linear-response peak
current at the default bias. Because the sensor *is* a single dot in transport, the same functions
serve the single-dot protocols (`coulomb-oscillations-single-dot`, `diamonds-single-dot`).

The sensor's potential is linear in the gate voltages and in the charge on the dots it watches:

    φ_s = α · V + β · ⟨N⟩ + δφ                                [V-equivalent]

so every dot charge transition shows as a step in sensor current.

**A sensor with leads** (ADR 0026). A compiled sensor dot has two tunnel barriers and a bath behind
each (`SensorLeads`): the rates Γ_S, Γ_D through the barriers follow the gate voltages through their
compiled exponential and saturate at the dot's attempt rate, and each bath conducts only where its
electron energy lies below the accumulation threshold. Its current is then the lifetime-broadened
resonant-tunnelling sum of `readout.broadening`,

    I = e Γ_SΓ_D/(Γ_S + Γ_D) · Σ_k [ n_γ(μ_k − eV_sd/2) − n_γ(μ_k + eV_sd/2) ],   γ = ħ(Γ_S + Γ_D),

so `peak_a` and `background_a` are not used: the peak height follows the barriers, opening them
raises the peaks until they broaden towards the e²/h bound, and the current between peaks is the
Lorentzian tails of the nearest levels (elastic cotunnelling through one level). Without leads the
sensor is the typed SET above, unchanged.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.special import expit

from phantom.constants import E_CHARGE_C, HBAR_EV_S
from phantom.readout.broadening import (
    lineshape_per_ev,
    lineshape_slope_per_ev2,
    occupation,
    wing_ev,
)
from phantom.readout.rates import saturating_rates_hz

MODEL_VERSION = 1
"""Bumped whenever any output of this module changes for the same parameters (ADR 0003); a sensor
with leads is new and leaves the output of a sensor without them unchanged."""

WING_WIDTHS = math.atanh(1.0 / math.sqrt(3.0))
"""Distance of the wing (the inflection of sech²) from a thermal peak, in peak widths."""
WING_SLOPE_FACTOR = 4.0 / (3.0 * math.sqrt(3.0))
"""|d sech²(x)/dx| at the wing."""


@dataclass(frozen=True)
class SensorLeads:
    """The sensor dot's two tunnel barriers and the bath behind each, compiled at a reference gate
    point (ADR 0026). Lead 0 is the source (the left barrier and bath), lead 1 the drain.

    The rates are the WKB transmissions through the barriers at the compiled tunnelling level times the
    dot's attempt rate (`readout.rates`, ADR 0023). At other gate voltages the exponent moves linearly
    with the compiled slopes d ln Γ/dV and a logistic saturates the rate at the attempt rate, anchored
    so the rate at the reference point is the compiled one:

        Γ(V) = f_a · σ( ln(Γ_ref / (f_a − Γ_ref)) + Σ_g s_g (V_g − V_ref,g) ),   σ(x) = 1/(1 + e⁻ˣ).

    A bath is the two-dimensional electron gas under a reservoir gate; it conducts with the occupation
    σ((U_th − U_bath)/w) of the soft accumulation step at the turn-on protocols' threshold, U_bath the
    electron energy at the gate's centroid, linear in the gate voltages. Calibration: none."""

    barrier_gates: tuple[str, str]
    reservoir_gates: tuple[str, str]
    v_ref_v: tuple[float, ...]
    """Carried gate voltages the rates are compiled at, one per gate, V."""
    gamma_ref_hz: tuple[float, float]
    """Tunnel rate through each barrier at `v_ref_v`, Hz."""
    attempt_rate_hz: float
    """The sensor dot's attempt frequency, Hz: the ceiling of each rate as its barrier opens."""
    log_rate_slopes_per_v: tuple[tuple[float, ...], tuple[float, ...]]
    """d ln Γ/dV of each lead with respect to every gate at the compiled tunnelling level, 1/V."""
    bath_lever_arms: tuple[tuple[float, ...], tuple[float, ...]]
    """Potential per volt of every gate at each bath, dimensionless."""
    bath_offset_ev: tuple[float, float]
    """Gate-independent electron energy at each bath (the draw's interface charges), eV."""
    bath_threshold_ev: float = -0.05
    """Electron energy below which a bath is accumulated, eV (the turn-on protocols' threshold)."""
    bath_width_ev: float = 0.004
    """Width of the accumulation step, eV."""
    tunnelling_level_ev: tuple[float, float] = (0.0, 0.0)
    """The energy each lead's electrons tunnel at, held fixed over gate voltage, eV (diagnostic and for
    recomputing a rate at the compile tier)."""
    barrier_height_ev: tuple[float, float] = (0.0, 0.0)
    """Barrier top above the tunnelling level at `v_ref_v`, eV (diagnostic)."""
    barrier_width_nm: tuple[float, float] = (0.0, 0.0)
    """Forbidden length at the tunnelling level at `v_ref_v`, nm (diagnostic)."""

    def __post_init__(self) -> None:
        n = len(self.v_ref_v)
        if any(len(row) != n for row in self.log_rate_slopes_per_v + self.bath_lever_arms):
            raise ValueError("slopes and bath lever arms need one entry per gate")
        if not self.attempt_rate_hz > 0 or self.bath_width_ev <= 0:
            raise ValueError("attempt_rate_hz and bath_width_ev must be positive")

    @classmethod
    def from_dict(cls, d: dict) -> SensorLeads:
        pairs = ("gamma_ref_hz", "bath_offset_ev", "tunnelling_level_ev", "barrier_height_ev", "barrier_width_nm", "barrier_gates", "reservoir_gates")
        rows = ("log_rate_slopes_per_v", "bath_lever_arms")
        return cls(**{
            **d,
            **{k: tuple(d[k]) for k in pairs if k in d},
            **{k: tuple(tuple(float(x) for x in row) for row in d[k]) for k in rows},
            "v_ref_v": tuple(float(x) for x in d["v_ref_v"]),
        })

    def _dv(self, v: np.ndarray | None) -> np.ndarray:
        if v is None:
            return np.zeros(len(self.v_ref_v))
        return np.asarray(v, dtype=float) - np.asarray(self.v_ref_v)

    def bath_energy_ev(self, v: np.ndarray | None = None) -> np.ndarray:
        """Electron energy at each bath (..., 2), eV, for carried voltages `v` (..., gates)."""
        v = np.asarray(self.v_ref_v) if v is None else np.asarray(v, dtype=float)
        return -(v @ np.asarray(self.bath_lever_arms).T) + np.asarray(self.bath_offset_ev)

    def bath_occupation(self, v: np.ndarray | None = None) -> np.ndarray:
        """How far each bath is accumulated (..., 2): 1 well above its turn-on, 0 depleted."""
        return expit((self.bath_threshold_ev - self.bath_energy_ev(v)) / self.bath_width_ev)

    def barrier_rates_hz(self, v: np.ndarray | None = None) -> np.ndarray:
        """Γ through each barrier (..., 2), Hz, with the baths full."""
        return saturating_rates_hz(np.asarray(self.gamma_ref_hz, dtype=float), self.attempt_rate_hz, self._dv(v) @ np.asarray(self.log_rate_slopes_per_v).T)

    def rates_hz(self, v: np.ndarray | None = None) -> np.ndarray:
        """Γ_S, Γ_D (..., 2), Hz: the barriers' rates times the baths' occupations."""
        return self.barrier_rates_hz(v) * self.bath_occupation(v)


@dataclass(frozen=True)
class SensorSETParams:
    """Parameters of the SET charge sensor and its coupling to the dots and gates.

    Defaults inherited from qdt-sim; a plausible Si-class device, not a measured one. Ranges: charging
    energy 0.3–3 meV, peak period 0.05–0.5 V, thermal width 5–50 mV, peak current 0.1–10 nA. A compiled
    sensor dot (`device.compiled_sensor`) sets φ_s to the potential at the dot in volts (lever arm
    1 eV/V: `period_v` = its charging energy, `width_v` = 2kT), β from the layout's capacitances, and
    `leads`; its peaks are then a few to tens of mV apart on its plunger (ADR 0026).
    """

    alpha_v_per_v: np.ndarray = field(default_factory=lambda: np.array([0.01, 0.01, 1.0]))
    """Sensor lever arm per gate: dφ_s/dV_gate, one entry per gate of the charge model."""
    beta_v_per_e: np.ndarray = field(default_factory=lambda: np.array([-0.004, -0.003]))
    """Sensor potential shift per electron on each dot (V-equivalent)."""
    period_v: float = 0.25
    """Spacing of the sensor's Coulomb peaks in φ_s (V)."""
    width_v: float = 0.04
    """Thermal width of a sensor peak in φ_s (V)."""
    peak_a: float = 1.0e-9
    """Linear-response peak current (A) at `bias_v`."""
    background_a: float = 0.05e-9
    """Sensor current far from a peak (A)."""
    phase_v: float = -0.06
    """φ_s of one sensor peak (V); the rest follow at `period_v`."""
    charging_energy_ev: float = 1.0e-3
    """Charging energy of the sensor SET (eV): e times the bias at a Coulomb diamond's apex."""
    bias_v: float = 50e-6
    """Source–drain bias the sensor is read at unless one is given (V)."""
    leads: SensorLeads | None = None
    """The sensor's barriers and baths (a compiled sensor dot); None for the typed SET."""

    def __post_init__(self) -> None:
        if self.beta_v_per_e.shape != (2,):
            raise ValueError("beta_v_per_e must be (2,)")
        if self.period_v <= 0 or self.width_v <= 0 or self.charging_energy_ev <= 0 or self.bias_v <= 0:
            raise ValueError("period_v, width_v, charging_energy_ev and bias_v must be positive")

    @property
    def lever_arm_ev_per_v(self) -> float:
        """α_s: sensor level shift per volt of φ_s (eV/V); one period in φ_s adds one electron."""
        return self.charging_energy_ev / self.period_v

    @property
    def kt_ev(self) -> float:
        """Thermal energy of the sensor SET (eV), set so its zero-bias peaks have width `width_v`."""
        return self.lever_arm_ev_per_v * self.width_v / 2

    def potential(self, v: np.ndarray, n: np.ndarray, delta_phi_v: float | np.ndarray = 0.0) -> np.ndarray:
        """φ_s (V-equivalent) for gate voltages `v` (..., gates) and dot occupation `n` (..., 2)."""
        return (
            np.asarray(v, dtype=float) @ self.alpha_v_per_v
            + np.asarray(n, dtype=float) @ self.beta_v_per_e
            + delta_phi_v
        )

    def levels_ev(self, phi: np.ndarray) -> np.ndarray:
        """μ_k (..., 5) in eV: the five sensor levels nearest the bias window at potential `phi`."""
        x = (np.asarray(phi, dtype=float) - self.phase_v) / self.period_v
        k = np.rint(x)[..., None] + np.arange(-2, 3)
        return self.charging_energy_ev * (k - x[..., None])

    def current_a(self, phi: np.ndarray, v_sd: float | np.ndarray | None = None, v: np.ndarray | None = None) -> np.ndarray:
        """Sensor current (A) at potential `phi` and source–drain bias `v_sd` (V; default `bias_v`).
        A sensor with leads takes its rates at carried gate voltages `v` (..., gates), at the compiled
        reference point when None; the typed sensor ignores `v`."""
        v_sd = self.bias_v if v_sd is None else np.asarray(v_sd, dtype=float)
        kt = self.kt_ev
        mu = self.levels_ev(phi)
        half = np.asarray(v_sd)[..., None] / 2  # eV numerically, for V_sd in V
        if self.leads is not None:
            gamma_eff, gamma_ev = self._coupling(v)
            window = occupation(mu - half, gamma_ev[..., None], kt) - occupation(mu + half, gamma_ev[..., None], kt)
            # far from every level both occupations round to 0 or 1 and their difference to ±1e-16:
            # the current keeps the sign of the bias
            window = np.where(half >= 0, np.maximum(window, 0.0), np.minimum(window, 0.0))
            return E_CHARGE_C * gamma_eff * window.sum(axis=-1)
        with np.errstate(over="ignore"):
            window = 1 / (1 + np.exp((mu - half) / kt)) - 1 / (1 + np.exp((mu + half) / kt))
        i_0 = 4 * kt * self.peak_a / self.bias_v
        return self.background_a + i_0 * window.sum(axis=-1)

    # --- the leads ---------------------------------------------------------------------------------

    def rates_hz(self, v: np.ndarray | None = None) -> np.ndarray | None:
        """Γ_S, Γ_D (..., 2) at carried voltages `v`, Hz; None for the typed sensor."""
        return None if self.leads is None else self.leads.rates_hz(v)

    def _coupling(self, v: np.ndarray | None) -> tuple[np.ndarray, np.ndarray]:
        """(e-less current scale Γ_SΓ_D/(Γ_S + Γ_D) in Hz, level width ħ(Γ_S + Γ_D) in eV) at `v`."""
        rates = self.leads.rates_hz(v)
        gs, gd = rates[..., 0], rates[..., 1]
        total = gs + gd
        return gs * gd / np.maximum(total, 1e-300), HBAR_EV_S * total

    def lifetime_width_v(self, v: np.ndarray | None = None) -> np.ndarray:
        """ħ(Γ_S + Γ_D) in φ_s volts: the Lorentzian full width of a level; 0 for the typed sensor."""
        if self.leads is None:
            return np.zeros(np.shape(v)[:-1] if v is not None else ())
        return self._coupling(v)[1] / self.lever_arm_ev_per_v

    def linear_response_peak_a(self, v_sd: float | None = None, v: np.ndarray | None = None) -> float:
        """Peak height of a Coulomb peak in linear response at bias `v_sd` (V), A: the typed sensor's
        `peak_a` scaled to the bias, or e·Γ_SΓ_D/(Γ_S + Γ_D)·eV_sd·L(0) with the broadened lineshape."""
        bias = self.bias_v if v_sd is None else float(v_sd)
        if self.leads is None:
            return float(self.peak_a * bias / self.bias_v)
        gamma_eff, gamma_ev = self._coupling(v)
        return float(E_CHARGE_C * gamma_eff * bias * lineshape_per_ev(0.0, gamma_ev, self.kt_ev))

    def wing_offset_v(self, v: np.ndarray | None = None) -> float:
        """φ_s from a peak to its wing, where the linear-response slope is steepest, V: the inflection
        of sech², `width_v`·atanh(1/√3), for the typed sensor; of the broadened lineshape with leads."""
        if self.leads is None:
            return float(WING_WIDTHS * self.width_v)
        return wing_ev(float(self._coupling(v)[1]), self.kt_ev) / self.lever_arm_ev_per_v

    def wing_slope_a_per_v(self, v_sd: float | None = None, v: np.ndarray | None = None) -> float:
        """|dI/dφ_s| at the wing in linear response at bias `v_sd`, A per V of φ_s."""
        bias = self.bias_v if v_sd is None else float(v_sd)
        if self.leads is None:
            return float(WING_SLOPE_FACTOR * self.linear_response_peak_a(bias) / self.width_v)
        gamma_eff, gamma_ev = self._coupling(v)
        eps = self.wing_offset_v(v) * self.lever_arm_ev_per_v
        slope = abs(float(lineshape_slope_per_ev2(eps, gamma_ev, self.kt_ev)))
        return float(E_CHARGE_C * gamma_eff * bias * slope * self.lever_arm_ev_per_v)

    def diamond_edges(self, k: int) -> tuple[tuple[float, float], tuple[float, float]]:
        """The two Coulomb-diamond edge lines through the k-th zero-bias degeneracy point, as
        (φ_s intercept at V_sd = 0, slope dV_sd/dφ_s): the rising edge and the falling one.
        An edge is where level k crosses a lead: μ_k = ±eV_sd/2 → V_sd = ±2 α_s (φ_s − φ_k)."""
        phi_k = self.phase_v + k * self.period_v
        s = 2 * self.lever_arm_ev_per_v
        return (phi_k, s), (phi_k, -s)

    def peak_positions_v(self, phi_min: float, phi_max: float) -> np.ndarray:
        """φ_s of every zero-bias Coulomb peak in [phi_min, phi_max]."""
        k0 = np.ceil((phi_min - self.phase_v) / self.period_v)
        k1 = np.floor((phi_max - self.phase_v) / self.period_v)
        return self.phase_v + np.arange(k0, k1 + 1) * self.period_v
