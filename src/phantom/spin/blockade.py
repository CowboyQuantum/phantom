"""Transport through the double dot at finite bias with Pauli spin blockade (ADR 0022).

At a finite source–drain bias a pair of bias triangles opens at each vertex of the honeycomb. Around the
(1,1)/(0,2) vertex and with electrons flowing from the left lead into dot 1 and out of dot 2 into the
right lead, the two cycles are

    electron triangle: (0,1) → (1,1) → (0,2) → (0,1)
    hole triangle:     (0,2) → (1,2) → (1,1) → (0,2)

and both pass through the interdot step (1,1) → (0,2). An electron loaded into (1,1) next to the one
already there forms a singlet one time in four and a triplet three times in four. The singlet relaxes
into (0,2); the triplet cannot, because the (0,2) ground state is a singlet and the (0,2) triplet lies
one valley splitting E_V higher (`valley`), so the current is blocked until the detuning exceeds E_V
or the triplet is mixed into the singlet (hyperfine or spin–orbit mixing, strongest at zero field where
the triplet states are degenerate with the singlet). Reversing the bias runs the cycles through
(0,2) → (1,1), where the pair always starts as a singlet, so the reverse direction carries no blockade.

Each cycle's current is e over the sum of the inverse rates of its steps, averaged over the singlet
and triplet branches: I = e / (¼ τ_S + ¾ τ_T). Lead steps carry Fermi factors at the electron
temperature; the interdot step is downhill-allowed with a thermal tail (a constant inelastic rate,
`calibration: none`). The two cycles add. This is the sequential-tunnelling picture, not a master
equation over all states (ADR 0017 covers the latter at zero bias).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phantom.charge.ci import ConstantInteractionParams
from phantom.constants import E_CHARGE_C, E_UNIT_V_PER_E_AF, K_B_EV_PER_K

MODEL_VERSION = 1

S01, S11, S02, S12 = (0, 1), (1, 1), (0, 2), (1, 2)


@dataclass(frozen=True)
class BlockadeRates:
    """Rates of the transport cycle (Hz) and the mixing that lifts the blockade. Placeholders."""

    gamma_lead_hz: float = 1e8
    """Dot-to-lead rate for loading and unloading."""
    gamma_interdot_hz: float = 1e7
    """Inelastic (1,1) → (0,2) relaxation rate once it is downhill."""
    gamma_mix_hz: float = 1e7
    """Triplet-to-singlet mixing rate at zero field (hyperfine and spin–orbit)."""
    b_mix_t: float = 1e-3
    """Field scale over which the mixing is quenched: Γ_mix(B) = Γ_mix / (1 + (B / b_mix_t)²)."""

    def mixing_hz(self, b_t: float) -> float:
        return self.gamma_mix_hz / (1.0 + (b_t / self.b_mix_t) ** 2)


def ci_energies_ev(ci: ConstantInteractionParams, v: np.ndarray, states: tuple[tuple[int, int], ...]) -> np.ndarray:
    """Constant-interaction energies (..., k) of the integer occupations `states` at gate voltages `v` (..., gates)."""
    n_g = ci.induced_charge(v)
    d = np.asarray(states, dtype=float) - n_g[..., None, :]
    c_inv = np.linalg.inv(ci.c_dot_af)
    return 0.5 * E_UNIT_V_PER_E_AF * np.einsum("...ki,ij,...kj->...k", d, c_inv, d)


def _fermi(x_ev: np.ndarray, kt_ev: float) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(np.clip(x_ev / kt_ev, -700, 700)))


def _cycle_current_a(steps_s: list[np.ndarray], steps_t: list[np.ndarray] | None) -> np.ndarray:
    """e / (¼ τ_S + ¾ τ_T) with τ = Σ 1/Γ over the cycle's steps; without a triplet branch, e / τ_S."""
    tau_s = sum(1.0 / np.maximum(g, 1e-300) for g in steps_s)
    if steps_t is None:
        return E_CHARGE_C / tau_s
    tau_t = sum(1.0 / np.maximum(g, 1e-300) for g in steps_t)
    return E_CHARGE_C / (0.25 * tau_s + 0.75 * tau_t)


def bias_triangle_current_a(
    ci: ConstantInteractionParams,
    v: np.ndarray,
    v_sd_v: float,
    e_valley_ev: float,
    b_t: float,
    rates: BlockadeRates | None = None,
    kt_ev: float | None = None,
) -> np.ndarray:
    """Current (..., ) through the pair at bias `v_sd_v` (positive: electrons enter dot 1 from the left
    lead at +eV/2 and leave dot 2 to the right lead at −eV/2) around the (1,1)/(0,2) vertex, with Pauli
    blockade lifted by the valley-excited (0,2) triplet beyond E_V and by field-dependent mixing."""
    rates = rates or BlockadeRates()
    kt = K_B_EV_PER_K * max(ci.temperature_k, 1e-6) if kt_ev is None else kt_ev
    e01, e11, e02, e12 = np.moveaxis(ci_energies_ev(ci, v, (S01, S11, S02, S12)), -1, 0)
    mu_l, mu_r = 0.5 * abs(v_sd_v), -0.5 * abs(v_sd_v)  # in eV for a unit charge
    gl, gm = rates.gamma_lead_hz, rates.gamma_interdot_hz
    eps = e11 - e02  # > 0: (1,1) above (0,2), the interdot step is downhill in the forward direction
    down = 1.0 / (1.0 + np.exp(np.clip(-eps / kt, -700, 700)))
    if v_sd_v >= 0:
        interdot_s = gm * down
        interdot_t = gm / (1.0 + np.exp(np.clip(-(eps - e_valley_ev) / kt, -700, 700))) + rates.mixing_hz(b_t)
        electron = _cycle_current_a(
            [gl * _fermi(e11 - e01 - mu_l, kt), interdot_s, gl * (1 - _fermi(e02 - e01 - mu_r, kt))],
            [gl * _fermi(e11 - e01 - mu_l, kt), interdot_t, gl * (1 - _fermi(e02 - e01 - mu_r, kt))],
        )
        hole = _cycle_current_a(
            [gl * _fermi(e12 - e02 - mu_l, kt), gl * (1 - _fermi(e12 - e11 - mu_r, kt)), interdot_s],
            [gl * _fermi(e12 - e02 - mu_l, kt), gl * (1 - _fermi(e12 - e11 - mu_r, kt)), interdot_t],
        )
        return electron + hole
    # reverse bias: the right lead sits at +eV/2 and loads dot 2; the pair starts in the (0,2) singlet
    up = 1.0 - down  # (0,2) → (1,1) is downhill when eps < 0
    interdot = gm * up
    electron = _cycle_current_a([gl * _fermi(e02 - e01 - mu_l, kt), interdot, gl * (1 - _fermi(e11 - e01 - mu_r, kt))], None)
    hole = _cycle_current_a([gl * _fermi(e12 - e11 - mu_l, kt), gl * (1 - _fermi(e12 - e02 - mu_r, kt)), interdot], None)
    return electron + hole
