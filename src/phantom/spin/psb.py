"""Pauli spin blockade readout with the valley draw (ADR 0021).

A pair separated in (1,1) is pulsed toward (0,2) for readout. The singlet tunnels into (0,2) once the
detuning ε passes the anticrossing (ε = 0 here, broadened by the tunnel coupling), because the (0,2)
ground state is a singlet in the lower valley. A triplet is blocked: the lowest (0,2) triplet puts one
electron in the upper valley of the receiving dot and lies E_V higher, so it needs ε > E_V. The charge
sensor therefore tells singlet from triplet only inside the window 0 < ε < E_V, and the window is the
valley splitting of the receiving dot, drawn per device (`valley`). A triplet that relaxes to the
singlet during the integration (T₁ at the read point, with the hot spot of `valley.signatures` when
the Zeeman energy sits near E_V) is read as a singlet.

Charge-readout error from the sensor (Step 4, `readout.single_shot`) is a separate factor that the
consumer multiplies in; here the visibility is the spin-to-charge part only.
"""

from __future__ import annotations

import numpy as np

MODEL_VERSION = 1


def charge_step(eps_ev: np.ndarray, threshold_ev: float, t_ev: float, kt_ev: float) -> np.ndarray:
    """Mean (0,2) occupation across an anticrossing at `threshold_ev` with coupling t (matrix element
    √2 t) at temperature kT: the ground-state weight ½(1 + x/√(x² + 8t²)) thermally averaged with the
    excited branch."""
    x = np.asarray(eps_ev, dtype=float) - threshold_ev
    w = np.sqrt(x**2 + 8 * t_ev**2)
    ground = 0.5 * (1 + x / w)
    if kt_ev <= 0:
        return ground
    p_exc = 1.0 / (1.0 + np.exp(np.minimum(w / kt_ev, 700.0)))
    return ground * (1 - p_exc) + (1 - ground) * p_exc


def psb_occupation(eps_ev: np.ndarray, e_valley_ev: float, t_ev: float, kt_ev: float) -> tuple[np.ndarray, np.ndarray]:
    """(0,2) occupation against detuning for a singlet and for a triplet: the singlet steps at ε = 0, the
    triplet at ε = E_V (one electron in the upper valley of the receiving dot)."""
    return charge_step(eps_ev, 0.0, t_ev, kt_ev), charge_step(eps_ev, e_valley_ev, t_ev, kt_ev)


def psb_contrast(eps_ev: np.ndarray, e_valley_ev: float, t_ev: float, kt_ev: float) -> np.ndarray:
    """Singlet-minus-triplet (0,2) occupation: the spin-to-charge contrast, nonzero inside the window."""
    s, t = psb_occupation(eps_ev, e_valley_ev, t_ev, kt_ev)
    return s - t


def readout_visibility(contrast: np.ndarray, t_integration_s: float, t1_read_s: float) -> np.ndarray:
    """Visibility F_S + F_T − 1 of the spin-to-charge map: the contrast times the probability that a
    triplet survives the integration, exp(−t_int / T₁)."""
    return np.asarray(contrast) * np.exp(-t_integration_s / t1_read_s)


def best_read_point(eps_ev: np.ndarray, e_valley_ev: float, t_ev: float, kt_ev: float) -> tuple[float, float]:
    """(ε, contrast) at the maximum of the contrast over `eps_ev`."""
    c = psb_contrast(eps_ev, e_valley_ev, t_ev, kt_ev)
    k = int(np.argmax(c))
    return float(eps_ev[k]), float(c[k])
