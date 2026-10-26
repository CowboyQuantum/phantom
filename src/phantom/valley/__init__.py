"""phantom.valley — Valley splitting per dot in Si as a statistical draw for the material and dot size. Atomistic models are out of scope.

Modules: `splitting` (Rice distribution from alloy disorder with the stated scalings; MOS deterministic
part; per-dot draws; interface-step hook), `signatures` (magnetospectroscopy kink, relaxation hot spot,
pulsed-gate spectroscopy).
"""

from phantom.valley.signatures import (
    MU_B_UEV_PER_T,
    kink_field_t,
    pulsed_spectroscopy_occupation,
    relaxation_rate_hz,
    second_electron_shift_uev,
)
from phantom.valley.splitting import ValleyParams, draw_valley_splitting_uev, rice_mean

MODEL_VERSION = 1

__all__ = [
    "MODEL_VERSION",
    "MU_B_UEV_PER_T",
    "ValleyParams",
    "draw_valley_splitting_uev",
    "kink_field_t",
    "pulsed_spectroscopy_occupation",
    "relaxation_rate_hz",
    "rice_mean",
    "second_electron_shift_uev",
]
