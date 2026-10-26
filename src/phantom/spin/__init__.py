"""phantom.spin — Few-qubit spin dynamics under pulses: Zeeman with gate-dependent g, exchange from `charge`, driving, hyperfine, Lindblad with dephasing from charge noise; Pauli spin blockade readout with the valley draw.

Modules: `exchange` (J from the Hubbard two-electron sector, the detuning axis, dJ/dε, quasi-static
averaging), `dynamics` (rotating-frame Hamiltonians, Lindblad propagator, Rabi and Ramsey limits, the
spin draw), `psb` (singlet and triplet charge steps, the readout window set by the valley splitting,
visibility).
"""

from phantom.spin.dynamics import (
    H_EV_S,
    MU_B_EV_PER_T,
    Lindblad,
    SpinDraw,
    SpinParams,
    draw_spin,
    probability,
    quasi_static_sigma_ev,
    rabi_probability,
    ramsey_probability,
    rotating_frame_hamiltonian_ev,
    singlet_return_probability,
    singlet_state,
    t2_star_s,
)
from phantom.spin.exchange import (
    detuning_noise_sigma_ev,
    detuning_voltages,
    exchange_j_ev,
    exchange_second_order_ev,
    exchange_sensitivity_ev_per_ev,
    quasi_static_average,
    symmetric_point_v,
)
from phantom.spin.psb import best_read_point, psb_contrast, psb_occupation, readout_visibility

MODEL_VERSION = 1

__all__ = [
    "H_EV_S",
    "MODEL_VERSION",
    "MU_B_EV_PER_T",
    "Lindblad",
    "SpinDraw",
    "SpinParams",
    "best_read_point",
    "detuning_noise_sigma_ev",
    "detuning_voltages",
    "draw_spin",
    "exchange_j_ev",
    "exchange_second_order_ev",
    "exchange_sensitivity_ev_per_ev",
    "probability",
    "psb_contrast",
    "psb_occupation",
    "quasi_static_average",
    "quasi_static_sigma_ev",
    "rabi_probability",
    "ramsey_probability",
    "readout_visibility",
    "rotating_frame_hamiltonian_ev",
    "singlet_return_probability",
    "singlet_state",
    "symmetric_point_v",
    "t2_star_s",
]
