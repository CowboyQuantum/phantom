"""phantom.spin — Few-qubit spin dynamics under pulses: Zeeman with gate-dependent g, exchange from `charge`, driving, hyperfine, Lindblad with dephasing from charge noise; Pauli spin blockade readout with the valley draw.

Modules: `exchange` (J from the Hubbard two-electron sector, the detuning axis, dJ/dε, quasi-static
averaging; J of the valence electrons in a compiled device's CI window, #60), `dynamics` (rotating-frame Hamiltonians, Lindblad propagator, Rabi and Ramsey limits, the
spin draw, EDSR drive strength), `psb` (singlet and triplet charge steps, the readout window set by the
valley splitting, visibility), `blockade` (bias triangles with Pauli blockade lifted by the valley
splitting and by mixing), `mixing` (the spin funnel), `elzerman` (spin-dependent tunnelling shots),
`noise_spectroscopy` (filter functions, CPMG coherence on the draw's spectrum), `gates` (exchange-gate
fidelity with the draw's noise).
"""

from phantom.spin.blockade import BlockadeRates, bias_triangle_current_a, ci_energies_ev
from phantom.spin.dynamics import (
    H_EV_S,
    MU_B_EV_PER_T,
    Lindblad,
    SpinDraw,
    SpinParams,
    draw_spin,
    drive_displacement_nm,
    edsr_rabi_energy_ev,
    probability,
    quasi_static_sigma_ev,
    rabi_probability,
    ramsey_probability,
    rotating_frame_hamiltonian_ev,
    singlet_return_probability,
    singlet_state,
    t2_star_s,
)
from phantom.spin.elzerman import elzerman_shots, empty_probability
from phantom.spin.exchange import (
    detuning_noise_sigma_ev,
    detuning_voltages,
    exchange_j_ev,
    exchange_second_order_ev,
    exchange_sensitivity_ev_per_ev,
    local_exchange_ev,
    quasi_static_average,
    singlet_exchange_ev,
    symmetric_point_v,
)
from phantom.spin.gates import average_gate_fidelity, entanglement_fidelity, swap_gate_fidelity
from phantom.spin.mixing import funnel_singlet_probability, two_level_transfer
from phantom.spin.noise_spectroscopy import (
    coherence,
    cpmg_filter,
    filter_frequency_hz,
    log_grid_hz,
    omega_spectrum,
    one_over_f_amplitude_sq,
    one_over_f_psd,
)
from phantom.spin.psb import best_read_point, psb_contrast, psb_occupation, readout_visibility

MODEL_VERSION = 1

__all__ = [
    "H_EV_S",
    "MODEL_VERSION",
    "MU_B_EV_PER_T",
    "BlockadeRates",
    "Lindblad",
    "SpinDraw",
    "SpinParams",
    "average_gate_fidelity",
    "best_read_point",
    "bias_triangle_current_a",
    "ci_energies_ev",
    "coherence",
    "cpmg_filter",
    "detuning_noise_sigma_ev",
    "detuning_voltages",
    "draw_spin",
    "drive_displacement_nm",
    "edsr_rabi_energy_ev",
    "elzerman_shots",
    "empty_probability",
    "entanglement_fidelity",
    "exchange_j_ev",
    "exchange_second_order_ev",
    "exchange_sensitivity_ev_per_ev",
    "filter_frequency_hz",
    "funnel_singlet_probability",
    "local_exchange_ev",
    "log_grid_hz",
    "omega_spectrum",
    "one_over_f_amplitude_sq",
    "one_over_f_psd",
    "probability",
    "psb_contrast",
    "psb_occupation",
    "quasi_static_average",
    "quasi_static_sigma_ev",
    "rabi_probability",
    "ramsey_probability",
    "readout_visibility",
    "rotating_frame_hamiltonian_ev",
    "singlet_exchange_ev",
    "singlet_return_probability",
    "singlet_state",
    "swap_gate_fidelity",
    "symmetric_point_v",
    "t2_star_s",
    "two_level_transfer",
]
