"""phantom.optimum — For a compiled device, the best achievable: sensor operating point, lattice (virtual-gate matrix), interdot coupling targets, and with `spin`, readout fidelity and coherence. Consumers that hide truth do not forward this.

Module: `targets` (`optimum(device)` and its parts: sensor flank, virtual-gate matrix, barrier voltage
for a target coupling, the PSB read point and Elzerman peak time, the exchange sweet spot with its T₂*
and SWAP fidelity).
"""

from phantom.optimum.targets import (
    CoherenceOptimum,
    CouplingOptimum,
    LatticeOptimum,
    Optimum,
    ReadoutOptimum,
    SensorOptimum,
    coherence_optimum,
    coupling_optimum,
    lattice_optimum,
    optimum,
    readout_optimum,
    sensor_optimum,
    virtual_gate_matrix,
)

MODEL_VERSION = 1

__all__ = [
    "MODEL_VERSION",
    "CoherenceOptimum",
    "CouplingOptimum",
    "LatticeOptimum",
    "Optimum",
    "ReadoutOptimum",
    "SensorOptimum",
    "coherence_optimum",
    "coupling_optimum",
    "lattice_optimum",
    "optimum",
    "readout_optimum",
    "sensor_optimum",
    "virtual_gate_matrix",
]
