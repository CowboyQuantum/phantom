"""phantom.compile — The offline tier: layout + material + catalog draw + seed → a serialisable, versioned `Device`. May import every layer; only the consumer surface (`phantom.api`, `phantom.contract`) imports it (ADR 0024).

Modules: `device` (`CompileParams`, `compile_device`; the record itself is `phantom.device.Device`),
`backends` (closed-form electrostatics; the FEM hook behind the `[fem]` extra).
"""

from phantom.compile.backends import Backend, basis_for, closed_form_basis, fem_basis
from phantom.compile.device import (
    CompileParams,
    compile_device,
    interdot_tunnel_coupling_ev,
    replace_params,
)
from phantom.device import DEVICE_FORMAT_VERSION, CompiledBarrier, CompiledDot, Device

MODEL_VERSION = 1

__all__ = [
    "DEVICE_FORMAT_VERSION",
    "MODEL_VERSION",
    "Backend",
    "CompileParams",
    "CompiledBarrier",
    "CompiledDot",
    "Device",
    "basis_for",
    "closed_form_basis",
    "compile_device",
    "fem_basis",
    "interdot_tunnel_coupling_ev",
    "replace_params",
]
