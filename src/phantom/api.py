"""phantom.api — the one import surface for consumers (plan §4, ADR 0003, 0024; documented in docs/CONTRACT.md).

A consumer imports from here and nowhere else in phantom; every name below is covered by the API
snapshot (`docs/api/surface.json`) and by the contract test (`python -m phantom.contract`). A change
to any of them is listed in `CHANGELOG.md`, and a breaking change is a breaking release.

    from phantom.api import CompileParams, Device, compile_device, optimum, two_dot_layout

    device = compile_device(two_dot_layout("SiGe", 80.0), CompileParams(operating_point_v=...), seed=2)
    device.to_json("device.json")                 # 8 KB, canonical, diffable; Device.from_json reloads it
    state = device.at(v)                          # noise-free, vectorised over leading axes
    run = device.session(noise_seed=7)            # a clock, the device's noise, the dot charge
    run.hold(v_pulse, 1e-6); m = run.sense(v, dwell_s=1e-3, samples=100, v_sd=50e-6)
    device.knowable()                             # what a user may be shown
    device.truth(), optimum(device)               # truth: a consumer that hides it does not forward it

Nothing in a layer imports this module; it sits above `compile` with `phantom.contract`.
"""

from __future__ import annotations

import importlib
import pkgutil

import phantom
from phantom import __version__
from phantom.compile import CompileParams, compile_device
from phantom.device import DEVICE_FORMAT_VERSION, CompiledBarrier, CompiledDot, Device
from phantom.device.session import Measurement, Session, SpinOutcome, State
from phantom.device.truth import Knowable, SensorTruth, Truth
from phantom.layout import Gate, Layer, Layout, Stack, two_dot_layout
from phantom.optimum import (
    CoherenceOptimum,
    CouplingOptimum,
    LatticeOptimum,
    Optimum,
    ReadoutOptimum,
    SensorOptimum,
    optimum,
)
from phantom.spin import SpinParams
from phantom.valley import ValleyParams


def model_versions() -> dict[str, int]:
    """`MODEL_VERSION` of every phantom module that defines one, by dotted name relative to `phantom`
    (the runtime models a measurement used, beside the compile-time ones in `Device.provenance`)."""
    out: dict[str, int] = {}
    for info in pkgutil.walk_packages(phantom.__path__, "phantom."):
        if info.name in ("phantom.api", "phantom.contract", "phantom.__main__"):
            continue
        module = importlib.import_module(info.name)
        version = getattr(module, "MODEL_VERSION", None)
        if isinstance(version, int):
            out[info.name.removeprefix("phantom.")] = version
    return dict(sorted(out.items()))


__all__ = [
    "DEVICE_FORMAT_VERSION",
    "CoherenceOptimum",
    "CompileParams",
    "CompiledBarrier",
    "CompiledDot",
    "CouplingOptimum",
    "Device",
    "Gate",
    "Knowable",
    "LatticeOptimum",
    "Layer",
    "Layout",
    "Measurement",
    "Optimum",
    "ReadoutOptimum",
    "SensorOptimum",
    "SensorTruth",
    "Session",
    "SpinOutcome",
    "SpinParams",
    "Stack",
    "State",
    "Truth",
    "ValleyParams",
    "__version__",
    "compile_device",
    "model_versions",
    "optimum",
    "two_dot_layout",
]
