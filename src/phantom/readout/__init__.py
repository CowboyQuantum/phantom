"""phantom.readout — Lead tunnel rates, master equation, latching, charge-sensor (SET) and RF response, noise processes (white, telegraph, 1/f, fluctuator ensembles) with explicit time steps. Seeded from qdt-sim's noise models (Step 1).

Modules: `sensor` (SET sensor and single-dot transport current), `noise` (seeded processes with a time
axis), `measure` (a CI double dot read by the sensor with noise), `lockin` (the `sensor_lockin`
hardware-derivative channel), `transport` (phenomenological turn-on and pinch-off shapes, Step 1 only), `loading` (rate-limited
loading through a lead barrier in a raster scan; Step 4 replaces it with a master equation).
"""

from phantom.readout.loading import LeadBarrier, LeadRateLandscape, latched_occupation
from phantom.readout.lockin import LockinParams, lockin_sweep_a_per_v, lockin_sweep_dynamic_a_per_v
from phantom.readout.measure import SensedDoubleDot, raster_times_s
from phantom.readout.noise import (
    Composite,
    NoiseProcess,
    NoNoise,
    OneOverFNoise,
    TelegraphNoise,
    WhiteNoise,
)
from phantom.readout.sensor import SensorSETParams
from phantom.readout.transport import (
    ChannelTurnOn,
    channel_transmission_current_a,
    gate_leakage_current_a,
    ohmic_current_a,
    pinch_off_current_a,
)

MODEL_VERSION = 2
"""Layer model version: bumped whenever any output of the layer changes for the same inputs (ADR 0003)."""

__all__ = [
    "MODEL_VERSION",
    "ChannelTurnOn",
    "Composite",
    "LeadBarrier",
    "LeadRateLandscape",
    "LockinParams",
    "NoNoise",
    "NoiseProcess",
    "OneOverFNoise",
    "SensedDoubleDot",
    "SensorSETParams",
    "TelegraphNoise",
    "WhiteNoise",
    "channel_transmission_current_a",
    "gate_leakage_current_a",
    "latched_occupation",
    "lockin_sweep_a_per_v",
    "lockin_sweep_dynamic_a_per_v",
    "ohmic_current_a",
    "pinch_off_current_a",
    "raster_times_s",
]
