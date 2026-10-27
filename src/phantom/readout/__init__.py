"""phantom.readout — Lead tunnel rates, master equation, latching, charge-sensor (SET) and RF response, noise processes (white, telegraph, 1/f, fluctuator ensembles) with explicit time steps. Seeded from qdt-sim's noise models (Step 1).

Modules: `sensor` (SET sensor and single-dot transport current), `noise` (seeded processes with a time
axis), `measure` (a CI double dot read by the sensor with noise), `lockin` (the `sensor_lockin`
hardware-derivative channel), `transport` (phenomenological turn-on and pinch-off shapes, Step 1 only), `loading` (Step 1 placeholder for
rate-limited loading; superseded for dynamics by `master_equation`), `rates` (WKB tunnel rates from barrier height
and width), `master_equation` (rate equation over the CI charge states with detailed balance; latching from the lag), `rf`
(lumped tank circuit with the sensor as load: reflection and resonance), `single_shot` (integrated shots,
histograms, signal-to-noise against integration time).
"""

from phantom.readout.loading import LeadRateLandscape, latched_occupation
from phantom.readout.lockin import LockinParams, lockin_sweep_a_per_v, lockin_sweep_dynamic_a_per_v
from phantom.readout.master_equation import RateEquation
from phantom.readout.measure import SensedDoubleDot, raster_times_s
from phantom.readout.noise import (
    Composite,
    NoiseProcess,
    NoNoise,
    OneOverFNoise,
    TelegraphNoise,
    WhiteNoise,
)
from phantom.readout.rates import (
    attempt_rate_hz,
    wkb_rate_from_action_hz,
    wkb_transmission,
    wkb_tunnel_rate_hz,
)
from phantom.readout.rf import TankCircuit, sensor_resistance_ohm
from phantom.readout.sensor import SensorSETParams
from phantom.readout.single_shot import ShotResult, noisy_device_for_sampling, single_shots
from phantom.readout.transport import (
    ChannelTurnOn,
    channel_transmission_current_a,
    gate_leakage_current_a,
    ohmic_current_a,
    pinch_off_current_a,
)

MODEL_VERSION = 4
"""Layer model version: bumped whenever any output of the layer changes for the same inputs (ADR 0003)."""

__all__ = [
    "MODEL_VERSION",
    "ChannelTurnOn",
    "Composite",
    "LeadRateLandscape",
    "LockinParams",
    "NoNoise",
    "NoiseProcess",
    "OneOverFNoise",
    "RateEquation",
    "SensedDoubleDot",
    "SensorSETParams",
    "ShotResult",
    "TankCircuit",
    "TelegraphNoise",
    "WhiteNoise",
    "attempt_rate_hz",
    "channel_transmission_current_a",
    "gate_leakage_current_a",
    "latched_occupation",
    "lockin_sweep_a_per_v",
    "lockin_sweep_dynamic_a_per_v",
    "noisy_device_for_sampling",
    "ohmic_current_a",
    "pinch_off_current_a",
    "raster_times_s",
    "sensor_resistance_ohm",
    "single_shots",
    "wkb_rate_from_action_hz",
    "wkb_transmission",
    "wkb_tunnel_rate_hz",
]
