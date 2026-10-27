# The consumer contract

What phantom promises a consumer (ADR 0003, 0024). Everything here is reachable from one import,
`phantom.api`, is recorded in `docs/api/surface.json`, and is checked by `python -m phantom.contract`,
which a consumer runs against the release it pins and against its own device files.

```python
from phantom.api import CompileParams, Device, compile_device, optimum, two_dot_layout

op = {"B_L": 0.15, "P1": 0.6, "B_M": 0.2, "P2": 0.6, "B_R": 0.15, "SG": 0.6, "R_L": 0.7, "R_R": 0.7}
device = compile_device(two_dot_layout("SiGe", 80.0), CompileParams(operating_point_v=op), seed=2)
device.to_json("device.json")             # canonical JSON, about 8 KB; Device.from_json reloads it

run = device.session(noise_seed=7)        # one run: a clock, the device's noise, the dot charge
run.hold(v_pulse, 1e-6)                   # sit at v_pulse for 1 µs; the charge evolves there
m = run.sense(v_read, dwell_s=1e-3, samples=100, v_sd=50e-6)
m.times_s, m.current_a, m.n_present       # (100,), (100,), (100, 2)
```

## The surface

| Call | Returns | Notes |
|---|---|---|
| `compile_device(layout, CompileParams, seed)` | `Device` | deterministic in the seed; the slow tier (about 0.15 s) |
| `Device.to_json`, `dumps`, `from_json`, `loads` | | text-identical round trip; `DEVICE_FORMAT_VERSION` |
| `Device.at(v, v_sd=None)` | `State` | noise-free, no time; vectorised over leading axes of `v` (…, gates) |
| `Device.session(noise=True, dynamics=True, charge="sampled", noise_seed=None)` | `Session` | one run |
| `Session.sense(v, dwell_s, samples=None, v_sd=None)` | `Measurement` | reads; advances the clock |
| `Session.hold(v, duration_s)` | | no read; advances the clock; the charge evolves at `v` |
| `Session.wait(duration_s)` | | `hold` at the last point given |
| `Session.spin_readout(p_triplet, detuning_ev, shots, t_integration_s, t1_read_s)` | `SpinOutcome` | PSB readout of a prepared pair |
| `Session.time_s`, `Session.last_point_v` | | the clock (s); the last programmed point |
| `Device.measure(v, dwell_s, samples, v_sd, noise, dynamics, charge, noise_seed)` | `Measurement` | a fresh session's `sense` |
| `Device.knowable()` | `Knowable` | what a user may be shown |
| `Device.truth(lattice_gates, near_v, v_sd, v_programmed)` | `Truth` | hidden parameters |
| `optimum(device, t_target_ev)` | `Optimum` | the best achievable: sensor, lattice, coupling, readout, coherence |
| `model_versions()` | `dict[str, int]` | every module's `MODEL_VERSION`, for a consumer's provenance |

Voltages are always the **programmed** gate voltages, one per entry of `Device.gate_names`, in volts.
phantom maps them through the device's wiring draw (open and shorted gates) to the voltages the gates
carry; `State.v_carried_v` shows the result.

## Time

A session's clock starts at zero and only moves forward. `sense` holds each point for `samples` dwells
(one when `samples` is None) and reads it at the end of each: sample k of a call is read at
t₀ + (k + 1)·dwell, t₀ the clock when the call starts, and the clock ends at the last sample. `hold`
and `wait` advance it by their duration and read nothing; `spin_readout` advances it by
`shots × t_integration_s`. A zero or negative dwell, or a negative hold, raises `ValueError`. Before its
first call the device has sat at the first point long enough to be in equilibrium there.

## What one sample is

At each sample's time the noise processes of the device's draw are sampled: a common-mode gate offset
added to the carried voltages, the switchers near the sensor shifting its potential, and white readout
noise added to the current. The white floor keeps its density: σ is the draw's value at 1 ms times
√(1 ms / dwell). The charge in the sensor depends on two settings, fixed per session:

| | `dynamics=True` (default) | `dynamics=False` |
|---|---|---|
| `charge="sampled"` (default) | one state of a Markov jump trajectory of the rate equation: from the state present, the next is drawn from exp(W·dwell); consecutive samples are correlated when the rates are slow against the dwell | an independent Boltzmann draw at each sample (qdt-sim's behaviour) |
| `charge="mean"` | ⟨N⟩ of the distribution propagated by exp(W·dwell): latching and lag | the thermal ⟨N⟩; without noise this is `State.sensor_current_a` exactly |

`Measurement.n_present` is the integer state in `sampled` mode and ⟨N⟩ in `mean` mode. A hold
samples the gate noise once at its start and keeps it through the hold. With dynamics a sample costs
one matrix exponential, about 1 ms; without, about 0.05 ms.

## Seeds and reproducibility

- The device seed is the device's identity: it fixes the disorder, wiring, valley, spin and noise
  draws. A device file pins all of them even when a later phantom compiles differently.
- `noise_seed=None` replays the device's own noise realisation, so two fresh sessions given identical
  calls return identical results. An integer `noise_seed` keeps every noise parameter of the draw and
  gives a new realisation; the same integer replays it.
- Identical calls give identical output bit for bit on the same platform and numpy build, and agree
  within 1e-9 relative across platforms (ADR 0010, 0014). A point evaluated alone and the same point
  inside a batch agree within 1e-12 relative, not bit for bit.

## Knowable and truth

`Device.knowable()` holds what a person running the device knows before measuring: the material, the
gate names, the layout, the electron temperature, the programmed point, and the sensor, barrier and
reservoir gates. **Everything else on a `Device` is truth**, including the seed (recompiling from it
reproduces every draw), `params`, the scenario draw and `operating_point_v`, which holds the carried
voltages and the sensor set point compilation chose. A consumer that hides truth from a user forwards
`knowable()` and measurements, never the `Device`, `truth()` or `optimum()`. phantom has no hidden mode.

`Device.truth()` gathers the hidden parameters, and states two quantities a tune-up recovers:

- `charge_lattice_e_per_v`: electrons per programmed volt on (dot 1, dot 2) [rows] from the two
  `lattice_gates` [columns], through the wiring draw.
- `sensor` (`SensorTruth`, None when the sensor gate does not move the sensor): with the other gates at
  `v_programmed` (the programmed point by default) and the dots' charge held at its thermal value there,
  the peak in programmed sensor-gate volts nearest `near_v`, the period, the sensor's charging energy,
  its lever arm (charging energy / period), the width w of I ∝ sech²((V − V_peak)/w), the wing
  V_peak + w·atanh(1/√3) (the inflection) and the sensitivity |dI/dV| = 4h/(3√3 w) there at bias
  `v_sd`, h the linear-response peak height.

## Versions and releases

- `phantom.__version__` follows semantic versioning in 0.x: the minor version rises for any breaking
  change to this surface or to `DEVICE_FORMAT_VERSION`, the patch version otherwise. Every change to
  the surface shows in `docs/api/surface.json` and is listed in `CHANGELOG.md`.
- A change of output for the same parameters, seed and calls bumps the `MODEL_VERSION` of the module
  that changed. `Device.provenance` records the compile-time versions; `model_versions()` gives them all.
- The owner tags releases; a consumer pins a tag. How often phantom releases and consumers move their
  pins is an open decision (Step 8).

## Running the contract test

```
python -m phantom.contract                      # the reference device
python -m phantom.contract devices/*.json       # and every device file a consumer commits
```

It prints one line per check and exits 0 when every check holds. The checks: the surface imports;
compile is deterministic; the round trip is text-identical and the reload answers bit for bit;
provenance names every model version; `at` shapes, finiteness and bias; `knowable` fields; `truth`
fields and the sensor's relations; fresh sessions replay, noise seeds replay and differ; the clock;
equilibrium; bias in `sense`; sample semantics; a hold moves the charge; `optimum` has five finite
parts; `spin_readout` shapes and clock.

## Not yet in the contract

- Rates are compiled at the operating point and are constant over voltage at runtime.
- `spin_readout` takes a prepared triplet probability; spin driven by gate voltages over time is
  issue #41.
- The wiring draw's sensor-coupling asymmetry is recorded in the scenario but not applied to the sensor.
