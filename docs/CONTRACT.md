# The consumer contract

What phantom promises a consumer (ADR 0003, 0024). Everything here is reachable from one import,
`phantom.api`, is recorded in `docs/api/surface.json`, and is checked by `python -m phantom.contract`,
which a consumer runs against the release it pins and against its own device files.

```python
from phantom.api import CompileParams, Device, compile_device, optimum, two_dot_layout

op = {"B_L": 0.15, "P1": 0.6, "B_M": 0.2, "P2": 0.6, "B_R": 0.15, "SG": 0.6, "R_L": 0.7, "R_R": 0.7,
      "SB_L": 0.4, "SB_R": 0.3, "SR_L": 0.7, "SR_R": 0.7}
device = compile_device(two_dot_layout("SiGe", 80.0), CompileParams(operating_point_v=op), seed=2)
device.to_json("device.json")             # canonical JSON, about 13 KB; Device.from_json reloads it

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

`Measurement.n_present` is the charge that entered the sensor. In `sampled` mode it is one eigenstate of
the dots' charge (ADR 0029): its total charge is an integer, and its split between the dots is the
eigenstate's ⟨N⟩, an integer up to the interdot admixture (t/Δ)² away from interdot lines and near half an
electron each way at the middle of one when t > kT. In `mean` mode it is ⟨N⟩. A hold
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
gate names, the layout, the electron temperature, the programmed point, the sensor, barrier and
reservoir gates, and the sensor's own barrier and bath gates (`sensor_barrier_gates`,
`sensor_reservoir_gates`, empty for a layout without them). **Everything else on a `Device` is truth**, including the seed (recompiling from it
reproduces every draw), `params`, the scenario draw and `operating_point_v`, which holds the carried
voltages and the sensor set point compilation chose. A consumer that hides truth from a user forwards
`knowable()` and measurements, never the `Device`, `truth()` or `optimum()`. phantom has no hidden mode.

`Device.truth()` gathers the hidden parameters, and states two quantities a tune-up recovers:

- `charge_lattice_e_per_v`: electrons per programmed volt on (dot 1, dot 2) [rows] from the two
  `lattice_gates` [columns], through the wiring draw.
- `sensor` (`SensorTruth`, None when the sensor gate does not move the sensor): with the other gates at
  `v_programmed` (the programmed point by default) and the dots' charge held at its thermal value there,
  the peak in programmed sensor-gate volts nearest `near_v`, the period, the sensor's charging energy,
  its lever arm (charging energy / period), the thermal width w of I ∝ sech²((V − V_peak)/w), the
  sensor's two tunnel rates at that peak and their lifetime width ħ(Γ_S + Γ_D) in sensor-gate volts,
  the wing (the steepest point of the linear-response peak above it) and the sensitivity |dI/dV| there
  at bias `v_sd`. The peak is the thermal sech² convolved with a Lorentzian of the lifetime width
  (ADR 0026); with a negligible lifetime width the wing is V_peak + w·atanh(1/√3) and the sensitivity
  4h/(3√3 w), h the linear-response peak height. The rates are those at the peak: along the sensor gate
  they change by a factor e every few mV, so the truth states each peak with its own rates.

The sensor of a compiled device is a dot with two barriers and a bath behind each (ADR 0026). Its peaks
are a few hundred microvolts wide and tens of millivolts apart on the sensor gate, so the programmed
sensor gate alone usually sits between peaks: tune it first, as on a device (the compiled set point in
`operating_point_v`, the truth's peak and wing, and `optimum(device).sensor` are where a tune-up lands).
`State.sensor_gamma_hz` gives its source and drain rates at any point: they follow the sensor's barrier,
bath and plunger gates and every other gate's cross-talk, and a bath gate below its turn-on empties that
side.

The dots' own rates follow the gates too (ADR 0027): `State.gamma_left_hz`, `gamma_right_hz` and
`gamma_interdot_hz` are arrays with the input's leading axes, the compiled rates moved by their slopes
against every gate and capped at the attempt rate. A session's dynamics use the rates at each sample's
voltages: closing a barrier gate slows its lead until a sweep latches and a parked point telegraphs.
`Truth` states the rates at the operating point.

The interdot tunnel coupling follows every gate as well (ADR 0029, 0034): `State.tunnel_coupling_ev` is
(ħω₀/π)·√T, with T the interdot barrier's one transmission (below), and the charge states of one total
charge hybridize through it. An interdot line is then the two-level polarization ½(1 − ε/Ω·tanh(Ω/2kT)),
Ω = √(ε² + 4t²), so opening the barrier gate widens it; `State.n_mean` and the sensor read the
eigenstates, and a session's dynamics run over them. A dot is isolated only when its own barrier and the
interdot barrier are both closed: with B_M open, dot 2's lead loads dot 1 in second order in t.
`optimum(device).coupling` inverts this t along the barrier gate.

The barrier gates also move and reshape the dots (ADR 0030). Compile re-derives the double dot along
each barrier gate, and at runtime the capacitances, lever arms and offsets at the dots, t and the rates
follow it. `State.dot_positions_nm`, `dot_radii_nm` and `addition_energies_ev` report the shape, and
`State.dots_merged` flags a point where a barrier gate is opened past the last point at which the
potential forms two dots; the shape there is held at that point, and the charge model is not meant to
describe it, but the interdot barrier's top keeps following the gates there (ADR 0034), so a far-open
B_M leaves one dot at every occupancy. Because the dots' charge is counted absolutely, a barrier gate
moves the plunger lines by more than its cross-capacitance.

The exchange follows them as well (ADR 0031). `State.exchange_ev` is the singlet–triplet splitting of the
last electron on each dot in the ground charge state, at the gated t and the shaped dots, and
`State.exchange_frequency_hz` is J/h. It is a qubit pair's exchange where both dots hold an odd number of
electrons (read the parity from `n_ground`), and it jumps where the ground charge state changes.
`Session.spin_readout` reads at the interdot coupling of the last point the session was given.

The interdot barrier has one transmission (ADR 0034): Kemble's, of its top above the interdot tunnelling
level, `State.interdot_barrier_top_ev`. The top is moved by every gate through two placeholder screening
factors (`CompileParams.interdot_barrier_screening` for the barrier gate, times `interdot_cross_screening`
for every other gate; 1 and 1 is the bare potential), and along the barrier gates through the shape
grid, continued past the merge. t, `gamma_interdot_hz` = f_a·T and the conductance all follow it: on the
reference device t changes by a factor e every 23 mV of B_M, and a plunger acts about four times more
weakly than B_M.

The dots grow with their occupancy (ADR 0035). Each dot is a classical puddle whose radius scales as
N^(1/3), anchored at the compiled dot: `State.puddle_radii_nm`. Its capacitances follow its gate-induced
charge, so `State.addition_energies_ev` falls as N^(−1/3) toward high occupancy and the honeycomb's
period shrinks with it (more electrons per volt), and the interdot barrier is seen from the puddles'
edges: at fixed gates it is tens of meV higher at the first electrons than at the compiled occupancy,
so t, the interdot rate and the conductance fall toward (1,1). At the operating point nothing changes.
`CompileParams.puddles=False` keeps the dots at their compiled size.

The interdot conductance sets the honeycomb's triple-point splitting (ADR 0032).
`State.interdot_conductance` is that transmission per channel, and `State.interdot_splitting` is the
fractional splitting f, the interdot segment's length as a fraction of the cell along anti-detuning.
Opening the barrier gate between the dots takes f from its geometric value (small hexagons) towards one
(the diagonal lines of a superdot) over about 65 mV of B_M on the reference device; closing it shrinks
the segments back.

With dynamics, the sensor dot's own charge is a process on the same clock (ADR 0028). Each sample's current
is the charge that left to the drain during the dwell over the dwell, so a fast sensor reads its
stationary current with shot noise shrinking as 1/√dwell, and a sensor slowed by its barriers or baths
reads single-electron bursts while its charge switches. `Measurement.sensor_charge` holds that charge at
the end of each dwell (only its changes are physical). The sensor acts back: its charge relative to its
stationary mean shifts each dot's levels by −β_i per electron. Without dynamics nothing of this applies,
and a long `mean`-mode dwell still reads `State.sensor_current_a` (within 1e-6).

## Versions and releases

- `phantom.__version__` follows semantic versioning in 0.x: the minor version rises for any breaking
  change to this surface or to `DEVICE_FORMAT_VERSION`, the patch version otherwise. Every change to
  the surface shows in `docs/api/surface.json` and is listed in `CHANGELOG.md`.
- A change of output for the same parameters, seed and calls bumps the `MODEL_VERSION` of the module
  that changed. `Device.provenance` records the compile-time versions; `model_versions()` gives them all.
- A release is an annotated tag `vX.Y.Z` on a green `main` commit, made by the owner or by a session the owner
  tells to (ADR 0025); a consumer pins a tag. How often phantom releases and consumers move their
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

- Every rate, the dots' and the sensor's, follows the gate voltages through a compiled exponential that
  holds within about 20 % of the compile tier's change over ±10 mV and drifts beyond (ADR 0026, 0027).
- Within one dwell the sensor sees the dots' charge at its end and the dots see the sensor's at its
  start (a first-order splitting, ADR 0028).
- `spin_readout` takes a prepared triplet probability; spin driven by gate voltages over time is
  issue #41.
