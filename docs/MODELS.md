# How the models fit together (ADR 0011)

Two diagrams, kept current with every change to a model. The first is drawn by hand and says how a
measurement is produced from a gate vector: which model feeds which, in the order the physics happens.
Whoever adds or rewires a model updates it in the same pull request. The second is generated from the
imports in `src/phantom` by `scripts/model_graph.py` and checked in CI, so it can never be stale.

## Measurement flow (hand-maintained)

Read left to right: a gate vector and a time go in, a measured quantity comes out. Boxes are models;
the layer each lives in is named in the box. Dashed boxes are Step 1 placeholders that a later step
replaces (the step is named).

```mermaid
flowchart LR
  V["gate vector V(t), source–drain bias, time t"]

  subgraph charge
    CI["ci.ConstantInteractionParams<br/>occupation N(V) at width kT; lever arms, triple points"]
    HB["hubbard.HubbardParams<br/>N(V) with tunnel coupling; t → 0 gives CI"]
    BAR["dot.Barrier t(V_B), dot.Dot E_c(radius)"]
  end

  subgraph readout
    RATES["loading.LeadRateLandscape / LeadBarrier<br/>dot-to-lead rate Γ(V); lifetime width ħΓ<br/>(placeholder → Step 4 master equation)"]
    LAT["loading.latched_occupation<br/>charge actually present given Γ and dwell time"]
    SENS["sensor.SensorSETParams<br/>φ_s = αV + βN + δφ → current I(φ_s, V_sd)"]
    NOISE["noise.WhiteNoise / TelegraphNoise / OneOverFNoise<br/>gate offset g(t), sensor offset s(t), readout r(t)"]
    MEAS["measure.SensedDoubleDot<br/>sweep_a, trace_a: raster times, noise applied"]
    LCK["lockin.lockin_sweep_a_per_v (quasi-static)<br/>lockin_sweep_dynamic_a_per_v (with Γ)"]
    TR["transport.ChannelTurnOn, pinch_off_current_a<br/>(placeholder → Step 2 potential)"]
  end

  OUT["arrays: current (A), demodulated dI/dV (A/V), N present"]

  V --> CI
  V --> HB
  BAR --> HB
  V --> RATES
  CI --> LAT
  RATES --> LAT
  LAT --> SENS
  CI --> SENS
  HB --> SENS
  NOISE --> MEAS
  SENS --> MEAS
  MEAS --> LCK
  RATES --> LCK
  V --> TR
  MEAS --> OUT
  LCK --> OUT
  TR --> OUT

  style RATES stroke-dasharray: 4 3
  style TR stroke-dasharray: 4 3
```

Which protocol uses which path (Step 1):

| protocol | path |
|---|---|
| `turn-on-single-gate`, `pinch-off-barrier` | V → transport → current |
| `coulomb-oscillations-single-dot`, `diamonds-single-dot`, `sensor-coulomb-peak` | V → CI (plungers at 0) → sensor → current; the sensor is itself the single dot in transport |
| `stability-coarse-2dot`, `stability-fine-2dot` | V(t) → CI → sensor, with noise through measure |
| `stability-lockin-p1`, `-p2`, `-common`, `-detuning` | V(t) chopped → rates → CI (broadened) → latched occupation → sensor → demodulation |
| `plunger-vs-barrier-map` | V → CI with a lead-barrier gate → LeadBarrier rate → latched occupation → sensor |
| `noise-time-trace`, `noise-psd` | fixed V, t → noise → measure → current |

## Module dependency graph (generated; do not edit by hand)

<!-- module-graph:start -->
```mermaid
flowchart LR
  subgraph charge
    charge_ci["ci"]
    charge_dot["dot"]
    charge_hubbard["hubbard"]
  end
  subgraph readout
    readout_loading["loading"]
    readout_lockin["lockin"]
    readout_measure["measure"]
    readout_noise["noise"]
    readout_sensor["sensor"]
    readout_transport["transport"]
  end
  charge_dot --> charge_ci
  charge_hubbard --> charge_ci
  charge_hubbard --> charge_dot
  readout_loading --> charge_ci
  readout_lockin --> readout_loading
  readout_lockin --> readout_measure
  readout_measure --> charge_ci
  readout_measure --> readout_noise
  readout_measure --> readout_sensor
  readout_transport --> charge_dot
```

| module | what it owns |
|---|---|
| `charge.ci` | Constant-interaction (CI) model of a few-electron double quantum dot. |
| `charge.dot` | Dots and barriers as the named parts of a charge model. |
| `charge.hubbard` | Two-site extended Hubbard model of a few-electron double dot. |
| `readout.loading` | Rate-limited loading of the dots from their reservoirs, and the lead-rate landscape over gate space. |
| `readout.lockin` | The `sensor_lockin` channel: a hardware derivative of the sensor signal. |
| `readout.measure` | A double dot read by an SET sensor, with noise on an explicit time axis. |
| `readout.noise` | Noise processes with an explicit time axis. |
| `readout.sensor` | A charge sensor: a single-electron transistor (SET) read in sequential tunnelling. |
| `readout.transport` | Channel turn-on and barrier pinch-off: the direct-current shapes of the electrical-health protocols. |
<!-- module-graph:end -->
