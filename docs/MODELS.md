# How the models fit together (ADR 0011)

Two diagrams and a runtime ledger, kept current with every change to a model. The first diagram is
drawn by hand and says how a measurement is produced from a gate vector: which model feeds which, in the
order the physics happens. Whoever adds or rewires a model updates it in the same pull request. The
second is generated from the imports in `src/phantom` by `scripts/model_graph.py` and checked in CI, so
it can never be stale. The ledger records what each model and each protocol costs to run: declared
budgets from the protocol files, and the last measurement recorded by `scripts/gallery.py --record`.

## Measurement flow (hand-maintained)

Read left to right: a layout and a gate vector go in, a measured quantity comes out. Boxes are models;
the layer each lives in is named in the box. Dashed boxes are placeholders or offline checks that a later
step replaces (the step is named).

```mermaid
flowchart LR
  V["gate vector V(t), source–drain bias, time t"]

  subgraph layout
    LAY["geometry.Layout<br/>gate polygons, Stack (effective depth), MaterialClass; JSON file"]
  end

  subgraph disorder
    SCN["scenario.draw_scenario / apply_scenario<br/>material params + seed + severity → perturbed layout and basis; serialisable; severity 0 = clean"]
    CHG["charges.draw_fixed_charges, charge_energy_ev<br/>Poisson interface charges → image-charge energy offset (placeholder density)"]
    MIS["geometry_faults.draw_misalignment, apply_misalignment<br/>gate-layer shift and per-gate edge error as polygon transforms"]
    WIR["wiring.WiringDraw<br/>open and shorted gates (M·V + c), gate-to-gate leaks, sensor-coupling factors"]
    TRP["trapping.effective_gate_voltage<br/>history-dependent trapped fraction → screened gate voltage"]
  end

  subgraph electrostatics
    BASIS["closed_form.GateBasis (+ offset_ev)<br/>solid angle of each gate polygon / 2π on a grid: potential per volt"]
    DOTS["dots.find_dots, barrier_between<br/>minima → dot positions, harmonic radii, lever arms; saddles → barrier heights"]
    SCH["schrodinger.dot_levels_ev<br/>offline check of the harmonic estimate (compile tier)"]
  end

  subgraph charge
    DER["derived.derive_double_dot<br/>C_dot from radii, C_gate = C_dot·α → ConstantInteractionParams; sensor α from the basis"]
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

  LAY --> BASIS
  LAY --> SCN
  SCN --> CHG
  SCN --> MIS
  SCN --> WIR
  V --> WIR
  V --> TRP
  WIR --> BASIS
  TRP --> BASIS
  WIR --> SENS
  WIR --> TR
  MIS --> BASIS
  CHG --> BASIS
  BASIS --> DOTS
  DOTS --> DER
  DOTS -.-> SCH
  DER --> CI
  V --> BASIS
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
  style SCH stroke-dasharray: 4 3
  style CHG stroke-dasharray: 4 3
  style MIS stroke-dasharray: 4 3
  style WIR stroke-dasharray: 4 3
  style TRP stroke-dasharray: 4 3
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
| `potential-map-layout`, `potential-cut-along-channel` | Layout → GateBasis → U(x, y) at the operating point (Step 2) |
| `stability-coarse-2dot-layout` | Layout → GateBasis → dots → derived CI and sensor α → sensor flank set → sweep with noise (Step 2) |
| `pinch-off-crosstalk` | Layout → GateBasis at the barrier → linear pinch-off line; phenomenological step (Step 2) |
| `turn-on-global`, `turn-on-single-gate-layout`, `gate-hysteresis` | Layout → GateBasis → highest energy along the channel → transport.channel_transmission_current_a (Step 2) |
| `gate-leakage-iv`, `ohmic-iv` | Layout gate areas and stack thickness → transport.gate_leakage_current_a; transport.ohmic_current_a (Step 2, phenomenological) |
| `crosstalk-lines` | Layout → derived CI and sensor α → noise-free maps against four gate pairs (Step 2) |
| `potential-map-disorder` | Layout → scenario draw (charges, misalignment) → perturbed basis → U(x, y), clean beside draws (Step 3) |
| `stability-coarse-2dot-disorder` | Layout → scenario draw → perturbed basis → dots → derived CI and sensor α → sensor flank → sweep, per draw (Step 3) |
| `unintentional-dots` | Layout (MOS) → many scenario draws → minima count in the channel vs. density (Step 3) |
| `fault-signatures-stability` | Layout → wiring draw applied to programmed V → derived CI and sensor (β scaled) → sweep (Step 3) |
| `gate-leakage-pair` | wiring leak network + transport.gate_leakage_current_a (Step 3) |
| `gate-hysteresis-trapping` | trapping → effective gate voltage → channel barrier → transport.channel_transmission_current_a (Step 3) |

## Module dependency graph (generated; do not edit by hand)

<!-- module-graph:start -->
```mermaid
flowchart LR
  subgraph layout
    layout_geometry["geometry"]
  end
  subgraph electrostatics
    electrostatics_closed_form["closed_form"]
    electrostatics_dots["dots"]
    electrostatics_schrodinger["schrodinger"]
  end
  subgraph disorder
    disorder_charges["charges"]
    disorder_geometry_faults["geometry_faults"]
    disorder_scenario["scenario"]
    disorder_trapping["trapping"]
    disorder_wiring["wiring"]
  end
  subgraph charge
    charge_ci["ci"]
    charge_derived["derived"]
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
  charge_ci --> constants
  charge_derived --> charge_ci
  charge_derived --> charge_dot
  charge_derived --> electrostatics_closed_form
  charge_derived --> electrostatics_dots
  charge_derived --> layout_geometry
  charge_dot --> constants
  charge_hubbard --> charge_ci
  charge_hubbard --> charge_dot
  disorder_charges --> constants
  disorder_geometry_faults --> layout_geometry
  disorder_scenario --> disorder_charges
  disorder_scenario --> disorder_geometry_faults
  disorder_scenario --> disorder_wiring
  disorder_scenario --> electrostatics_closed_form
  disorder_scenario --> layout_geometry
  electrostatics_closed_form --> layout_geometry
  electrostatics_dots --> constants
  electrostatics_dots --> electrostatics_closed_form
  electrostatics_schrodinger --> electrostatics_dots
  readout_loading --> constants
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
| `charge.derived` | Charge-model parameters derived from a layout instead of typed (Step 2, ADR 0013). |
| `charge.dot` | Dots and barriers as the named parts of a charge model. |
| `charge.hubbard` | Two-site extended Hubbard model of a few-electron double dot. |
| `constants` | Physical constants (CODATA 2018), below every layer so any layer may import them (ADR 0004). |
| `disorder.charges` | Random fixed charges at an interface, as a perturbation of the electron energy landscape. |
| `disorder.geometry_faults` | Gate misalignment and width error as transforms of the layout polygons. |
| `disorder.scenario` | A scenario draw: material + catalog parameters + seed → a perturbed device, serialisable. |
| `disorder.trapping` | Charge trapping with gate history: the hysteresis loop in a gate sweep. |
| `disorder.wiring` | Gate wiring faults: open gates, shorted gates, gate-to-gate leakage, and sensor-coupling asymmetry. |
| `electrostatics.closed_form` | Closed-form potential of polygonal surface gates over a two-dimensional electron gas. |
| `electrostatics.dots` | Dots and barriers as features of the potential-energy landscape. |
| `electrostatics.schrodinger` | Two-dimensional effective-mass Schrödinger solve on a potential-energy window: dot levels and spacings. |
| `layout.geometry` | Gate polygons, the dielectric stack and the material class: data only, no physics. |
| `readout.loading` | Rate-limited loading of the dots from their reservoirs, and the lead-rate landscape over gate space. |
| `readout.lockin` | The `sensor_lockin` channel: a hardware derivative of the sensor signal. |
| `readout.measure` | A double dot read by an SET sensor, with noise on an explicit time axis. |
| `readout.noise` | Noise processes with an explicit time axis. |
| `readout.sensor` | A charge sensor: a single-electron transistor (SET) read in sequential tunnelling. |
| `readout.transport` | Channel turn-on and barrier pinch-off: the direct-current shapes of the electrical-health protocols. |
<!-- module-graph:end -->

## Runtimes (generated from `docs/gallery/measured.json` and the protocol files; do not edit by hand)

<!-- runtimes:start -->
Last recorded on 2026-10-26 at commit `eb0100e` on Windows-11-10.0.26200-SP0 (python 3.12.10, numpy 2.5.4), by `scripts/gallery.py --record`. Wall times are machine-dependent; counts are exact. The reference platform's numbers are in the CI artifact `gallery-ubuntu-latest` (ADR 0010).

### Per model (microseconds per point or sample)

| model | what is timed | points | us per point |
|---|---|---|---|
| `charge.ci.ConstantInteractionParams.mean_occupation` | thermal occupation per gate point | 10,000 | 1.43 |
| `charge.hubbard.HubbardParams.mean_occupation` | 16x16 diagonalisation per gate point | 10,000 | 9.872 |
| `readout.sensor.SensorSETParams.current_a` | sensor current per point | 100,000 | 0.118 |
| `readout.noise.TelegraphNoise.trace` | one fluctuator per sample | 100,000 | 0.658 |
| `readout.noise.OneOverFNoise.trace` | 24 fluctuators per sample | 10,000 | 16.021 |
| `readout.measure.SensedDoubleDot.sweep_a (noisy)` | occupation + sensor + three noise traces per point | 10,000 | 18.406 |
| `readout.lockin.lockin_sweep_a_per_v` | 8 evaluations per point, noisy | 2,500 | 145.759 |
| `readout.lockin.lockin_sweep_dynamic_a_per_v` | 8 evaluations + rate dynamics per point, noisy | 2,500 | 154.833 |
| `readout.loading.latched_occupation` | per scan step | 100,000 | 0.936 |

### Per protocol (declared budget beside the last measurement)

| protocol | points | compile budget s / measured | run budget ms per point / measured | evaluations per point |
|---|---|---|---|---|
| `coulomb-oscillations-single-dot` | 1,201 | 0.5 / 0.0 | 0.05 / 0.00191 | 1 |
| `crosstalk-lines` | 57,600 | 0.5 / 0.0939 | 0.05 / 0.00195 | 1 |
| `diamonds-single-dot` | 64,521 | 0.5 / 0.0 | 0.02 / 0.00014 | 1 |
| `fault-signatures-stability` | 57,600 | 2.0 / 0.0662 | 0.1 / 0.01905 | 1 |
| `gate-hysteresis-trapping` | 722 | 0.5 / 0.0424 | 1.0 / 0.07654 | 1 |
| `gate-hysteresis` | 361 | 0.5 / 0.0405 | 1.0 / 0.07208 | 1 |
| `gate-leakage-iv` | 1,806 | 0.5 / 0.0 | 0.02 / 9e-05 | 1 |
| `gate-leakage-pair` | 201 | 0.5 / 0.0 | 0.02 / 0.00027 | 1 |
| `noise-psd` | 30,000 | 0.5 / 0.0003 | 0.1 / 0.01878 | 1 |
| `noise-time-trace` | 30,000 | 0.5 / 0.0005 | 0.1 / 0.01861 | 1 |
| `ohmic-iv` | 201 | 0.5 / 0.0 | 0.02 / 6e-05 | 1 |
| `pinch-off-barrier` | 601 | 0.5 / 0.0 | 0.02 / 2e-05 | 1 |
| `pinch-off-crosstalk` | 24,321 | 0.5 / 0.0419 | 0.01 / 2e-05 | 1 |
| `plunger-vs-barrier-map` | 28,800 | 0.5 / 0.0001 | 0.05 / 0.00287 | 1 |
| `potential-cut-along-channel` | 1,203 | 0.5 / 0.0415 | 0.01 / 0.00035 | 1 |
| `potential-map-disorder` | 138,312 | 1.0 / 0.1577 | 0.001 / 0.0 | 1 |
| `potential-map-layout` | 34,578 | 0.5 / 0.0401 | 0.001 / 0.0 | 1 |
| `sensor-coulomb-peak` | 601 | 0.5 / 0.0 | 0.02 / 0.00187 | 1 |
| `stability-coarse-2dot-disorder` | 57,600 | 2.0 / 0.188 | 0.1 / 0.02261 | 1 |
| `stability-coarse-2dot-layout` | 22,500 | 0.5 / 0.0798 | 0.1 / 0.01913 | 1 |
| `stability-coarse-2dot` | 40,000 | 0.5 / 0.0004 | 0.1 / 0.0208 | 1 |
| `stability-fine-2dot` | 22,500 | 0.5 / 0.0005 | 0.1 / 0.0189 | 1 |
| `stability-lockin-common` | 40,000 | 0.5 / 0.0004 | 0.5 / 0.16262 | 8 |
| `stability-lockin-detuning` | 40,000 | 0.5 / 0.0004 | 0.5 / 0.16237 | 8 |
| `stability-lockin-p1` | 40,000 | 0.5 / 0.0006 | 0.5 / 0.16375 | 8 |
| `stability-lockin-p2` | 40,000 | 0.5 / 0.0004 | 0.5 / 0.1631 | 8 |
| `turn-on-global` | 201 | 0.5 / 0.0363 | 1.0 / 0.07725 | 1 |
| `turn-on-single-gate-layout` | 361 | 0.5 / 0.0386 | 1.0 / 0.07483 | 1 |
| `turn-on-single-gate` | 601 | 0.5 / 0.0 | 0.02 / 7e-05 | 1 |
| `unintentional-dots` | 125 | 1.0 / 0.0 | 400.0 / 82.70662 | 1 |
<!-- runtimes:end -->
