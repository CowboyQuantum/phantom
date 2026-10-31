# Changelog

phantom follows semantic versioning in 0.x (ADR 0024): the minor version rises for any breaking change
to the consumer surface (`phantom.api`, recorded in `docs/api/surface.json`) or to the device file
format, the patch version otherwise. A change of output for the same parameters, seed and calls bumps
the `MODEL_VERSION` of the module that changed and is listed here. A release is tagged `vX.Y.Z` on a green
`main` commit, by the owner or by a session the owner tells to (ADR 0025).

## 0.4.0 — the barrier gates set the interdot coupling, the dots' shape and the exchange (#58, #59, #60; untagged)

The interdot tunnel coupling t follows every gate at runtime, and the charge states of one total charge
hybridize through it (ADR 0029): opening B_M widens an interdot line from the thermal step to the
two-level polarization, and closing it is what isolates a dot. The barrier gates move and reshape the
dots through a grid compiled along each of them (ADR 0030): positions, sizes, lever arms, capacitances,
t and the rates follow, until the dots merge. The exchange of the dots' valence electrons follows them
(ADR 0031): J comes from the CI window's local sector at the gated t and the shaped dots.

### Breaking
- `State.tunnel_coupling_ev` follows every gate, not the barrier gate alone: the compiled t moved by half
  the slopes of ln Γ_m (Γ_m ∝ t²) and capped at ħω₀/π, compile's t at a transmission of one; along the
  barrier gates, half the change of the interdot action on the shape grid. On the reference device it
  runs from below 0.1 µeV to above 100 µeV as B_M moves by tens of mV around the programmed point.
- `State.n_mean`, `State.n_ground`, `State.sensor_current_a` and every session read the hybridized
  eigenstates of the CI window: near an interdot line the charge is no longer a thermal step, and away
  from one it carries an admixture (t/Δ)² of the neighbouring charge state.
- `Measurement.n_present` in `sampled` mode is the charge of one eigenstate: its total is an integer, its
  split between the dots is not. The contract checks the total (`sample_semantics`); a consumer that
  labels charge states rounds it.
- With dynamics, dot 1 loads through dot 2's lead in second order in t, Γ_R t²(1/Δ + 1/Δ')²: closing B_L
  alone no longer latches dot 1 on the reference device (about 10⁶ Hz through dot 2 against 10⁴ Hz
  through B_L). Closing B_L and B_M together does. The contract's hold check compares the total charge.
- `optimum(device).coupling.v_barrier_v` inverts the runtime coupling along the barrier gate (Brent's method on ln t)
  (NaN outside the range t takes within ±1 V).
- `DEVICE_FORMAT_VERSION` 4: `Device.shape`, the double dot's changes on a grid along B_L, B_M and B_R.
  Files from 0.3.0 do not load; recompile them from their seed.
- The barrier gates move the dots (#59): the runtime CI model's capacitances and induced charge, t and
  the lead and interdot rates follow them, so every map with a barrier gate off its operating value
  changes. A growing dot takes electrons at a fixed potential (the charge is counted absolutely), so a
  barrier gate moves the plunger lines by more than its cross-capacitance.
- Along a barrier gate a rate levels off at f_a/(2 − T₀) once the barrier falls below its tunnelling
  level (its action vanishes), not at the attempt rate; the barrier gates' changes add.

- `State.exchange_ev` is J of the dots' valence electrons (#60): the singlet secular equation on the CI
  window's local sector, (N₁, N₂) with (N₁±1, N₂∓1), at the gated t and the shaped dots. It is a qubit
  pair's exchange where both dots hold an odd number of electrons, and it jumps where the ground charge
  state changes. Before, J came from a 0–2 electron Hubbard model and read about 13 meV at the reference
  device's programmed point, whatever t was.
- `Session.spin_readout` reads at the interdot coupling of the last point the session was given, not at
  the compiled t.

### Added
- `charge.hybrid`: `InterdotCoupling` (t at the gates) and the CI window's eigenstates, batched per total
  charge (`window_eigenstates`, `eigenbasis`); `ConstantInteractionParams.coupling`, `tunnel_coupling_ev`,
  `eigenstates`; `phantom.charge.InterdotCoupling`.
- `Device.runtime_coupling()`; `RateEquation.eigen_rate_matrix`, and `readout.master_equation.project`
  and `carried_states` for a distribution over eigenstates.
- Gallery protocol `interdot-line-vs-barrier` (unsigned): t and Γ_m along B_M, detuning cuts at three B_M
  beside the two-level formula, and the interdot line over B_M.
- `charge.shape.ShapeGrid`; `ConstantInteractionParams.shape`, `local`, `addition_energies_at_ev`;
  `InterdotCoupling.shape`, `GateRates.shape`; `charge.derived.capacitance_matrix_af`;
  `electrostatics.dots.refined_dot`; `CompileParams.shape_gates`, `shape_span_v`, `shape_points`.
- `State.dot_positions_nm`, `dot_radii_nm`, `addition_energies_ev` and `dots_merged`.
- `State.exchange_frequency_hz` (J/h); `spin.exchange.local_exchange_ev` and `singlet_exchange_ev` (the
  secular solver, now returning J = −min(0, Δ_R, Δ_L) where t = 0 without dividing zero by zero).
- Gallery protocol `exchange-frequency-vs-barrier` (unsigned): t and J/h along B_M, and J/h across the
  nearest odd–odd cell at three B_M beside its second-order form.
- Gallery protocol `dot-shape-vs-barrier` (unsigned): positions, radii, addition energies, C_m, t and the
  merged flag along each barrier gate, and charge maps over P1 and P2 at three B_M.

### Changed
- `session-closing-barrier` sweeps three ways: barriers open, B_L closed alone (loaded through dot 2),
  and B_L and B_M closed (latching, telegraph); its budget rises from 1,665 to 2,067 points.
- The hybridized window is ±4 in the statics and follows t in the dynamics, so the window's truncation no
  longer leaves a seam at the middle of an interdot segment, where rint(n_g) changes (ADR 0030).
- `Session` propagates the dots' charge with `readout.master_equation.propagator`, exp(W·dt) accurate in
  every entry relative to itself (non-negative Taylor steps and squarings), not with `expm`: a slow
  transition beside rates 10⁵ times faster keeps its digits on every platform (ADR 0030).
- Model versions: `charge.hybrid` 2 (new), `charge.shape` 1 (new), `readout.master_equation` 3,
  `compile.device` 5, `device` 6, `device.session` 7, `optimum.targets` 3, `spin.exchange` 3. Without a coupling and a
  shape, `ConstantInteractionParams`, `RateEquation` and `GateRates` give their 0.3.0 output bit for bit.
- The device file grows from about 14 kB to 62 kB (the shape grid); compile takes about 0.2 s more.

## 0.3.0 — the sensor dot's barriers and bath, and gate-dependent charge dynamics (#47, #53, #54; `v0.3.0` at `6920f19`)

The charge sensor becomes a dot like the qubit dots, with its own tunnel barriers and a bath behind each
that can be tuned (ADR 0026). Every compiled device's sensor changes. The qubit dots' lead and interdot
rates follow the gates at runtime (ADR 0027), so a closed barrier latches and telegraphs in a session,
and with dynamics the sensor dot's own charge evolves in time and acts back on the dots (ADR 0028).

### Breaking
- `State.gamma_left_hz`, `gamma_right_hz` and `gamma_interdot_hz` are arrays (...) at the carried
  voltages, no longer the compiled constants (#53). A session with dynamics uses the rates at each
  sample's voltages, so its output changes wherever the gates leave the operating point.
- A session with dynamics and a sensor with leads reads the sensor's counted charge per dwell: shot
  noise, single-electron bursts when the sensor is slow, and the back-action of the sensor's charge on
  the dots (#54). Its currents and charges change; without dynamics nothing changes.
- `DEVICE_FORMAT_VERSION` 3: each `CompiledBarrier` carries `attempt_rate_hz`, `tunnelling_level_ev` and
  `log_rate_slopes_per_v`. Files compiled on earlier 0.3.0 commits (format 2) do not load; recompile
  them from their seed.
- `two_dot_layout` gains four gates after the first eight: sensor barriers `SB_L`, `SB_R` beside `SG` and
  sensor bath (reservoir) gates `SR_L`, `SR_R` beyond them; gate vectors are 12 long and every operating
  point needs the four. The sensor's search region spans its barriers. The gallery's and the contract's
  programmed point adds `SB_L: 0.4, SB_R: 0.3, SR_L: 0.7, SR_R: 0.7`, set for seed 2; at that point seeds
  0, 3, 6 and 9 of 0–11 form no device (6 as before).
- `DEVICE_FORMAT_VERSION` 2: `Device` gains `sensor_charging_energy_ev`, `sensor_beta_v_per_e` and
  `sensor_leads`. Device files from 0.2.0 do not load; recompile them from their seed.
- A compiled sensor with leads reads φ_s as the potential at its dot in volts (lever arm 1 eV/V): its
  peaks are its charging energy apart (about 35 mV on SG for the gallery device, against about 1 V
  before), its thermal width is the device's 2kT, and its per-electron coupling comes from the layout
  (about −0.35 mV, against the typed −4 and −3 mV). Its current comes through its leads: no typed
  `peak_a` or `background_a`. The programmed sensor gate alone usually sits between peaks; tune it (the
  compiled set point and `optimum(...).sensor` do).
- The noise draw's switcher amplitudes are converted to the sensor dot's units (times its mean coupling
  over the typed sensor's), so a switcher keeps its size against a dot electron's signal; the gate 1/f
  noise now moves the sharper sensor visibly.
- `SensorTruth.wing_v` and `sensitivity_a_per_v` describe the lifetime-broadened peak (its steepest
  point and the slope there) and equal the sech² statements only when the lifetime width is negligible.

### Added
- `readout.broadening`: the Fermi function convolved with a Lorentzian level through the complex
  digamma function, the linear-response lineshape, its wing by bisection, the peak conductance, and
  `polygamma` for complex arguments.
- `readout.sensor.SensorLeads`: the sensor's rates compiled from the WKB action through its barriers at a
  fixed tunnelling level, their slopes against every gate, the attempt-rate ceiling, and each bath's
  occupation at the accumulation threshold; `SensorSETParams.leads`, `rates_hz`, `lifetime_width_v`,
  `linear_response_peak_a`, `wing_offset_v`, `wing_slope_a_per_v`; `current_a(phi, v_sd, v)` takes the
  carried gate voltages. `phantom.api.SensorLeads`.
- `charge.derived.sensor_coupling`: the sensor's charging energy and coupling from the three-dot
  capacitance matrix; `DerivedDoubleDot.sensor_charging_energy_ev`, `sensor_beta_v_per_e`.
- `CompileParams.sensor_barrier_gates`, `sensor_reservoir_gates`, `bath_threshold_ev`, `bath_width_ev`.
- `State.sensor_gamma_hz`; `Knowable.sensor_barrier_gates`, `sensor_reservoir_gates`;
  `SensorTruth.lifetime_width_v`, `gamma_source_hz`, `gamma_drain_hz`; `Device.switcher_scale`.
- `readout.rates.GateRates` and `saturating_rates_hz` (the sensor's leads use it, bit for bit as before);
  `RateEquation.gate_rates` and `rates_hz(v)`; `Device.runtime_rates()` gives the gated rate equation.
- `readout.sensor_dynamics` (the sensor's charge chain, counting cumulants, exact trajectories);
  `Measurement.sensor_charge`; `RateEquation.rate_matrix(..., level_shift_ev)`.
- Gallery protocol `session-slow-sensor` (unsigned): the sensor's rates as its barriers close, the
  parked current with the barriers open and closed, the sensor's charge, and a dot transition read by
  a fast and by a slow sensor.
- Gallery protocol `session-closing-barrier` (unsigned): the rates along B_L, a P1 sweep up and down
  with B_L open and closed, and a telegraph trace at the closed transition.
- `electrostatics.tunnelling_level_ev` (the level rule `barrier_action` already used, now shared) and
  `GateBasis.offset_at`.
- Gallery protocols `sensor-bath-turn-on`, `sensor-barrier-barrier-map`, `sensor-barrier-rates`,
  `sensor-peak-broadening` (unsigned); `compiled-device-draws` records the sensor rates; the session
  protocols tune the sensor gate first; `optimum-sensor-flank` sweeps ±30 mV.
- Reference entries `beenakker1991-coulomb-blockade`, `foxman1993-quantum-levels-coulomb-island`
  (`proposed`).

### Changed
- Compile's sensor set point: one step and a re-derivation as before, then a last step with the
  re-derived sensor and its leads; the double dot is derived again at the final point.
- `session-hold-relaxation` overlays exp(−Γ_out t), Γ_out the total rate out of the held charge state at
  the read point: the interdot path through dot 2's lead doubled the relaxation rate on the gallery
  device, also before #53.
- Model versions: `layout.geometry` 2, `charge.derived` 2, `compile.device` 4, `device` 4,
  `device.session` 4, `readout.sensor_dynamics` 1 (new), `device.truth` 2, `readout.broadening` 1 (new); layer versions `layout` 2,
  `charge` 4. The typed sensor (no leads) gives the same output as in 0.2.0.
- `gate-leakage-iv` evaluates twelve gates: its budget rises from 2,408 to 3,612 points.

## 0.2.0 — the consumer contract (Step 8; `v0.2.0` at `f72daac`)

First release a consumer can pin. Only the fix below changes the output of a 0.1.0 model.

### Fixed
- The compiled device applies the wiring draw's sensor-coupling asymmetry (ADR 0016, catalog #11; #43):
  `device.compiled_sensor` scales each dot's per-electron sensor coupling by the draw's factor, and both
  `Device.runtime_sensor` and compile's sensor set point use it. Every compiled device with a wiring draw
  reads a different sensor current and sits at a different compiled sensor set point; device files
  compiled before this change keep their old set point and are off the flank of the new sensor.
  `device` MODEL_VERSION 2, `compile` MODEL_VERSION 2; `Device.provenance` records the `device` version.

### Added
- `phantom.api`: the one import surface for consumers, and `model_versions()`, the `MODEL_VERSION` of
  every module (ADR 0024, `docs/CONTRACT.md`).
- `Device.at(v, v_sd)` → `State` (`device.session.state_at`): carried voltages, thermal and ground-state
  charge, sensor potential and current at a source–drain bias, interdot t, exchange J, compiled rates;
  noise-free and vectorised.
- `Device.session(noise, dynamics, charge, noise_seed)` → `Session` with a forward-only clock:
  `sense(v, dwell_s, samples, v_sd)` → `Measurement(times_s, current_a, n_present)`, `hold(v, duration_s)`,
  `wait(duration_s)`, `spin_readout(p_triplet, detuning_ev, shots, t_integration_s, t1_read_s)` →
  `SpinOutcome`; `time_s`, `last_point_v`. `Device.measure` is a fresh session's `sense`.
  - the sensor's source–drain bias reaches `sense`, per point or as a scalar;
  - `noise_seed` gives a new noise realisation with the device's noise parameters; `None` replays the
    device's own;
  - `hold` keeps the gates at a point and evolves the charge there without reading;
  - one sample's charge is stated: a Markov jump trajectory of the rate equation or an independent
    Boltzmann draw (`charge="sampled"`, default), or ⟨N⟩ (`charge="mean"`).
- `Device.truth(lattice_gates, near_v, v_sd, v_programmed)` → `Truth` (`device.truth`): the hidden
  parameters, the charge lattice in electrons per programmed volt through the wiring draw, and
  `SensorTruth` (peak, period, charging energy, lever arm, width, wing, sensitivity at a bias).
- `Device.knowable()` → `Knowable`: the fields a consumer may show a user; everything else is truth.
- `phantom.contract`: `run(device)` and `python -m phantom.contract [device.json ...]`, exit status 0 or 1.
- `scripts/api_surface.py` and `docs/api/surface.json`, checked by `scripts/check_all.sh`.
- Gallery protocols `session-charge-trace` and `session-hold-relaxation` (unsigned).
- Model versions: `device.session` 1, `device.truth` 1; `Device.sensor_beta_factors()`.

### Changed
- Rule 2: `phantom.api` and `phantom.contract` sit above `compile` and may import it; no layer or root
  module imports them (`scripts/check_imports.py` enforces it, ADR 0024).
- `readout.master_equation.RateEquation.rate_matrix` builds the matrix one move type at a time instead
  of one element at a time; bit-identical output, about twice as fast.
- `Device.provenance["phantom_version"]` of newly compiled devices reads `0.2.0`.

## 0.1.0 — development line, untagged (Steps 0–7)

The physics layers and the compile tier, built step by step; see `docs/steps/` and ADRs 0001–0023.
`layout` and closed-form `electrostatics` (Step 2); `disorder` with interface charges, misalignment,
wiring faults, trapping and per-device noise draws (Step 3); `charge` (constant interaction, Hubbard)
and `readout` (SET sensor, noise processes, WKB rates, master equation, lock-in, RF, single shot) (Steps
1, 4); `valley` (Step 5); `spin` (exchange, Lindblad dynamics, PSB, blockade, funnel, Elzerman, noise
spectroscopy, gate fidelity) (Step 6); `compile`, `device` and `optimum` (Step 7). The qdt-sim core was
migrated bit for bit apart from the empty-corner fix in `charge.ci` (MODEL_VERSION 2).
