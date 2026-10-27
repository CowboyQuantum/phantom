# Changelog

phantom follows semantic versioning in 0.x (ADR 0024): the minor version rises for any breaking change
to the consumer surface (`phantom.api`, recorded in `docs/api/surface.json`) or to the device file
format, the patch version otherwise. A change of output for the same parameters, seed and calls bumps
the `MODEL_VERSION` of the module that changed and is listed here. The owner tags releases.

## 0.2.0 — the consumer contract (Step 8; to be tagged `v0.2.0` by the owner)

First release a consumer can pin. Nothing below changes the output of any 0.1.0 model.

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
- Model versions: `device.session` 1, `device.truth` 1.

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
