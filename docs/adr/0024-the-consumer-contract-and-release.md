# 0024 — The consumer contract: one import surface above `compile`, a session with a forward clock, truth apart from what a user may see, a contract test consumers run, a checked surface snapshot, and 0.x releases

**Status:** Accepted
**Date:** 2026-10-27
**Step:** 8

## Context
ADR 0003 promised consumers plan §4 (`compile`, `Device.at`, `Device.measure`, `Device.truth`,
`optimum`), a contract test the qdt adapter also runs, and a changelog. Step 7 built the compiled
`Device` and `optimum` but no runtime surface with time: `runtime_sensed()` hands out the migrated
qdt-sim core, which a consumer would have to drive, seed and clock itself. qdt's Analysis 0007
(issue #40) found six gaps between the planned surface and what its adapter needs: the sensor bias in
`sense`, a noise seed apart from the device seed, holding voltages without reading, a statement of
what one sample is, truth under the names a tune-up recovers, and the line between knowable fields and
truth. Rule 2 said nothing imports `compile`, yet a consumer's one import has to reach `compile_device`.

## Decision
- **One import surface.** `phantom.api` re-exports every name a consumer uses (`compile_device`,
  `CompileParams`, `Device`, `State`, `Session`, `Measurement`, `SpinOutcome`, `Truth`, `SensorTruth`,
  `Knowable`, `optimum` and its parts, the layout types, `ValleyParams`, `SpinParams`,
  `DEVICE_FORMAT_VERSION`, `__version__`) and adds `model_versions()`. `phantom.api` and
  `phantom.contract` sit above `compile`: they may import every layer, `contract` imports `api`, and no
  layer or root module imports either (rule 2 amended; `scripts/check_imports.py`).
- **The runtime surface lives in `device`** (`device.session`, `device.truth`), so it reads the record
  and the layers below and never `compile`. `Device.at(v, v_sd)` is the noise-free `State`, vectorised.
  `Device.session(noise, dynamics, charge, noise_seed)` is a `Session` with a clock that only moves
  forward: `sense(v, dwell_s, samples, v_sd)` reads each point at the end of each dwell, sample k at
  t₀ + (k+1)·dwell; `hold(v, duration_s)` keeps the gates at `v` and evolves the charge there;
  `wait(duration_s)` holds at the last point; `spin_readout` samples the PSB readout of a prepared
  pair as a protocol step (voltage-driven spin is #41). `Device.measure` is a fresh session's `sense`.
- **What one sample is**, fixed per session: `charge="sampled"` (default) puts one charge state in the
  sensor per sample, the state of a Markov jump trajectory of the rate equation (from the present
  state, the next is drawn from the column of exp(W·dwell)) with `dynamics=True`, an independent
  Boltzmann draw (qdt-sim's behaviour) without; `charge="mean"` puts ⟨N⟩ in, propagated or thermal.
  The noise processes are sampled at each sample's time; the white floor keeps the draw's density at
  1 ms (σ ∝ √(1 ms / dwell)). A hold samples the gate noise once at its start.
- **Seeds.** The device seed is the device's identity. `noise_seed=None` replays the device's own noise
  realisation, so identical calls on two fresh sessions are identical; an integer replaces only the
  noise draw's realisation seed (a hash of both seeds) and keys the session's charge and spin streams,
  keeping every noise parameter.
- **Knowable and truth.** `Device.knowable()` is what a user may be shown: material, gate names, layout,
  temperature, programmed point, sensor, barrier and reservoir gates. Everything else on a `Device` is
  truth, including the seed and `operating_point_v` (it holds the compiled sensor set point).
  `Device.truth()` gathers it under the quantities a tune-up recovers: the charge lattice in electrons
  per programmed volt (through the wiring matrix) and the sensor peak in programmed sensor-gate volts
  (peak, period, charging energy, lever arm, width, the wing at atanh(1/√3) widths and the linear-
  response sensitivity there at a given bias), so qdt's adapter only renames.
- **Contract test.** `phantom.contract.run(device)` checks every promise above on a device;
  `python -m phantom.contract [device.json ...]` runs it on the reference device and on each file, exit
  status 0 or 1. phantom's CI runs it (`tests/test_contract.py`); qdt's CI runs it on the pinned tag
  and its committed device files.
- **Surface snapshot.** `scripts/api_surface.py` writes `docs/api/surface.json` (every public name, its
  signature, dataclass fields, public methods); `check_all.sh` fails when the surface differs from it,
  so an API change is a reviewed diff with a `CHANGELOG.md` entry.
- **Releases.** Semantic versioning in 0.x: the minor version rises for any breaking change to the
  surface or to `DEVICE_FORMAT_VERSION`, the patch version otherwise; a change of output for the same
  parameters, seed and calls bumps the module's `MODEL_VERSION` and is listed in the changelog. The
  owner tags releases (`v0.2.0` first); how often phantom releases and qdt moves its pin is the owner's
  open decision (Step 8).

## Consequences
- Easier: qdt's adapter becomes renaming and plumbing over `Session` (its `hold`/`sense` protocol maps
  one to one); a consumer gets a new noise realisation on the same committed device; truth and the
  knowable fields are stated by phantom, not re-derived by a consumer; an API change cannot slip in
  unreviewed.
- Harder: every name in `phantom.api` is a promise; a session with dynamics costs about 1 ms per
  sample (one matrix exponential), so long maps belong in `dynamics=False` or `charge="mean"`; a point
  evaluated alone and inside a batch agree within 1e-12 relative, not bit for bit (numpy's kernels).
- Must now: the owner tags v0.2.0 after this merges; the qdt side (#140, #141, #142) pins it; the
  compiled rates are constant at runtime (`State` says so) until a later step makes them voltage-
  dependent; the wiring draw's sensor-coupling asymmetry is not yet applied by `runtime_sensor`.

## Alternatives considered
- `Device` methods only, no `phantom.api` — rejected: `compile_device` lives above `device`, so a
  consumer would import from three layers and the surface would have no single place to snapshot.
- `sense` returns the equilibrium average only — rejected: qdt samples a charge state per shot, and a
  trace at a degeneracy must switch; the average stays available as `charge="mean"`.
- A noise seed that redraws the noise parameters — rejected: the amplitudes and switchers are the
  device's draw (ADR 0019); only the realisation belongs to a run.
- Sampling the charge per sample from the propagated distribution — rejected: it loses the correlation
  between consecutive samples that a slow lead produces; the trajectory keeps it.
- `wait(duration_s, v)` instead of `hold(v, duration_s)` — `hold` names the gate point first like
  `sense`, and `wait` stays the hold at the last point.
