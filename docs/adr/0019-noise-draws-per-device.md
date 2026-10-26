# 0019 — Noise is drawn per device: a 1/f amplitude, a Poisson number of dominant switchers, a white floor

**Status:** Accepted
**Date:** 2026-10-26
**Step:** 4

## Context
Step 1 typed one set of noise processes for every device (150 µV rms 1/f, one 4 mV switcher at 0.5 Hz,
15 pA white). A phantom instance should carry its own noise as a draw (ADR 0002): devices differ in
their 1/f amplitude by a factor of a few, in how many dominant two-level fluctuators sit near the
sensor, and in those fluctuators' rates over decades. The catalog entries that will set these
distributions (#1 oxide traps and fluctuators, #12 slow drift) are not written yet.

## Decision
- `disorder.noise_draw` draws, per device and seed: the 1/f amplitude (log-normal around a material
  median), a Poisson number of dominant switchers with log-normal amplitudes and log-uniform up and
  down rates, and the white floor. Severity scales the amplitude and the switcher count; severity 0
  leaves only the white floor. The draw is part of the scenario draw and serialises with it.
- `readout.SensedDoubleDot.from_noise_draw` turns a draw into processes: the 1/f ensemble on the gates,
  the switchers summed at the sensor, the white floor rescaled to the sampling interval. `readout`
  imports `disorder` (downward); `disorder` holds only numbers.
- The 1/f band reaches down to 0.1 mHz so that cuts repeated over hours drift (`drift-repeated-cuts`);
  the dominant switcher's dwell-time statistics are the `switcher-histogram` protocol.
- `plunger-vs-barrier-map` (Step 1) is re-rendered on the master equation with WKB lead rates per
  barrier row (ADR 0017). The lock-in dynamics stay on the Step 1 rate-landscape placeholder for now:
  one master-equation step per half-cycle would cost minutes per map; this is recorded as a gap.
- All numbers are placeholders (`calibration: none`) until the catalog entries and the charge-noise
  reference are confirmed; the Step 4 P1 row on the charge-noise spectrum then fixes the median.

## Consequences
- Easier: noise realism is a draw like disorder; the SNR and histogram protocols vary by device.
- Harder: one more set of placeholders awaiting the owner; the lock-in placeholder lingers.
- Must now: `readout.loading` keeps `LeadRateLandscape` and `latched_occupation` for the lock-in only,
  and says so; `LeadBarrier` and the old plunger-vs-barrier path are gone.

## Alternatives considered
- One typed noise set per material — rejected: ADR 0002.
- A physical trap ensemble with capture/emission energetics — deferred to when entry #1 exists.
