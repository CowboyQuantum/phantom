# 0018 — RF reflectometry as a lumped tank loaded by the sensor; single-shot readout as integrated noisy shots

**Status:** Accepted
**Date:** 2026-10-26
**Step:** 4

## Context
Step 4 owes the readout response beyond the DC sensor current: the radio-frequency reflectometry that
real setups use for fast charge readout, and the single-shot statistics (histograms, signal-to-noise,
fidelity against integration time) that decide whether a readout works. The step file's open decision
was how much of RF to model: a transfer function or the tank plus sensor impedance.

## Decision
- **RF** (`readout.rf.TankCircuit`): a series inductor into the sensor's small-bias resistance in
  parallel with a parasitic capacitance, driven through a 50 Ω line. Γ = (Z − Z₀)/(Z + Z₀) with
  Z = iωL + 1/(iωC_p + 1/R). This gives the resonance, the matching resistance L/(C_p Z₀) at which the
  reflection vanishes, and the phase flip across it, from three numbers. The sensor resistance is
  V_sd / I of the SET model, so the RF flank mirrors the DC Coulomb peak. Cables, amplifier noise and a
  tunable match are not modelled; a consumer that needs them puts a transfer function in front.
- **Single shot** (`readout.single_shot`): a shot integrates the sensor current for τ at a readout
  point with the double dot in one of two charge states; shots are consecutive in time on the device's
  noise processes, so 1/f drift and telegraph switching move the distributions between shots as they do
  in a real histogram. White readout noise keeps a fixed spectral density: its per-sample standard
  deviation scales as 1/√dt from the 15 pA stated at 1 ms sampling, so its shot contribution falls as
  1/√τ. SNR is the separation of the two means over the pooled standard deviation; fidelity is the
  fraction of shots on the correct side of the midpoint threshold.
- **Protocols**: `rf-resonance`, `rf-sensor-flank`, `single-shot-histogram`, `snr-vs-integration-time`
  (the last shows the √τ law and its breakdown when slow noise takes over).

## Consequences
- Easier: RF and DC readout share one sensor model; the SNR curve exposes the noise model's low-frequency
  content directly, which the P1 charge-noise row will constrain.
- Harder: amplifier noise temperature, the dominant limit in practice, is absent; the SNR here is the
  device-side bound. The noise amplitudes are placeholders (`calibration: none`).
- Must now: the catalog entry for traps and fluctuators (#1) sets the telegraph and 1/f parameters the
  histogram depends on; until then the picture is a placeholder shape.

## Alternatives considered
- A transfer function only — rejected: matching and the phase flip are what a tune-up looks for, and
  they need the impedance.
- Full circuit with cable and amplifier — deferred; adds parameters without changing the device physics.
