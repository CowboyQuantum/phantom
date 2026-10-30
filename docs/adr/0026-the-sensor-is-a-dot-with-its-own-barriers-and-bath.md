# 0026 — The sensor is a dot with its own barriers and bath: compiled leads, a lifetime-broadened current, its charging energy and coupling from the layout

**Status:** Accepted
**Date:** 2026-10-29
**Step:** 7 (issue #47, after the Step 8 release)

## Context
The charge sensor was a potential minimum under one gate (SG) with nothing on either side of it: in the
textbook layout the potential rose back to the ungated level, so the sensor had no barriers and no
leads, and nothing could tune its coupling. Its Coulomb-peak shape was qdt-sim's typed SET: a peak
period of 0.25 V of a "V-equivalent" sensor potential φ_s with a lever arm of 4·10⁻³ eV/V, a thermal
energy of 80 µeV, a peak current of 1 nA and a per-electron coupling of −4 and −3 mV. Its peaks were
about 1 V apart on SG, a hundred times a real sensor dot's spacing.

The owner asked for barriers and a bath that can be tuned. Giving the typed sensor rates alone does not
work: compile moves SG onto the flank of the nearest peak, up to half a period (0.5 V) away, and a
sensor dot whose barriers come from the landscape is pinched off or wide open after such a move (rates
from 10⁻²⁰ to 10¹² Hz in a scan). Real barriers only make sense on a sensor whose peaks are as far apart
as its charging energy says. The closed-form potential has no absolute Fermi level and no screening
(ADR 0013), so the absolute energies of a dot and a lead do not say where electrons tunnel; the
qubit dots' lead rates use the higher of the dot's ground level and the lead's band bottom plus 2 meV,
and are constant at runtime. The noise draw's switcher amplitudes (ADR 0019) are stated in the typed
sensor's φ_s units.

## Decision
Give the textbook layout a sensor channel of the qubit channel's pattern (barrier gates SB_L, SB_R and
bath gates SR_L, SR_R), and make a compiled sensor with leads a dot like the qubit dots: its rates
compiled from its barriers with their gate slopes and an attempt-rate ceiling, its baths conducting above
an accumulation threshold, its current the lifetime-broadened single-level sum, and its charging energy,
coupling and temperature from the layout and the device. Specifically:

- **Geometry.** The four gates are listed after the first eight, which keep their order; the sensor's
  search region spans its barriers. A layout without the named sensor reservoirs keeps the typed sensor.
- **Charging energy and coupling.** From the three-dot capacitance matrix (the double dot's 2×2 block,
  the sensor's total capacitance by the qubit dots' rule, sphere mutual capacitances):
  E_c,s = e²(C⁻¹)_ss and β_i = −(e/1 V)(C⁻¹)_si. φ_s is the potential at the dot in volts (lever arm
  1 eV/V), the peak period is E_c,s, the thermal width 2kT at the device temperature, and the wiring
  draw's asymmetry multiplies β (ADR 0016). The sensor's charge does not act back on the dots.
- **Rates.** Each lead's rate is f_attempt·exp(−2S) through its barrier to the centroid of its bath gate,
  at a tunnelling level placed once by the qubit dots' rule at the wired programmed point and held fixed:
  on a Coulomb peak electrons tunnel at the leads' electrochemical potential, which no gate moves. The
  slopes d ln Γ/dV against every gate are central differences of ±10 mV at that level. At runtime the
  exponent moves with the slopes and a logistic, anchored to give the compiled rate exactly, saturates it
  at the attempt rate.
- **Baths.** A bath conducts with the occupation of a soft step where the electron energy at its gate's
  centroid (linear in the gate voltages, with the draw's offset) crosses −50 meV over 4 meV, the turn-on
  protocols' placeholder (calibration none); it multiplies that lead's rate.
- **Current.** I = e·Γ_SΓ_D/(Γ_S + Γ_D)·Σ_k [n_γ(μ_k − eV/2) − n_γ(μ_k + eV/2)] with n_γ the Fermi function
  convolved with a Lorentzian of full width γ = ħ(Γ_S + Γ_D), in closed form through the complex digamma
  function (`readout.broadening`, entries `beenakker1991-coulomb-blockade` and
  `foxman1993-quantum-levels-coulomb-island`, both `proposed`). No typed background: between peaks the
  current is the Lorentzian tails. Spin degeneracy is not counted.
- **Set point.** One step onto the flank and a re-derivation, as before, then one more step with the
  re-derived sensor and its leads, without deriving again, the peak width taken at the step's own target
  (three fixed passes). The sensor's linearisation and leads are the ones the set point used; the double
  dot is derived again at the final point.
- **Noise.** The draw's switcher amplitudes are multiplied by the ratio of the sensor dot's mean coupling
  to the typed sensor's, so a switcher keeps its size against a dot electron's signal.
- **Consumers.** `State.sensor_gamma_hz`; `Knowable.sensor_barrier_gates`, `sensor_reservoir_gates`;
  `SensorTruth.lifetime_width_v`, `gamma_source_hz`, `gamma_drain_hz`, with the wing and slope of the
  broadened lineshape (the steepest point, by bisection on Re ψ‴); `SensorLeads` in `phantom.api`;
  `DEVICE_FORMAT_VERSION` 2; release 0.3.0.

## Consequences
- Easier: a consumer tunes the sensor as on a device: a bath turn-on, a barrier pinch-off, a
  barrier–barrier map with a diagonal Coulomb line, peaks tens of mV apart on SG that grow and then
  broaden towards e²/h as the barriers open; the sensor's charging energy and coupling come from the same
  geometry as the qubit dots'; the truth states the rates and the broadened wing.
- Harder: the rates inherit the closed-form potential's sensitivity, a factor e every 3–5 mV on a sensor
  barrier or SG, and a wide bath gate reaches under the adjacent barrier about as strongly as the
  barrier gate itself; the runtime exponential stays within 20 % of the compile tier over ±10 mV and
  drifts beyond (a factor 12 at −30 mV on SG). A sharp, strongly coupled sensor makes raw stability maps
  look like a real uncompensated sensor's, and the gate 1/f noise now moves it visibly. The common
  programmed point's sensor barriers are set for seed 2: four of seeds 0–11 form no sensor dot there.
- Must now: qdt moves its pin to 0.3.0 and adds the four gates; device files from 0.2.0 do not load
  (recompile from the seed); a consumer tunes SG before reading, since the programmed SG alone usually
  sits between peaks; screening (plan §8 #1) is what would temper the sensitivity.

## Alternatives considered
- Rates on the typed sensor, its peak shape kept — rejected: the set point's half-period move of SG
  (0.5 V) pinches the barriers or removes them, so the compiled rates are a near-random function of the
  draw; the typed 80 µeV thermal energy and the device's 13 µeV also disagree.
- Rates constant at runtime, as the qubit dots' leads — rejected: the request is barriers and a bath that
  tune; the sensor's rates are what a sensor tune-up moves.
- The tunnelling level following the dot's ground state at every gate point — rejected: raising SG then
  closes the sensor and the compile tier and runtime disagree in sign between devices; the level pinned
  at the leads' electrochemical potential is the physical choice for a many-electron sensor on a peak.
- Sequential tunnelling without lifetime broadening — rejected: opening the barriers would raise the
  current without bound, past the conductance quantum.
- A bath as a series resistance — rejected for now: it needs a 2DEG mobility and a self-consistent bias
  division; the occupation step uses the threshold the turn-on protocols already use.
- Iterating the set point with re-derivations to a tolerance — rejected: the sensor potential α·V jumps
  by several peak widths when a re-derivation moves the dot minimum by a grid step, so it does not
  converge, and a tolerance test differs between platforms (ADR 0014).
