# Protocol roster — proposed, for the owner to edit

A candidate list of measurement protocols for the gallery (ADR 0008), gathered from the usual course of
bringing up a gate-defined Si/SiGe or Si-MOS spin-qubit device: electrical health, forming dots and a
sensor, double-dot charge stability, noise and drift, spin readout and control. Every row is `proposed`
until the owner has reviewed it. Rows become one `protocols/<id>.yaml` each when their step opens.

Conventions: `P` plunger, `B` barrier, `S` sensor plunger, `V_sd` source–drain bias, `B_ext` external
field, `f` microwave frequency, `t` time. Readout channels: `current`, `sensor_conductance`,
`rf_response`, `charge_state`, `spin_outcome`. "Step" is where the physics first exists; a protocol may be
re-rendered at every later step. Any protocol can be rendered for a clean instance and for N disorder
draws side by side once Step 3 exists; that is a render option, not a separate protocol.

## A. Electrical health (layout, electrostatics) — Step 2
| id | sweep | readout | must show |
|---|---|---|---|
| `gate-leakage-iv` | each gate vs. ground, wide V range | leakage current | ~0 inside the safe window; onset of leakage at the extremes; an open gate shows no accumulation anywhere |
| `ohmic-iv` | V_sd, all gates accumulated | current | linear; series resistance |
| `turn-on-global` | all gates together | current | threshold, then saturation; MOS vs. SiGe threshold scale |
| `turn-on-single-gate` | one P, others fixed accumulated | current | threshold and saturation for that gate |
| `pinch-off-barrier` | one B, others fixed | current | exponential pinch-off; pinch-off voltage |
| `pinch-off-crosstalk` | B_i pinch-off repeated vs. V of gate j | pinch-off voltage of i | straight line; slope is the cross-capacitance ratio |
| `gate-hysteresis` | one gate up then down | current | overlap for a clean device; a loop when charge trapping is present (Step 3) |

## B. Single dot and sensor formation (charge, readout) — Step 1
| id | sweep | readout | must show |
|---|---|---|---|
| `coulomb-oscillations-single-dot` | P with barriers set | current | periodic peaks; spacing from the addition energy |
| `diamonds-single-dot` | (P, V_sd) | current or dI/dV | Coulomb diamonds; slopes equal the lever arms; height equals the addition energy |
| `barrier-barrier-map` | (B1, B2) at fixed P | current | the region where a dot forms; diagonal Coulomb lines; pinch-off corners |
| `sensor-coulomb-peak` | S | sensor_conductance | one peak; the steep flank where the sensor is operated |
| `sensor-compensation-map` | (S, P of a dot) | sensor_conductance | tilted peak lines broken at the dot's charge transitions; slope gives the compensation |

## C. Double-dot charge stability (charge, readout) — Step 1, extended in Step 4
| id | sweep | readout | must show |
|---|---|---|---|
| `stability-coarse-2dot` | (P1, P2), wide | sensor_conductance or charge_state | honeycomb; transitions stop toward (0,0); lines of both dots visible with different slopes |
| `stability-fine-2dot` | (P1, P2) around one pair of triple points | sensor_conductance | interdot line; broadening grows with t_c; the two triple points |
| `stability-derivative` | as above, rendered as dS/dP | sensor_conductance | transitions as sharp lines; a render variant worth its own sign-off |
| `interdot-detuning-cut` | ε across (1,0)–(0,1) at fixed P1+P2 | sensor_conductance | a step whose width gives t_c (and electron temperature at small t_c) |
| `tunnel-coupling-vs-barrier` | detuning cut repeated vs. B_interdot | step width | exponential t_c(V_B) |
| `detuning-cut-vs-temperature` | detuning cut vs. electron temperature | step width | thermal broadening floor; where t_c dominates |
| `bias-triangles` | (P1, P2) at finite V_sd | current | triangle pairs; edges give lever arms; excited-state lines inside |
| `stability-vs-sweep-rate` | (P1, P2) at several sweep rates | sensor_conductance | latching and missing transitions as the rate exceeds the lead rate (Step 4) |
| `lead-rate-time-trace` | t at a transition, after a gate step | sensor_conductance | exponential tunnelling into or out of the dot; the lead rate (Step 4) |
| `crosstalk-lines` | (P_i, P_j) for every pair incl. barriers | sensor_conductance | transition-line slopes; the full cross-capacitance (virtual-gate) matrix |

## D. Noise and drift (readout, disorder) — Step 1, extended in Steps 3 and 4
| id | sweep | readout | must show |
|---|---|---|---|
| `noise-time-trace` | t at a fixed sensor flank, seconds to minutes | sensor_conductance | white noise, telegraph switching, slow drift, each visible |
| `noise-psd` | derived from the trace | PSD | 1/f slope and amplitude at 1 Hz; white floor; Lorentzian bumps from fluctuators |
| `drift-repeated-cuts` | one P cut repeated over hours | transition position vs. t | slow drift and occasional jumps |
| `stability-repeat-diff` | stability map twice, hours apart | difference | rigid shifts vs. local changes; which lines moved |
| `switcher-histogram` | long trace | dwell-time histograms | two-level statistics; switching rates |

## E. Readout hardware (readout) — Step 4
| id | sweep | readout | must show |
|---|---|---|---|
| `rf-resonance` | carrier f | rf_response amplitude and phase | the tank resonance; its shift with sensor conductance |
| `rf-sensor-flank` | S at the resonance | rf_response | the same flank as the DC peak, in reflected signal |
| `single-shot-histogram` | many repeats at a readout point | counts vs. signal | bimodal histogram; separation vs. integration time |
| `snr-vs-integration-time` | integration time | SNR | square-root growth until 1/f takes over |

## F. Valley and excited states (valley, charge) — Step 5
| id | sweep | readout | must show |
|---|---|---|---|
| `magnetospectroscopy-transition` | charge-transition position vs. B_ext | sensor_conductance | kink where the Zeeman energy equals the valley splitting |
| `pulsed-gate-spectroscopy` | pulse amplitude vs. P | sensor_conductance | excited-state lines; orbital and valley spacings |
| `valley-splitting-draws` | N disorder draws | valley splitting | the distribution across dots and the fraction below a usable threshold |
| `t1-hot-spot` | relaxation time vs. B_ext | T1 | the dip where the Zeeman and valley energies cross |

## G. Spin readout and control (spin) — Step 6
| id | sweep | readout | must show |
|---|---|---|---|
| `psb-bias-triangles` | (P1, P2) at finite V_sd across (1,1)–(0,2) | current | current suppressed at the triangle base; lifted at B_ext = 0 or by mixing |
| `psb-readout-window` | ε after loading a random spin state | sensor_conductance | the detuning window where singlet and triplet differ; latched vs. unlatched variants |
| `spin-funnel` | (ε, B_ext) with a short wait | singlet probability | the funnel of S–T₋ mixing |
| `elzerman-trace` | t after pulsing to the read level | charge_state | spin-dependent tunnel-out blips; the read window |
| `esr-line` | (B_ext, f) | spin_outcome | one resonance line per dot; slope gives the g-factor; splitting between dots |
| `rabi-chevron` | (f − f₀, burst time) | spin_outcome | the chevron; Rabi frequency vs. drive |
| `ramsey-fringes` | free evolution time at fixed detuning | spin_outcome | fringes with a decaying envelope; T₂* |
| `exchange-oscillations` | (ε, evolution time) | singlet probability | oscillation frequency rising toward zero detuning; J(ε) |
| `exchange-vs-barrier` | J extracted vs. B_interdot | J | exponential |

## Candidates not included, for the owner to decide
- Multi-dot arrays (linear 3+ dots): per-pair stability maps and sensor coverage. Deferred until a layout step defines arrays.
- Charge-noise spectroscopy through the qubit (T₂ echo vs. filter order). Deferred to Step 6 or later.
- Anything for germanium holes (ADR 0002 lists it as later): g-tensor anisotropy, spin–orbit field.
