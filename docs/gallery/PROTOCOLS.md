# Protocol roster — proposed, for the owner to edit

A candidate list of measurement protocols for the gallery (ADR 0008), gathered from the usual course of
bringing up a gate-defined Si/SiGe or Si-MOS spin-qubit device: electrical health, forming dots and a
sensor, double-dot charge stability, noise and drift, spin readout and control. Every row is `proposed`
until the owner has reviewed it. Rows become one `protocols/<id>.yaml` each when their step opens.

Conventions: `P` plunger, `B` barrier, `S` sensor plunger, `V_sd` source–drain bias, `B_ext` external
field, `f` microwave frequency, `t` time, `ε` detuning, `τ` evolution or burst time. Readout channels:
`current`, `sensor_conductance`, `sensor_lockin` (demodulated sensor response to a square-wave
modulation of a named gate, with modulation amplitude and frequency as parameters), `rf_response`,
`charge_state`, `spin_outcome`. "Step" is where the physics first exists; a protocol may be re-rendered
at every later step. Any protocol can be rendered for a clean instance and for N disorder draws side by
side once Step 3 exists; that is a render option, not a separate protocol. A numerical derivative of any
map is likewise a render option; the hardware derivative (`sensor_lockin`) is a protocol of its own
because it has its own physics (finite modulation amplitude, noise at the modulation frequency).

## A. Electrical health (layout, electrostatics) — Step 2
| id | sweep | readout | must show |
|---|---|---|---|
| `potential-map-layout` | (x, y) in the 2DEG plane at a named gate setting, gate polygons overlaid | potential energy (meV) | two wells and the barriers between them; the dot positions the layout implies |
| `potential-cut-along-channel` | position along the dot-to-dot axis, repeated for several barrier voltages | potential energy (meV) | the barrier height falling as the gate opens; the well depths |
| `gate-leakage-iv` | each gate vs. ground, wide V range | leakage current | ~0 inside the safe window; onset of leakage at the extremes; an open gate shows no accumulation anywhere |
| `ohmic-iv` | V_sd, all gates accumulated | current | linear; series resistance |
| `turn-on-global` | all gates together | current | threshold, then saturation; MOS vs. SiGe threshold scale |
| `turn-on-single-gate` | one P, others fixed accumulated | current | threshold and saturation for that gate |
| `turn-on-single-gate-layout` | one barrier swept with every other gate open, from the layout | current | the channel opens where the energy under that barrier, the highest point along the channel, crosses the threshold; threshold and slope from the lever arms (Step 2) |
| `pinch-off-barrier` | one B, others fixed | current | exponential pinch-off; pinch-off voltage |
| `pinch-off-crosstalk` | B_i pinch-off repeated vs. V of gate j | pinch-off voltage of i | straight line; slope is the cross-capacitance ratio |
| `gate-hysteresis` | one gate up then down | current | overlap for a clean device; a loop when charge trapping is present (Step 3) |

## B. Single dot and sensor formation (charge, readout) — Step 1
| id | sweep | readout | must show |
|---|---|---|---|
| `coulomb-oscillations-single-dot` | P with barriers set | current | periodic peaks; spacing from the addition energy; peak width from temperature and lead coupling |
| `diamonds-single-dot` | (P, V_sd) | current or dI/dV | Coulomb diamonds; slopes equal the lever arms; height equals the addition energy |
| `barrier-barrier-map` | (B1, B2) at fixed P | current | the region where a dot forms; diagonal Coulomb lines; pinch-off corners |
| `plunger-vs-barrier-map` | (P1 fast, lead barrier BL of dot 1 slow); dot 1 is the only dot with a reservoir, dot 2 loads only through dot 1; P2 fixed | sensor_conductance, plus the charge present and the charge at equilibrium | the dot-1 transition lines tilt with the lead barrier's lever arm; as the barrier closes the tunnel rate drops below the inverse dwell time and the lines fade, latch along the row and carry across the flyback; the interdot-barrier variant is a separate protocol when Step 4 adds tunnel-coupling broadening (owner's configuration, 2026-10-25) |
| `sensor-coulomb-peak` | S | sensor_conductance | one peak; the steep flank where the sensor is operated |
| `sensor-compensation-map` | (S, P of a dot) | sensor_conductance | tilted peak lines broken at the dot's charge transitions; slope gives the compensation |

## C. Double-dot charge stability (charge, readout) — Step 1, extended in Step 4
| id | sweep | readout | must show |
|---|---|---|---|
| `stability-coarse-2dot` | (P1, P2), wide | sensor_conductance or charge_state | honeycomb; transitions stop toward (0,0); lines of both dots visible with different slopes |
| `stability-coarse-2dot-layout` | as above, but every parameter except the sensor's per-electron coupling derived from a layout file; the sensor gate set on its flank first | sensor_conductance | the same honeycomb, with period e/C_gate and slopes from the derived cross lever arms (Step 2) |
| `stability-fine-2dot` | (P1, P2) around one pair of triple points | sensor_conductance | interdot line; broadening grows with t_c; the two triple points |
| `stability-lockin-p1` | (P1, P2) with a 0.5 mV square wave on P1 | sensor_lockin, plus the lead-rate map | dot-1 lines as signed lines; where the lead rate is below the chop frequency the charge cannot follow and lines latch or vanish; where it is far above, lines are lifetime-broadened and fainter |
| `stability-lockin-p2` | same with the chop on P2 | sensor_lockin, plus N1 present | dot-2 lines strong, dot-1 lines weak |
| `stability-lockin-common` | same with the chop on P1 and P2 together | sensor_lockin | both dots' lead lines; the interdot line weak |
| `stability-lockin-detuning` | same with the chop on P1 up and P2 down | sensor_lockin, plus N1 present | the interdot line strongest; lead lines of the two dots with opposite sign |
| `interdot-detuning-cut` | ε across (1,0)–(0,1) at fixed P1+P2 | sensor_conductance | a step whose width gives t_c (and electron temperature at small t_c) |
| `tunnel-coupling-vs-barrier` | detuning cut repeated vs. B_interdot | step width | exponential t_c(V_B) |
| `detuning-cut-vs-temperature` | detuning cut vs. electron temperature | step width | thermal broadening floor; where t_c dominates |
| `bias-triangles` | (P1, P2) at finite V_sd | current | triangle pairs; edges give lever arms; excited-state lines inside |
| `stability-vs-sweep-rate` | (P1 slow, P2 fast) at several dwell times, lead rates from the layout's WKB barriers | sensor_conductance | the equilibrium honeycomb at long dwell; lines fade and shift along the fast axis as the dwell falls below one over the lead rate; frozen at the shortest (Step 4) |
| `stability-direction` | the same map at a lagging dwell, fast axis forward and backward | sensor_conductance, charge difference | latched transitions on opposite sides of their equilibrium position: the hysteresis of a latched diagram (Step 4) |
| `stability-direction` | the same map at a lagging dwell, fast axis forward and backward | sensor_conductance, charge difference | latched transitions on opposite sides of their equilibrium position: the hysteresis of a latched diagram (Step 4) |
| `lead-rate-time-trace` | t after a plunger step across a transition, for several lead-barrier settings | charge on the dot | exponential relaxation whose rate is the WKB lead rate the layout predicts; decades of rate over a few hundred mV of barrier (Step 4) |
| `crosstalk-lines` | (P_i, P_j) for every pair incl. barriers | sensor_conductance | transition-line slopes; the full cross-capacitance (virtual-gate) matrix |

## D. Noise and drift (readout, disorder) — Step 1, extended in Steps 3 and 4
| id | sweep | readout | must show |
|---|---|---|---|
| `noise-time-trace` | t at a fixed sensor flank, seconds to minutes | sensor_conductance | white noise, telegraph switching, slow drift, each visible |
| `noise-psd` | derived from the trace | PSD | 1/f slope and amplitude at 1 Hz; white floor; Lorentzian bumps from fluctuators |
| `charge-noise-sensor-referred` | PSD at the sensor flank, divided by the flank slope and by the sensor's lever arm to a dot gate | equivalent gate-voltage noise (V/√Hz) and energy noise (eV/√Hz) vs. f | the charge-noise amplitude at 1 Hz and the 1/f exponent, in gate-referred units that can be compared with published values and with the qubit-side `charge-noise-spectroscopy-cpmg`; the same disorder draw must give both |
| `drift-repeated-cuts` | one P cut repeated over hours | transition position vs. t | slow drift and occasional jumps |
| `stability-repeat-diff` | stability map twice, hours apart | difference | rigid shifts vs. local changes; which lines moved |
| `switcher-histogram` | long trace | dwell-time histograms | two-level statistics; switching rates |

## E. Readout hardware (readout) — Step 4
| id | sweep | readout | must show |
|---|---|---|---|
| `rf-resonance` | carrier f, several sensor resistances | rf_response amplitude and phase | the tank resonance; the dip to zero at the matched resistance; the phase flipping by half a turn across the match (Step 4) |
| `rf-sensor-flank` | S at the resonance | DC current, sensor resistance, reflected amplitude | the same flank as the DC peak in reflected signal, steepest where the resistance crosses the match (Step 4) |
| `single-shot-histogram` | many consecutive shots at a readout point | counts vs. integrated signal, shot sequence | bimodal histogram, one electron apart; both peaks moved by 1/f and telegraph noise between shots (Step 4) |
| `snr-vs-integration-time` | integration time | SNR, fidelity | square-root growth while white noise dominates, saturation and decline once slow noise moves the distributions (Step 4) |

## F. Valley and excited states (valley, charge) — Step 5
| id | sweep | readout | must show |
|---|---|---|---|
| `magnetospectroscopy-transition` | charge-transition position vs. B_ext | sensor_conductance | kink where the Zeeman energy equals the valley splitting |
| `pulsed-gate-spectroscopy` | pulse amplitude vs. P | sensor_conductance | excited-state lines; orbital and valley spacings |
| `valley-splitting-draws` | N disorder draws | valley splitting | the distribution across dots and the fraction below a usable threshold |
| `t1-hot-spot` | relaxation time vs. B_ext | T1 | the dip where the Zeeman and valley energies cross |

## G. Spin readout, exchange and control (spin) — Step 6
| id | sweep | readout | must show |
|---|---|---|---|
| `psb-bias-triangles` | (P1, P2) at finite V_sd across (1,1)–(0,2) | current | current suppressed at the triangle base; lifted at B_ext = 0 or by mixing |
| `psb-readout-window` | ε after loading a random spin state | sensor_conductance | the detuning window where singlet and triplet differ; latched vs. unlatched variants |
| `spin-funnel` | (ε, B_ext) with a short wait | singlet probability | the funnel of S–T₋ mixing |
| `elzerman-trace` | t after pulsing to the read level | charge_state | spin-dependent tunnel-out blips; the read window |
| `esr-line` | (B_ext, f) | spin_outcome | one resonance line per dot; slope gives the g-factor; splitting between dots |
| `rabi-chevron` | (f − f₀, τ) | spin_outcome | the chevron; Rabi frequency vs. drive |
| `ramsey-fringes` | τ at fixed detuning | spin_outcome | fringes with a decaying envelope; T₂* |
| `exchange-fingerprint` | (ε, τ) after a singlet is separated into (1,1) | singlet probability | the fingerprint: oscillations whose frequency J(ε)/h rises toward the (0,2) anticrossing and whose visibility decays with charge noise, fastest where dJ/dε is largest |
| `exchange-house-cells` | (P1, P2) around the evolve coordinate, with J held at a high value by the interdot barrier, at fixed τ | singlet probability | the house-cell pattern: cells of constant singlet return bounded by the charge-transition lines of the honeycomb and by contours of constant J; the evolve point sits inside a cell; cell edges move with the barrier setting |
| `exchange-vs-detuning` | J extracted from the fingerprint vs. ε | J | the J(ε) curve; its slope dJ/dε, which sets charge-noise sensitivity |
| `exchange-vs-barrier` | J extracted vs. B_interdot at the symmetric point | J | exponential |
| `charge-noise-spectroscopy-cpmg` | number of π pulses at fixed total τ, on a charge-sensitive qubit state (exchange or detuning) | coherence vs. filter frequency | S_ε(f) over the kHz–MHz decades, with the exponent of the 1/f^α law; must join the sensor `noise-psd` at low frequency with the same amplitude scaling |
| `charge-noise-spectroscopy-j-vs-ramsey` | Ramsey T₂* vs. ε (through dJ/dε) | T₂* | T₂* ∝ 1/(dJ/dε) × noise amplitude; a direct check that the same charge-noise draw governs the qubit and the sensor |

## H. Disorder draws (layout, electrostatics, disorder) — Step 3
| id | sweep | readout | must show |
|---|---|---|---|
| `potential-map-disorder` | (x, y) landscape, clean beside three seeded draws | potential energy (meV), gate outlines | the gate layer shifted and resized; charge dips and bumps (visible at MOS densities, below a meV for SiGe) |
| `stability-coarse-2dot-disorder` | (P1, P2) coarse map, clean beside three draws, each re-derived | sensor_conductance | lines tilt and shift, periods change between draws; an extra family of lines when a draw has an unintentional minimum |
| `unintentional-dots` | charge-density severity, many draws | extra minima in the channel (count) | the mean and 90th percentile of extra minima and the fraction of affected draws rising with density; the P1 row compares with a published MOS/SiGe comparison |
| `fault-signatures-stability` | (P1, P2) coarse map, programmed voltages, clean wiring and named faults | sensor_conductance | one dot's lines faint (sensor asymmetry); a barrier following a plunger (short) changes that dot's slopes and spacing; an open barrier floating at a stuck potential shifts the lines |
| `gate-leakage-pair` | one gate swept, others grounded, with a finite resistance to a neighbour | leakage current of both gates | equal and opposite linear currents on the pair on top of the dielectric leakage; the healthy curve flat |
| `gate-hysteresis-trapping` | one gate up then down with trapping | current, trapped fraction | a loop: the downward trace turns off at a shifted voltage; the trapped fraction charges above capture and releases below release |

## Candidates not included, for the owner to decide
- Multi-dot arrays (linear 3+ dots): per-pair stability maps and sensor coverage. Deferred until a layout step defines arrays.
- Anything for germanium holes (ADR 0002 lists it as later): g-tensor anisotropy, spin–orbit field.
