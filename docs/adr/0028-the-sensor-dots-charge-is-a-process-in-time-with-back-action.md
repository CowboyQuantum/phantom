# 0028 — The sensor dot's charge is a process in time: a birth–death chain over its levels on the session clock, counted electrons per dwell, and back-action on the dots against its mean

**Status:** Accepted
**Date:** 2026-10-30
**Step:** 8 follow-up (issue #54, after ADR 0026 and 0027)

## Context
After ADR 0026 the sensor is a dot with barriers and baths. Its rates follow its gates. Its current,
though, was the stationary current at each sample's voltages and dot charge, so only the two qubit
dots' charge evolved in a session.

A real sensor with a barrier nearly pinched off does not give a smooth small current. Its own
occupation switches, the current arrives one electron at a time, and a slow sensor reads late or
misses a transition. The owner asked for this, with the gates setting the regime. They also chose to
include back-action, always on with dynamics: the sensor's charge shifts the dots' levels through the
mutual capacitance by which the dots shift the sensor.

Three forces shaped the design:

- **Cost.** A session pays one matrix exponential over the dots' 25-state window per sample (about
  0.2–0.6 ms). A joint chain over dots × sensor, or exact counting statistics by a counting field,
  would cost tens of times more.
- **Rates.** The sensor's rates span 10⁻⁵ to 10¹¹ Hz with the gates. No single algorithm samples
  both ends cheaply.
- **Coupling.** The coupling is strong: |β| ≈ 0.35 meV ≫ kT = 13 µeV, so a sensor electron moves a dot
  level by tens of kT.
- **Absolute charges.** The closed-form model counts the dots' charge absolutely (about 25 electrons
  each). A back-action measured from an absolute reference shifts the dots by millivolts. Seed 7 put
  the sensor 3 electrons below the charge the gates alone induce, a shift of 1.5 meV.

## Decision
With dynamics, the sensor dot's charge M is a birth–death chain over its levels on the session clock.
Each sample's current is the charge that leaves to the drain during the dwell, divided by the dwell.
The sensor acts back on the dots through −β_i·(M − M̄), where M̄ is the mean of its stationary state
at the present voltages and dot charge. Specifically:

- **Chain.** From M the sensor gains an electron at Γ_S n_S(μ_M) + Γ_D n_D(μ_M) and loses one at
  Γ_S(1 − n_S(μ_{M−1})) + Γ_D(1 − n_D(μ_{M−1})), with μ_M = E_c(M − x). The leads' occupations are
  convolved with the level's Lorentzian (`readout.broadening`). Through one level the stationary
  current is exactly ADR 0026's. The chain holds both sides of every level within the bias window
  plus 20 kT, at least both sides of the nearest level, and every charge it still carries. The other
  levels add their ADR 0026 current as a constant.
- **Counted charge.**
  - When the chain makes at most 200 jumps in a dwell, the counted charge and the end state come from
    an exact trajectory.
  - Otherwise the count is Gaussian. Its mean, given the start state, comes from one augmented matrix
    exponential. Its variance is the long-time variance of full counting statistics
    (1ᵀJ₂p − 2·1ᵀJ₁W^D J₁p), and the end state is drawn from e^{WT}.
  - In `mean` mode the distribution is propagated and the expected count is read.
  - Once the dwell spans more than 50 of the chain's relaxation times, the end state is the stationary
    p and the integral of e^{Ws}q is T·p + W^D(p − q), in closed form.
- **Order within a sample.** The dots step first, with the back-action of the sensor's charge at the
  start of the dwell. The sensor then steps at the dots' new charge and the sample's switcher shift.
  A `hold` evolves both, at the last bias given.
- **Back-action.** The shift is −β_i·(M − M̄) per electron on dot i, with the β the sensor reads
  (wiring asymmetry included). The closed-form potential has no electrons and stands for the sensor
  at its mean. What it leaves out, and what acts back, is the sensor's switching and its lag. Within
  a dwell the dots' rates keep detailed balance in the shifted energies.
- **Surface.**
  - `Measurement.sensor_charge`: the sensor's charge at the end of each dwell, or None without
    dynamics or for a typed sensor.
  - `RateEquation.rate_matrix(..., level_shift_ev)`.
  - The new module `readout.sensor_dynamics`, with its own random stream [seed, 8, …, 2].
  - Without dynamics, or with a typed sensor, the output is unchanged. `device.session`
    MODEL_VERSION 4. Release 0.3.0.

## Consequences
- Easier:
  - A consumer closes the sensor's barriers, or lowers a bath gate towards its turn-on, and reads
    single-electron bursts and the sensor's own telegraph.
  - A slow sensor misses a dot transition that a fast one reads cleanly (`session-slow-sensor`).
  - The shot noise of a fast sensor shrinks as 1/√dwell, with the Fano factor of the chain.
  - A sensor that lags or switches shifts the dots by its coupling.
- Harder:
  - The splitting is first order in the dwell. Within one dwell the sensor sees the dots' charge at
    its end, and the dots see the sensor's charge at its start. A dot jump therefore reaches the
    sensor's rates in the sample where it happens, and the sensor's lag, not a sample's delay, is
    what makes a slow sensor late.
  - Putting the broadened occupations into sequential rates is phenomenological. It is exact in its
    stationary mean and in the thermal limit, but its noise is the sequential one.
  - The chain has no joint energy with the dots: detailed balance holds per step, not over the joint
    process.
  - At ΓT ≈ 10⁸ the augmented exponential keeps only about eight digits, and they differ between
    platforms. The first Windows comparison differed by 1.5·10⁻⁹. Beyond 50 relaxation times,
    `propagate` therefore uses the relaxed closed form T·p + W^D(p − q) (ADR 0014).
  - A session with a sensor with leads costs about 0.3–0.7 ms per sample, against 0.2.
- Must now: a consumer that reads `n_present` alone still gets the dots. `sensor_charge` counts on the
  level numbering of φ_s, so only its changes are physical.

## Alternatives considered
- **A joint chain over dots × sensor states.** Rejected for cost: about 50–75 states and a matrix
  exponential per sample, plus counting fields for the current.
- **A back-action reference at the induced charge (x₀ + ½, or x(N) + ½).** Rejected: with the model's
  absolute dot charge it shifts the dots by millivolts (seed 7: 1.5 meV) or adds a (β·N)² term.
- **Exact trajectories at every rate.** Rejected: a fast sensor makes 10⁸ jumps in a millisecond.
- **The Gaussian count at every rate.** Rejected: it cannot show single-electron bursts.
- **Back-action off, or opt-in.** Rejected by the owner, who chose it always on with dynamics.
- **A sample's current as the instantaneous stationary current of the sensor's state.** Rejected: it
  shows neither bursts nor shot noise, and its long-dwell mean is not the counted charge.
