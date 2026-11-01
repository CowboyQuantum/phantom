# 0034 — The interdot barrier has one transmission: its screened top, continued past the merge, sets t, the interdot rate and the conductance

**Status:** Accepted
**Date:** 2026-11-01
**Step:** 8 follow-up (issue #70, with #69; after ADR 0029, 0030, 0032)

## Context
Running qdt on `v0.4.0`, the owner saw the honeycomb turn into one superdot "basically instantly" along
the plungers, at a place that did not move with B_M, and asked for the transition to happen "over a
larger space" (2026-10-31). The measurements on the reference device (#70):

- **A steep transition.** The fractional splitting f went from 0.2 to 0.9 within 15 mV along the plunger
  diagonal and 30 mV along B_M. Kemble's transmission changes by a factor e for every ħω_b/2π ≈ 0.7 meV
  of barrier height, and in the closed-form potential the barrier top moves at −152 meV/V with B_M and
  −128 and −110 meV/V with the plungers: the plungers move the interdot barrier nearly as much as the
  barrier gate does, because nothing screens them (ADR 0013, #48). In a measured device a plunger's
  lever arm on the interdot barrier is several times smaller than the barrier gate's.
- **B_M stopped acting.** Above about 0.35 V every B_M gave the same map. The shape grid (ADR 0030) ends
  where the potential stops forming two dots and holds every quantity there, the barrier top included,
  while the plungers' part of the top kept growing at about 240 meV/V as they went down. So past the
  merge node the plungers always won at low enough occupancy and the double dot re-formed at the same
  place for every larger B_M, and B_M kept adding charge through its cross-capacitance while its effect
  on the barrier was frozen. The same hold is #69: below the compiled plungers, opening B_M could not
  raise t.
- **Three models of one barrier.** t moved with half the compiled slopes of the WKB action (ADR 0029),
  the interdot rate Γ_m with those slopes (ADR 0027), the conductance g with Kemble's transmission of the
  top (ADR 0032). Along the barrier gates t and Γ_m took the shape grid's action column, which is zero
  once the top sinks below the level, so t held at about half its ceiling while g went on to 1.

Four forces shaped the design:

- **The owner's requirement.** The honeycomb-to-superdot transition spans tens of millivolts, about 50 to
  100 mV of B_M on the reference device, and a far-open B_M leaves one dot at every occupancy.
- **No solver yet.** Screening by the electrons is the solver's business (plan §8 #1). Until then a
  placeholder has to say how much of a gate's bare action reaches the barrier, as the bath threshold and
  the lead Fermi level already do (`calibration: none`).
- **One barrier, one number.** The WKB action of a parabolic top is πΔU/ħω_b, so e^(−2S) is Kemble's
  deep limit; on the reference device the compiled WKB slopes of ln Γ_m agree with 2π·dΔU/dV/ħω_b within
  10–30 % on every gate. Kemble's form also passes through T = ½ and on to 1 continuously, which the
  action cannot (ADR 0032 rejected it for g for that reason).
- **Bit-exact anchor and cross-platform digits (ADR 0014).** At the operating point the compiled t, Γ_m
  and g must equal the runtime's, and the new pieces must be closed forms: a weighted sum, a linear
  continuation, a logistic.

## Decision
The barrier between the dots has one transmission, Kemble's of its top ΔU(V), and t, Γ_m and g are
functions of it: g = T, t = (ħω₀/π)·√T, Γ_m = f_a·T. ΔU(V) is the compiled top moved by every gate through
screened slopes and, along the barrier gates, through the shape grid's `barrier_top_ev` column continued
linearly past each end. Specifically:

- **The top.** `charge.barrier.InterdotBarrier`: ΔU(V) = ΔU_ref + Σ_g σ_g·δ_g(V_g), with δ_g the shape
  grid's change along a barrier gate (PCHIP inside, the interpolant's end slope beyond:
  `ShapeGrid.quantity_along(..., beyond="linear")`) and the compiled slope times ΔV for every other gate.
- **The screening.** σ_g is `CompileParams.interdot_barrier_screening` (0.4) for the barrier gate and that
  times `interdot_cross_screening` (0.3) for every other gate; both placeholders, `calibration: none`,
  and 1 and 1 is the bare potential. `Device.interdot_barrier_screening` records σ per gate
  (`DEVICE_FORMAT_VERSION` 6); the compiled slopes stay the bare potential's.
- **One transmission.** `InterdotSplitting` reads g from the barrier. `InterdotCoupling` with a `barrier`
  gives t_max·√T. `GateRates` with an `interdot` barrier gives the third rate's exponent as
  ln(T/(1 − T)) − ln(T_ref/(1 − T_ref)) = −2π(ΔU − ΔU_ref)/ħω_b, so Γ_m = f_a·T in `saturating_rates_hz`'s
  own form. `Device.runtime_interdot_barrier()` builds the barrier; `runtime_splitting`, `runtime_coupling`
  and `runtime_rates` share it.
- **The compiled record.** `Device.tunnel_coupling_ev` is (ħω₀/π)·√T_ref, the interdot `CompiledBarrier`'s
  `rate_hz` is f_a·T_ref and its `log_rate_slopes_per_v` are −2πσ_g·dΔU/dV_g/ħω_b, and
  `barrier_v_scale_v` is 2 over the barrier gate's slope (t's e-fold scale while closed). The WKB
  estimates stand only when the potential has no top between the dots at the operating point.
- **Past the merge.** The shape (positions, radii, capacitances, lever arms, the lead actions) is still
  held at the grid's end and flagged `dots_merged` (ADR 0030); only the barrier top continues, because it
  is the potential at a point of the held dots' cut, which the gates keep moving linearly.
- **The consumer surface.** `State.interdot_barrier_top_ev`.
- **Release.** 0.5.0 (breaking). Model versions: `charge.barrier` 1, `charge.shape` 3, `charge.splitting` 2,
  `charge.hybrid` 3, `readout.rates` 2, `compile.device` 7, `device` 8, `device.session` 9.

## Consequences
- Easier:
  - On the reference device f goes from 0.2 to 0.9 over 64 mV of B_M (was about 30) and over 126 mV of
    the plunger diagonal (was 15); each plunger moves the top a quarter as much as B_M; t changes by a
    factor e every 23 mV of B_M (was 9). Within ±0.15 V of the programmed plungers nothing merges at
    the programmed B_M (before, about 40 % of the window was one dot).
  - With B_M at 5 V the top is about 260 meV below the level at every plunger setting: one dot at every
    occupancy, while the plungers still remove electrons.
  - Below the compiled plungers, opening B_M raises t from 2 µeV to its ceiling and f from 0.16 to 0.92
    between 0.2 and 0.4 V (#69).
  - t, Γ_m and g agree everywhere: Γ_m/t² = f_a/t_max² exactly, and t reaches ħω₀/π where g reaches 1.
  - B_M's cross-capacitance and its action on the barrier now agree about what B_M does past the merge.
- Harder:
  - **Two more placeholders.** The screening factors are set by hand to meet the owner's requirement, not
    from a published lever arm; a `literature` row would need a digitised f(V_B) or t(V_B) curve.
  - **Only the interdot barrier is screened.** The lead and sensor barriers keep the bare potential's
    slopes (#48): a plunger still moves a lead rate by a factor e every 5 mV while it moves t by a factor
    e every 90 mV. Extending the two factors to those barriers is #48's work.
  - **A larger ceiling in use.** t reaches ħω₀/π (1.3 meV on the reference device) along B_M, where it held
    at about half that before; the hybridized window's truncation there is the ±4 window's.
  - **t at the operating point moved.** Kemble's T_ref is 0.0018 against the WKB action's 0.0012, so the
    reference device's compiled t is 57 µeV, not 46; the WKB estimate stands only without a top.
  - **A continued top is a linear extrapolation.** Far past the merge ΔU is the end slope times the
    distance, a few hundred meV below the level at 5 V; T is 1 to double precision there, so nothing
    depends on the value, but it is not a potential anyone computed.
- Must now:
  - Device files from 0.4.0 do not load (format 6); recompile them from their seed. The stored gallery
    device is re-recorded.
  - A consumer that read the honeycomb's merged region at B_M +30 mV finds it at about +110 mV; the
    protocols' B_M offsets moved accordingly.
  - A consumer that wants the bare potential sets both screening factors to 1.

## Alternatives considered
- **Scale only the non-barrier gates' slopes (the first of #70's suggestions alone).** Rejected as
  insufficient: along B_M the transition would stay at about 30 mV, and the hold past the merge would
  stand.
- **Widen Kemble's step by hand (a parameter on ħω_b).** Rejected. ħω_b is the fitted curvature of a
  computed top, and widening the step while the gates move the top at full speed misstates which quantity
  is uncertain. The screening factors scale how far the gates move the top, which is what the solver will
  change.
- **Hold B_M past the merge and expose the held range to the consumer (the third suggestion's second
  half).** Rejected by the physics: a barrier gate that far open removes the barrier at every occupancy.
  `State.dots_merged` still marks the held shape.
- **Continue every shape quantity past the merge.** Rejected (ADR 0030): a two-dot shape beyond the point
  where the potential forms one dot is meaningless; the top is the one quantity that is a potential at a
  fixed point.
- **Keep t anchored on the WKB action and let only its slopes follow Kemble.** Rejected: two transmissions
  for one barrier, differing by a factor 1.5 at the operating point, with T = ½ at two different tops.
- **Screen the lead and sensor barriers in the same change.** Deferred to #48: it changes the latching,
  pinch-off and sensor protocols for reasons the owner has not yet asked about.
