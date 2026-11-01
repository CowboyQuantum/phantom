# 0035 — The dots grow with their occupancy: classical puddles whose radius scales as N^(1/3), their capacitances following the gate-induced charge and the interdot barrier seen from their edges

**Status:** Accepted
**Date:** 2026-11-01
**Step:** 8 follow-up (issue #72; after ADR 0034)

## Context
Reviewing 0.5.0 (#70), the owner judged the honeycomb-to-superdot transition still "practically binary"
across a plunger window, and asked for the physics of the dots' size: as B_M opens, the triple-point
width should grow first at the highest occupancies, where the dots' wave functions reach furthest, and
move toward the (1,1) cell (2026-11-01). Measured on the reference device:

- **The front existed, with one mechanism.** Walking both plungers from the programmed point to the
  (2,3) cell, the superdot (f > 0.9) began at (25,29) with B_M at 0.30 V, at (19,19) at 0.50 V and
  everywhere at 0.80 V. The only thing driving it was the plungers' cross-coupling to the barrier top,
  about 0.5 meV per electron after ADR 0034's screening placeholder.
- **The dots did not grow.** Along the same walk the dot radius stayed at 9.8 nm and the addition
  energy at 9.4 meV from (2,3) to (25,31). The closed-form dot is a harmonic minimum with its
  one-electron radius (ADR 0013), and the constant-interaction model adds electrons without changing
  anything geometric.
- **The bare profile is shallow.** On the cut between the dots at the operating point the bare
  potential sits 4.1 meV (dot 1) and 25.1 meV (dot 2) below the tunnelling level; the forbidden stretch
  is 30 nm wide; an electron at the level turns around 15 nm from dot 1's centre and 37 nm from dot 2's.
- **Published.** A two-dimensional puddle of charge in a parabolic well has the semi-elliptic density
  of Shikin, Nazin, Heitmann and Demel (1991), which gives N ∝ R³, so R ∝ N^(1/3) and e²/C ∝ N^(−1/3).
  In an etched graphene double dot the interdot coupling rose monotonically from about 70 to 300 µeV as
  electrons were added to one dot (Sci. Rep. 3, 3175); in laterally coupled vertical dots it depends on
  the orbital the added electron occupies (PRL 93, 066806). No general law of coupling against electron
  number for gate-defined dots was found.

Four forces shaped the design:

- **The owner's requirement.** At a fixed barrier gate the honeycomb should read differently at
  different occupancies: a grid near (1,1), hexagons, then a superdot, and the front should move with
  the barrier gate.
- **Honest mechanism.** In a WKB picture at the Fermi level the "overlap of the wave functions" is the
  transmission of the barrier above that level. A self-consistent puddle raises the potential inside
  the dot to the Fermi level and sees the barrier from its edge; its growth with N and its repulsion of
  the saddle are the solver's business (plan §8 #1). A closed form can take the geometric part, the
  edge moving along the bare profile, and must say what it leaves out.
- **Bit-exact anchor.** At the operating point nothing may change: the compiled capacitances, t, the
  rates, the splitting and every test there.
- **No double counting of the gates.** ADR 0034's screened slopes already carry how the gates move the
  top. The occupancy term must depend on the occupancy alone.

## Decision
Each dot is a classical puddle whose radius scales as N^(1/3) from the compiled dot; its capacitances
follow its gate-induced charge, and the interdot barrier is seen from the puddles' edges. Specifically:

- **Scale.** At a gate point the geometric model (compiled lever arms and capacitances, the shape grid)
  gives each dot's gate-induced charge n_i. Its puddle's scale is s_i = √(max(n_i, ½)/N_ref,i), with
  N_ref,i its gate-induced charge at the operating point: C ∝ R ∝ N^(1/3) and N = Cφ/e give C ∝ √n.
  `charge.puddle.PuddleGrowth`.
- **Capacitances.** Each dot's self-capacitance (its total minus the mutual one, ADR 0013) is multiplied
  by s_i and the mutual capacitance by s₁s₂, written as changes from the compiled matrix so that s = 1
  returns it bit for bit. The induced charge is that of the scaled capacitances, s·n ∝ n^(3/2), the
  classical N ∝ φ^(3/2). `ConstantInteractionParams.puddle`; `local`, `induced_charge`, the window,
  the energies, the lever arms and `addition_energies_at_ev` all follow; `occupancy_scale(v)` exposes s.
- **Radius.** Each puddle's radius at the operating point is the bare potential's turning point at the
  interdot tunnelling level on the compiled cut toward the other dot: the puddle of N_ref electrons
  fills the bare allowed region. Elsewhere R_i = s_i·R_ref,i. `State.puddle_radii_nm`; `dot_radii_nm`
  keeps the shape grid's harmonic radius.
- **The barrier seen from the edge.** Compile stores the bare potential above the level along the cut
  at the operating point. The barrier's top is raised by the mean over the two dots of
  u(R_ref,i) − u(R_i): zero at the anchor, positive for a smaller puddle (it sees the profile from
  further in, where it is lower, a rigid shift of the profile outside the puddle as a flat
  self-consistent potential meets the level at its edge), negative for a larger one, −ΔU_ref once both
  reach the saddle. `InterdotBarrier.puddle` and a `scale` argument on `barrier_top_ev`, `transmission`
  and `log_odds_change`; `InterdotSplitting`, `InterdotCoupling`, `GateRates` and `RateEquation` pass
  the scale through, so t, Γ_m and g share it (ADR 0034).
- **Compile.** `CompileParams.puddles` (True). After the device is built, N_ref is read through the
  runtime CI model at the operating point, so the runtime finds exactly 1 there; the cut, the turning
  points (`puddle_turning_points_nm`, the saddle when the cut never rises above the level) and N_ref
  are stored in `Device.puddle` (`DEVICE_FORMAT_VERSION` 7).
- **Release.** 0.6.0 (breaking). Model versions: `charge.puddle` 1, `charge.ci` 3, `charge.barrier` 2,
  `charge.splitting` 3, `charge.hybrid` 4, `readout.rates` 3, `readout.master_equation` 4,
  `compile.device` 8, `device` 9, `device.session` 10. Reference entry
  `shikin1991-parabolic-dot-classical-density` (`proposed`). Protocol `honeycomb-vs-occupancy`.

## Consequences
- Easier:
  - On the reference device, at the programmed gates, from (1,1) to (28,35): the puddle radii grow from
    5/12 nm to 17/40 nm, the addition energies fall from 32/28 meV to 6/6 meV, the geometric splitting
    f_C grows from 0.06 to 0.16, the barrier top falls from 42 meV above the level to below it, and t
    rises from 10⁻¹⁰ µeV to its ceiling.
  - The superdot front moves with the barrier gate across occupancy, about 25 electrons per 0.5 V of
    B_M: (26,29) at 0.3 V, (19,19) at 0.5 V, (10,7) at 0.8 V, the first electrons at 1.0 V. A plunger
    map at one B_M shows a grid near (1,1), hexagons in the middle and a superdot at high occupancy
    (`honeycomb-vs-occupancy`).
  - The classical N ∝ φ^(3/2): more electrons per volt of plunger at high occupancy than at low.
- Harder:
  - **The anchor's inconsistency.** The compiled radius is the one-electron harmonic extent (ADR 0013),
    read here as the puddle of N_ref ≈ 23 electrons; scaling it down to one electron gives a 5 nm
    puddle and a 32 meV addition energy, large for a Si/SiGe dot. Anchoring at N = 1 instead would
    change the operating point (about 110 electrons at the same gates) and every protocol.
  - **Classical and geometric only.** The N^(1/3) law ignores the kinetic energy (shells) and the
    puddles' own screening of the gates and of each other; the edge shift uses the bare profile, not a
    self-consistent one; the compiled cut is held where the shape grid moves the dots. All
    `calibration: none`.
  - **The mutual capacitance's share.** The scaled induced charge is s·n only up to the mutual term,
    which scales as s₁s₂ (within 2 % on the reference device).
  - **Two placeholders act along the plungers.** The gates' screened cross-coupling (ADR 0034) and the
    puddles' growth both raise the barrier toward low occupancy; the transition along the plunger
    diagonal is 106 mV wide (was 126), along B_M 62 mV (was 64).
  - **Cost.** The scale needs the shape grid at every evaluation of t, the rates and the splitting, so
    the CI model keeps its last gate point's geometry and the barrier its last top (as the eigenstates
    were already kept): a session sample now costs about 30 % less than before this change rather than
    40 % more. The first Windows CI run of the change, before that, put three session protocols at
    3.3× their budget on a slow runner.
- Must now:
  - Device files from 0.5.0 do not load (format 7); recompile them from their seed. The stored gallery
    device is re-recorded.
  - Tests that compare the runtime capacitances with the compile tier along the barrier gates turn the
    puddles off (`puddle=None`), as they turn the splitting off.
  - A consumer that wants fixed dots sets `CompileParams.puddles=False`.

## Alternatives considered
- **Shorten the WKB action by the puddles' growth.** Rejected: the bare WKB at the fixed level already
  starts its forbidden region where the bare potential crosses the level, so taking it from the puddle
  edge as well would count the growth twice; the edge shift on the bare profile is the one consistent
  statement.
- **A per-electron shift of the barrier top as a parameter.** Rejected in favour of the edge on the
  compiled cut, which has no free number and gives each dot its own shift from its own depth.
- **Anchor the puddle at one electron (R(1) = ℓ).** Rejected for now: physically the cleaner anchor,
  but it moves the operating point to about 110 electrons and every protocol with it; left for the
  solver, which will give the radius at every N.
- **An N-dependent capacitance per window state.** Rejected: the quadratic-form energies with a
  per-state matrix cannot anchor both the line positions and the addition energy without an extra
  linear term, and every consumer of `local(v)` would need the state. Scaling at the gate point's
  continuous occupancy keeps one entry point.
- **The Hartree repulsion of the puddles on the saddle.** Deferred to the solver: it is the screening of
  #48, and its net with the gates' action is the number the screening placeholders stand in for.
