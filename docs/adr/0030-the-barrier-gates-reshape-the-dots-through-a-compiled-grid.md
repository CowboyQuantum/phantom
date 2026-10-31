# 0030 — The barrier gates reshape the dots through a compiled grid: the double dot re-derived along each barrier gate, its changes interpolated and added at runtime, held where the dots merge

**Status:** Accepted
**Date:** 2026-10-31
**Step:** 8 follow-up (issue #59, the second of #58–#60; after ADR 0029)

## Context
After ADR 0029 the interdot coupling followed the gates, but the dots themselves did not. Compile placed
them once at the operating point. Their positions, sizes, lever arms, the disorder offsets at them, the
capacitance matrix and the actions of the barriers around them then held at every gate point. A barrier
gate only moved the dots' levels through a fixed cross-capacitance, and through the linear slopes of
the rates and of t. The owner asked for the barrier gates to change the dots' "potential
locations/shape" (2026-10-31), and decided on a grid compiled along the barrier gates and interpolated
at runtime.

Five forces shaped the design:

- **Compile once, run cheap (ADR 0001).** A derivation of the double dot costs 1.5 ms plus three barrier
  actions, so a runtime point cannot call one. A grid of 17 nodes per gate costs about 0.2 s of compile.
- **A nonlinear potential map.** The closed-form potential is linear in the gates, but the minimum is
  not. On the reference device, B_M at +0.15 V draws the dots 18 nm together, C_m grows by 40 % and
  dot 1's addition energy falls by 1 meV. At +0.18 V the potential forms one dot. B_L and B_R move their
  own dot by several nm and merge it with its lead when opened far enough.
- **Pixel steps.** `find_dots` places a minimum on a 2 nm grid point. Taken directly, positions and lever
  arms jump in steps.
- **Bit-exact anchor.** At the operating point nothing may change: the compiled device and every test
  there stay exact.
- **File size.** The device file must stay small, and the basis is not stored.

## Decision
Compile re-derives the double dot along each barrier gate and records each quantity's change from the
operating point. The runtime interpolates each gate's changes and adds them. Specifically:

- **The grid.**
  - `CompileParams.shape_gates` (B_L, B_M, B_R), `shape_span_v` (±0.2 V) and `shape_points` (17).
  - At each node only that gate moves.
  - Each dot's minimum is refined between grid points by a parabola (`refined_dot`). Its lever arms and
    disorder offset are taken there.
  - The capacitance matrix is the derivation's own formula (`capacitance_matrix_af`) at the refined
    dots.
  - The three barrier actions are taken at the tunnelling levels held fixed (ADR 0027), from the
    refined dots.
  - Each grid runs outward from the operating point until the potential stops forming two dots, or a
    dot jumps by more than 12 nm. In that case it merged and another minimum took its place, and the
    end is marked merged.
  - `charge.shape.ShapeGrid` holds the nodes and the changes. `Device.shape`, `DEVICE_FORMAT_VERSION` 4.
- **The runtime.**
  - Each gate's changes are interpolated with PCHIP at its carried voltage and the gates' changes are
    added.
  - Beyond a grid's end the changes are held. `ShapeGrid.merged` flags a point beyond a merged end.
  - `ConstantInteractionParams.shape` gives C(V) and n_g = C(V)·(α(V)·V − Δo(V))/E_UNIT. This is the
    potential at the moved minimum, by the envelope theorem, and it is exact to first order in the
    plungers, which stay linear.
  - `InterdotCoupling.shape` and `GateRates.shape` take the exponents along the barrier gates from the
    grid (t from half the interdot change). The compiled slopes still serve every other gate.
- **The consumer surface.** `State` gains `dot_positions_nm`, `dot_radii_nm`, `addition_energies_ev` and
  `dots_merged`. `optimum(device).coupling` inverts t along the barrier gate by Brent's method.
- **The hybridized window (amends ADR 0029).**
  - The window centres on rint(n_g), which changes at the middle of every interdot segment
    (n_g = N + ½), the very point a tune-up reads. With t near a meV the ±2 chains of ADR 0029 left a seam
    of 2·10⁻⁵ e there. The first Windows comparison of this change differed by 1.8·10⁻⁵ when last-digit
    rounding put a point on the other side.
  - The statics now use a ±4 window. Its truncation, about (t/E_add)^(2w), stays below 10⁻¹² up to t's
    ceiling, at about 25 µs per point instead of 7.
  - The dynamics use the smallest window whose truncation is below 10⁻¹³ at each sample's t
    (`dynamics_window`): ±2 below t/E_add = 10⁻³, ±3 below 1.5·10⁻², ±4 above.
  - `charge.hybrid` 2, `readout.master_equation` 3.
- **Release.** 0.4.0. Model versions: `charge.shape` 1, `compile.device` 5, `device` 6,
  `device.session` 6, `optimum.targets` 3.

## Consequences
- Easier:
  - A barrier gate moves its dots, changes their size and their addition energies, couples them more
    strongly and reshapes the honeycomb (`dot-shape-vs-barrier`).
  - Opening B_M ends in a flagged merge instead of an unbounded extrapolation.
  - The rates and t follow the moved dots. When a barrier falls below its tunnelling level, its action
    is zero and the rate levels off at f_a/(2 − T₀), about half the attempt rate. Other gates' linear
    slopes still take it towards f_a.
- Harder:
  - **Separable changes.** The gates' changes are added, so cross terms (B_L with B_M on dot 1) are left
    out. Opening all three barriers together pushes the rates towards the attempt rate.
  - **The model's own limits.** The capacitances still come from the disk and sphere formulae of ADR
    0013. C_m grows only to about 2.1 aF before the merge, so a merged double dot does not become one dot
    in the CI model. Points beyond the merge hold the last two-dot shape and are flagged.
  - **Pixel noise.** The compile tier takes a dot's curvature at a grid point, so the radius and the
    capacitances move in steps of about 1 % as the minimum changes pixel. The interpolant runs through
    them, and halfway between nodes it is within 2 % of a direct derivation.
  - **Absolute charge.** The CI model counts the charge absolutely (n_g = Cφ/E ≈ 23 e). A dot that grows
    by 10 % takes about two more electrons at the same potential, so a barrier gate moves the lines
    along the plungers by more than its cross-capacitance. Tests and protocols that follow one
    transition name its electron count.
  - **Cost and size.** Compile takes 0.2 s more. The device file grows from 14 kB to 62 kB. A `State`
    point costs about 25 µs, and a session sample a few interpolations more.
  - **Exponent between nodes.** An exponent that saturates within one 25 mV step is interpolated coarsely
    there, by up to 1 in ln Γ.
  - **Not yet shaped.** The sensor's coupling to the dots (β) and the attempt rates stay at the operating
    point.
- Must now:
  - Device files from 0.3.0 and from earlier 0.4.0 commits do not load (format 4); recompile them from
    their seed.
  - A consumer that cached C_dot, the lever arms or t from the compiled device reads them per point.

## Alternatives considered
- **A tensor-product grid over the three barrier gates.** Rejected for size and compile time: 17³
  derivations and several MB. The cross terms it would capture are second order on the reference device.
- **Linear slopes of each quantity (first derivatives).** Rejected by the owner. They cannot show the
  merge or the saturation, and they drift within tens of mV.
- **Grid-point minima without refinement.** Rejected: positions and lever arms would step by 2 nm and the
  honeycomb would jump.
- **Including the plungers in the grid.** Deferred. Their lever arms are the largest, but at the
  operating point they move the minima less than the barriers do, and the envelope theorem keeps the
  energies right to first order.
- **Crossing into a single-dot model past the merge.** Deferred to the solver (plan §8 #1). The CI model
  has no state for one dot of N₁ + N₂ electrons. Flagging the merge and holding the last two-dot shape
  keeps the model honest about what it does not describe.
