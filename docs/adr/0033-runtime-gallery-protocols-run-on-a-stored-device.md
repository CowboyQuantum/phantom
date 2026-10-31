# 0033 — Runtime gallery protocols run on a stored device: compiled once, committed, checked against a fresh compile

**Status:** Accepted
**Date:** 2026-10-31
**Step:** 8 follow-up (found on #66's PR, #67)

## Context
The Windows gallery comparison (ADR 0010, 0014) failed `session-closing-barrier` on some runners: by
1.1·10⁻⁹ on #63 and by 1.1–2.2·10⁻⁹ on every run of #67. Other runners of the same commits agreed to
7·10⁻¹³. The Windows arrays (kept since this change) located the cause:

- **Compile differs between platforms.** Compile's numbers differ between Windows and Ubuntu by about
  10⁻¹². `compiled-device-draws` shows 10⁻¹²–10⁻¹¹, and before #58 did too. Even Γ_R along a plain
  `State` sweep differs by 2.4·10⁻¹²: about 600 ulps of its exponent.
- **This protocol amplifies it about a thousand-fold.** A latched dot sampled at the middle of its charge
  transition turns the energies' 10⁻¹² into 10⁻¹⁰ of charge: the charge moves by an electron over 4kT.
  The sensor reads that charge on its steepest flank, about another hundredfold.
- **No single runtime step is the source.** Perturbing each runtime step by its own last digits moves
  the protocol by at most 10⁻¹¹. The steps tried were the CI model, the shape grid, the rate matrix, the
  propagator, the sensor's chain, its occupation and its back-action.

The owner chose to remove the compile noise from the runtime protocols, not to desensitise the protocol
or hunt the noise in compile (2026-10-31).

## Decision
The runtime gallery protocols run on a device compiled once and committed. Specifically:

- **Storage.** `gallery.protocols.step7._compiled` loads `docs/gallery/devices/<key>.json`, with key a hash
  of the instance fields that set the device. One device serves all 17 runtime protocols today.
  `scripts/gallery.py --record-devices` compiles and stores them.
- **Stale or missing devices.** A stored device whose provenance (format and compile model versions, the
  model versions of every module compile draws on, now including `charge.shape`, from
  `compile_provenance()`) differs from the running code stops the gallery with that instruction. So does
  a missing one.
- **A P0 test.** `tests/test_gallery_devices.py` requires every stored device to belong to a protocol, to
  be current, and to equal a fresh compile on the test's platform within 1e-9 in every number. That
  catches a compile change made without a version bump.
- **Compile still tested.** The protocols that test compile itself compile fresh: `compiled-device-
  roundtrip` and `compiled-device-draws`. Their own comparisons stay within tolerance (10⁻¹¹).
- **Windows arrays kept.** CI uploads the Windows arrays even when the comparison fails, so a difference
  can be located.

## Consequences
- Easier:
  - Every platform runs the runtime protocols on the same compiled numbers, so the comparison measures
    the runtime alone. The amplifying protocols are no longer at the mercy of the runner's CPU.
  - A runtime protocol's compile phase loads JSON instead of compiling.
- Harder:
  - **Re-recording.** A change to compile, or to any model version compile records, needs
    `--record-devices` before the gallery runs. The gallery then refuses to run with a stale device.
  - **A 64 kB file in git** per distinct compile instance.
  - **A different reference.** The gallery's arrays are now those of the stored device, compiled on the
    machine that recorded it. They are not those of each platform's own compile, which differ by about
    10⁻¹².
- Must now:
  - The session protocol (`CLAUDE.md`, End) re-runs `scripts/gallery.py --record-devices` when compile
    changes, beside `--record`.

## Alternatives considered
- **Find and remove the compile noise.** Deferred as open-ended. Candidates are the closed-form basis's
  transcendental functions, pairwise sums of differing order, and finite-difference slopes, each a few
  ulps amplified.
- **Desensitise the protocol (report the charge, read the sensor off its flank).** Rejected by the owner.
  The noise would stay and could reach another protocol.
- **A looser cross-platform tolerance for session protocols.** Rejected: it breaks hard rule 3.
