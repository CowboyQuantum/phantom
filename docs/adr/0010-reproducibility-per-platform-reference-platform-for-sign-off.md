# 0010 — Bit-for-bit reproducibility holds per platform; across platforms within 1e-12; ubuntu-latest is the reference

**Status:** Accepted (amends 0003 rule 3 and 0008)
**Date:** 2026-10-25
**Step:** 1

## Context
Rule 3 (ADR 0003) promised identical output, bit for bit, on every platform CI runs. The first
cross-platform gallery run (PR #19) showed that ten of eleven protocols hash differently on
ubuntu-latest and windows-latest, including two that evaluate nothing but `exp`, `tanh` and
`logaddexp` of a gate voltage. numpy's elementary functions are dispatched to different SIMD and
C-library implementations on the two builds and differ in the last unit of precision. No model
written in numpy can promise bit identity across operating systems; the promise as written was false.

## Decision
- **Per platform, bit for bit.** For the same parameters, seed, `model_version`, operating system and
  numpy build, every output is identical byte for byte. This is what a consumer relies on: determinism
  by seed. CI asserts it by the Windows job and the Ubuntu job each reproducing their own fixtures.
- **Across platforms, within tolerance.** Arrays produced on different platforms agree to a relative
  tolerance of 1e-12, measured per array as the largest absolute difference divided by the largest
  magnitude in the array. The tolerance is a constant in `gallery/runner.py` and is stated here.
- **One reference platform for sign-off.** Sign-off records (ADR 0008) carry the array hash computed
  on `ubuntu-latest` in CI. That job is the reference. Every other platform's gallery job downloads the
  reference arrays and passes only if its own arrays agree within the tolerance; it reports a protocol
  as signed when the reference hash matches the record and its arrays agree with the reference.
- **Locally**, on a machine that is not the reference platform, a hash mismatch against a sign-off is
  reported as "unverifiable here" rather than as a failure; the verdict belongs to CI.

## Consequences
- Easier: the contract is true and checkable; sign-off has one unambiguous hash.
- Harder: a sign-off can only be verified exactly by the reference job; a change that moves an array by
  less than 1e-12 on another platform is invisible there, which is the intended meaning of the tolerance.
- Must now: CLAUDE.md rule 3 reads as above; the `gallery` CI jobs are `gallery` (reference) and
  `gallery-windows` (compared to it); `tests/test_qdt_sim_migration.py` keeps its 1e-12 tolerances, which
  were already the honest ones.

## Alternatives considered
- Disable SIMD dispatch in numpy on CI — rejected: the C libraries' `exp` still differ, and consumers
  would not run with that setting.
- Hash arrays rounded to a fixed number of digits — rejected: a value on a rounding boundary flips the
  hash at random; a tolerance comparison has no such edge.
- Keep the claim and skip the check — rejected: an untestable promise is worse than a weaker one.
