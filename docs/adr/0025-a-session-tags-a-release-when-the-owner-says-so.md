# 0025 — A session tags a release when the owner says so, on a green `main` with the changelog entry; the owner can always tag

**Status:** Accepted
**Date:** 2026-10-27
**Step:** 8

## Context
ADR 0024 decided that the owner tags releases, and issue #40's checklist said the same. After PR #44 the
owner asked that a session be able to tag when told to (#45), so that a release does not wait on the owner
running git. A tag is what consumers pin (qdt pins phantom by tag, ADR 0003), so the tagged commit must be
one that CI has passed and whose changelog names the version. Trying it on 2026-10-27, the session's
branch pushes succeeded but its push of `v0.2.0` was refused with HTTP 403, and its GitHub tools have no
call that creates a tag: whether a session can tag depends on the credentials of its environment.

## Decision
- **Release bullet of ADR 0024, amended:** a session tags a release only when the owner says so in the
  conversation. The tag is annotated, named `vX.Y.Z` after `pyproject.toml`'s version, and points at a commit
  on `main` whose CI run is green and whose `CHANGELOG.md` has the `X.Y.Z` entry. The session confirms
  all three before tagging and says on the step's issue which commit it tagged.
- The owner can always tag directly, on GitHub's release page or with git.
- A session whose environment refuses the tag push says so and leaves the tag to the owner. It never
  works around the refusal through another channel.
- How often phantom releases and qdt moves its pin stays the owner's open decision (Step 8).

## Consequences
- Easier: a release can be cut in the same conversation that finishes it.
- Harder: a session has three preconditions to check before a tag, and a tag pushed by mistake has to be
  deleted by hand on GitHub, where consumers may already have fetched it.
- Must now: `CLAUDE.md` carries the rule; the documents that said "the owner tags releases"
  (`CHANGELOG.md`, `docs/CONTRACT.md`, the Step 8 file) say this instead. `v0.2.0` itself is tagged by the
  owner, because this session's environment refused the push.

## Alternatives considered
- Keep tagging the owner's alone (ADR 0024 as written): rejected at the owner's request.
- A CI workflow that tags when the version in `pyproject.toml` changes on `main`: rejected for now. It tags
  without anyone saying so, and the owner wants the decision to release kept explicit.
- Let a session tag whenever a step closes: rejected, since cadence is the owner's open decision.
