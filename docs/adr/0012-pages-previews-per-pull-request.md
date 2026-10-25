# 0012 — Every pull request publishes a preview of the docs and the rendered gallery; main publishes the root

**Status:** Accepted (amends 0008)
**Date:** 2026-10-25
**Step:** 1 onward

## Context
ADR 0008 publishes the rendered gallery to GitHub Pages from `main`. The owner pointed out that this
defeats the process it serves: review and sign-off are supposed to happen on the branch, before the
merge, and a zip artifact or a hand-rebuilt page is a poor way to look at fourteen pictures and a
budget table. The published site should carry the branch under review as well as the settled state.

## Decision
- A `pages` workflow builds the documentation site plus the rendered gallery for every push to a pull
  request and publishes it under `/pr-<n>/` on the Pages site, served from the `gh-pages` branch. A
  sticky comment on the pull request carries the links (docs, rendered gallery, models and runtimes),
  updated on every push. When the pull request closes, its directory is removed.
- Pushes to `main` publish the root of the site the same way. The root therefore always shows the
  settled state, and `/pr-<n>/` shows what is being reviewed.
- Pages is served from the `gh-pages` branch (deploy-from-branch), not from the Actions deployment,
  because only a branch can hold several directories at once.
- Review and sign-off (ADR 0008) happen on the preview: the owner looks at `/pr-<n>/gallery/rendered/`,
  the discussion and the sign-off records land on the same branch, CI verifies the hashes, then the
  pull request merges.

## Consequences
- Easier: the gallery under review is one click from the pull request; nothing to download or rebuild.
- Harder: the `gh-pages` branch accumulates deploy commits; previews of closed pull requests are removed
  by the workflow but the branch history keeps them.
- Must now: the repository's Pages source is set to the `gh-pages` branch, root folder; the `pages` job
  leaves `ci.yml`; the gallery README points reviewers at the preview.

## Alternatives considered
- Keep publishing from `main` only and review the CI zip — rejected by the owner.
- Separate hosting for previews — rejected: one site, one address scheme, nothing external.
