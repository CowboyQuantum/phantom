# phantom

A stand-in spin-qubit device: the physics behind gate-defined quantum dots, from gate layout to spin
readout, with the disorder and fabrication defects that make real devices hard to tune. `phantom` is a
plain Python library (numpy/scipy); it knows nothing about instruments or measurement software. qdt uses it
as a simulated backend; agents and products use it to generate scenarios and grade tune-ups.

**Documentation:** https://cowboyquantum.github.io/phantom/ — the development plan, the step files, the ADRs,
[how the models fit together](https://cowboyquantum.github.io/phantom/MODELS/) with their runtimes, and the
[rendered gallery](https://cowboyquantum.github.io/phantom/gallery/rendered/index.html) of measurement protocols.
Every pull request publishes its own preview of the same site at `/pr-<n>/` (ADR 0012).

- Humans: start at the documentation site, or at `docs/DEVELOPMENT_PLAN.md` (roadmap §7).
- Agents: start at `CLAUDE.md`, then `docs/STATUS.md`.

```
uv sync --all-extras && uv run pytest && bash scripts/check_all.sh
```
