"""Load protocols, run them with timing and counters, hash arrays, check sign-offs and budgets, render."""

from __future__ import annotations

import dataclasses
import hashlib
import importlib
import json
import pathlib
import platform
import subprocess
import time
from dataclasses import dataclass, field
from typing import Any, Self

import numpy as np
import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROTOCOLS_DIR = ROOT / "docs/gallery/protocols"
SIGNOFF_DIR = ROOT / "docs/gallery/signoff"
MEASURED = ROOT / "docs/gallery/measured.json"
"""Last recorded protocol costs and model benchmarks (ADR 0011), written by `scripts/gallery.py --record`."""
CI_TIME_FACTOR = 3.0
"""Wall time may exceed its budget by this factor on shared CI runners (ADR 0009); `--strict` uses 1."""
CROSS_PLATFORM_RTOL = 1e-12
"""Arrays from another platform must agree with the reference platform's to this relative tolerance,
measured per array as max|a − b| / max|a| (ADR 0010)."""
REFERENCE_PLATFORM = "ubuntu-latest"


@dataclass(frozen=True)
class Protocol:
    id: str
    name: str
    description: str
    layers: list[str]
    runner: str
    """`module:function` under `gallery.protocols`, called as function(ctx, protocol) -> dict of arrays."""
    instance: dict[str, Any]
    sweep: list[dict[str, Any]]
    readout: dict[str, Any]
    render: dict[str, Any]
    budget: dict[str, Any]
    must_show: list[str] = field(default_factory=list)
    step: int | None = None

    @classmethod
    def load(cls, path: pathlib.Path) -> Protocol:
        d = yaml.safe_load(path.read_text(encoding="utf-8"))
        if d.get("id") != path.stem:
            raise ValueError(f"{path.name}: id must equal the file name")
        for key in ("budget", "runner", "sweep", "readout", "render"):
            if key not in d:
                raise ValueError(f"{path.name}: missing {key!r}")
        for phase, keys in (("compile", ("wall_s", "solver_calls")), ("run", ("wall_ms_per_point", "points", "model_evals_per_point", "solver_calls"))):
            missing = [k for k in keys if k not in d["budget"].get(phase, {})]
            if missing:
                raise ValueError(f"{path.name}: budget.{phase} missing {missing}")
        if d["budget"]["run"]["solver_calls"] != 0:
            raise ValueError(f"{path.name}: run.solver_calls must be 0 (rule 4)")
        names = {f.name for f in dataclasses.fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in names})


class Context:
    """Timing and counters a protocol function reports into (ADR 0009)."""

    def __init__(self) -> None:
        self.counts: dict[str, int] = {"model_evals": 0, "solver_calls": 0, "compile_solver_calls": 0}
        self.compile_s = 0.0
        self.run_s = 0.0
        self._phase: str | None = None
        self._t0 = 0.0
        self.points = 0

    def count(self, **kwargs: int) -> None:
        for k, n in kwargs.items():
            if k not in self.counts:
                raise KeyError(f"unknown counter {k!r}")
            self.counts[k] += int(n)

    def compile(self) -> _Phase:
        return _Phase(self, "compile")

    def run(self, points: int) -> _Phase:
        self.points = int(points)
        return _Phase(self, "run")


class _Phase:
    def __init__(self, ctx: Context, name: str) -> None:
        self.ctx, self.name = ctx, name

    def __enter__(self) -> Self:
        self.t0 = time.perf_counter()
        return self

    def __exit__(self, *exc: object) -> None:
        dt = time.perf_counter() - self.t0
        if self.name == "compile":
            self.ctx.compile_s += dt
        else:
            self.ctx.run_s += dt


@dataclass
class Result:
    protocol: Protocol
    arrays: dict[str, np.ndarray]
    sha256: str
    compile_s: float
    run_s: float
    points: int
    counts: dict[str, int]
    signoff: dict[str, Any] | None
    signoff_status: str  # signed | unsigned | MISMATCH | unverifiable-here
    budget_ok: bool
    budget_notes: list[str]
    reference_sha256: str | None = None
    reference_max_rel: float | None = None
    reference_ok: bool | None = None  # None when no reference was given


def compare_arrays(a: dict[str, np.ndarray], b: dict[str, np.ndarray], rtol: float = CROSS_PLATFORM_RTOL) -> tuple[bool, float]:
    """Whether every array in `a` agrees with its counterpart in `b` within `rtol` (ADR 0010), and the
    worst relative difference seen. Integer and boolean arrays must be equal exactly."""
    if set(a) != set(b):
        return False, float("inf")
    worst = 0.0
    for key, value in a.items():
        x, y = np.asarray(value), np.asarray(b[key])
        if x.shape != y.shape or x.dtype != y.dtype:
            return False, float("inf")
        if x.dtype.kind in "fc":
            if not np.array_equal(np.isfinite(x), np.isfinite(y)):
                return False, float("inf")
            fx, fy = x[np.isfinite(x)], y[np.isfinite(y)]
            if fx.size:
                scale = max(float(np.max(np.abs(fx))), np.finfo(float).tiny)
                worst = max(worst, float(np.max(np.abs(fx - fy))) / scale)
        elif not np.array_equal(x, y):
            return False, float("inf")
    return worst <= rtol, worst


def load_reference(reference_dir: pathlib.Path | None, pid: str) -> dict[str, np.ndarray] | None:
    if reference_dir is None:
        return None
    path = reference_dir / f"{pid}.npz"
    if not path.is_file():
        return None
    with np.load(path) as z:
        return {k: z[k] for k in z.files}


def canonical_hash(arrays: dict[str, np.ndarray]) -> str:
    """sha256 over the arrays in sorted key order: key, dtype, shape, raw bytes (ADR 0008)."""
    h = hashlib.sha256()
    for key in sorted(arrays):
        a = np.ascontiguousarray(arrays[key])
        h.update(key.encode())
        h.update(a.dtype.str.encode())
        h.update(str(a.shape).encode())
        h.update(a.tobytes())
    return h.hexdigest()


def load_protocols(only: list[str] | None = None) -> list[Protocol]:
    paths = sorted(PROTOCOLS_DIR.glob("*.yaml"))
    protocols = [Protocol.load(p) for p in paths]
    if only:
        protocols = [p for p in protocols if p.id in only]
    return protocols


def load_signoff(pid: str) -> dict[str, Any] | None:
    path = SIGNOFF_DIR / f"{pid}.yaml"
    return yaml.safe_load(path.read_text(encoding="utf-8")) if path.is_file() else None


def resolve_runner(spec: str):
    module_name, func_name = spec.split(":")
    module = importlib.import_module(f"gallery.protocols.{module_name}")
    return getattr(module, func_name)


def check_budget(p: Protocol, ctx: Context, time_factor: float) -> tuple[bool, list[str]]:
    notes: list[str] = []
    b = p.budget
    ok = True
    if ctx.points != b["run"]["points"]:
        ok = False
        notes.append(f"points {ctx.points} != declared {b['run']['points']}")
    declared_evals = b["run"]["model_evals_per_point"] * b["run"]["points"]
    if ctx.counts["model_evals"] != declared_evals:
        ok = False
        notes.append(f"model_evals {ctx.counts['model_evals']} != declared {declared_evals}")
    if ctx.counts["solver_calls"] != 0:
        ok = False
        notes.append(f"run-phase solver_calls {ctx.counts['solver_calls']} != 0")
    if ctx.counts["compile_solver_calls"] != b["compile"]["solver_calls"]:
        ok = False
        notes.append(f"compile solver_calls {ctx.counts['compile_solver_calls']} != declared {b['compile']['solver_calls']}")
    per_point_ms = 1e3 * ctx.run_s / max(ctx.points, 1)
    if per_point_ms > time_factor * b["run"]["wall_ms_per_point"]:
        ok = False
        notes.append(f"run {per_point_ms:.4f} ms/pt > {time_factor:g}x budget {b['run']['wall_ms_per_point']}")
    if ctx.compile_s > time_factor * b["compile"]["wall_s"]:
        ok = False
        notes.append(f"compile {ctx.compile_s:.3f} s > {time_factor:g}x budget {b['compile']['wall_s']}")
    return ok, notes


def run_protocol(
    p: Protocol,
    time_factor: float = CI_TIME_FACTOR,
    reference_dir: pathlib.Path | None = None,
    is_reference_platform: bool = True,
    rtol: float = CROSS_PLATFORM_RTOL,
) -> Result:
    """Run `p`. On the reference platform the sign-off is checked by exact hash. Elsewhere, with
    `reference_dir` holding the reference platform's `.npz` files, the arrays are compared within `rtol`
    and the sign-off is judged by the reference hash; without a reference, a differing hash is
    `unverifiable-here` (ADR 0010)."""
    ctx = Context()
    arrays = resolve_runner(p.runner)(ctx, p)
    arrays = {k: np.asarray(v) for k, v in arrays.items()}
    sha = canonical_hash(arrays)
    signoff = load_signoff(p.id)
    reference = load_reference(reference_dir, p.id)
    ref_sha = canonical_hash(reference) if reference is not None else None
    ref_ok, ref_rel = (compare_arrays(arrays, reference, rtol) if reference is not None else (None, None))
    judged_sha = ref_sha if reference is not None else sha
    if signoff is None:
        status = "unsigned"
    elif signoff.get("array_sha256") == judged_sha and (ref_ok is None or ref_ok):
        status = "signed"
    elif not is_reference_platform and reference is None:
        status = "unverifiable-here"
    else:
        status = "MISMATCH"
    ok, notes = check_budget(p, ctx, time_factor)
    return Result(
        p, arrays, sha, ctx.compile_s, ctx.run_s, ctx.points, dict(ctx.counts), signoff, status, ok, notes,
        reference_sha256=ref_sha, reference_max_rel=ref_rel, reference_ok=ref_ok,
    )


def git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, cwd=ROOT, check=False
        ).stdout.strip()
    except OSError:
        return "unknown"


def ledger_row(r: Result) -> dict[str, Any]:
    b = r.protocol.budget
    return {
        "id": r.protocol.id,
        "step": r.protocol.step,
        "points": r.points,
        "compile_s": round(r.compile_s, 4),
        "compile_budget_s": b["compile"]["wall_s"],
        "run_ms_per_point": round(1e3 * r.run_s / max(r.points, 1), 5),
        "run_budget_ms_per_point": b["run"]["wall_ms_per_point"],
        "model_evals_per_point": r.counts["model_evals"] / max(r.points, 1),
        "declared_model_evals_per_point": b["run"]["model_evals_per_point"],
        "solver_calls": r.counts["solver_calls"],
        "budget_ok": r.budget_ok,
        "budget_notes": r.budget_notes,
        "sha256": r.sha256,
        "signoff": r.signoff_status,
        "reference_sha256": r.reference_sha256,
        "reference_max_rel": r.reference_max_rel,
        "reference_ok": r.reference_ok,
    }


def write_outputs(results: list[Result], out: pathlib.Path, render: bool) -> None:
    out.mkdir(parents=True, exist_ok=True)
    rows = [ledger_row(r) for r in results]
    meta = {"commit": git_commit(), "platform": platform.platform(), "python": platform.python_version(), "numpy": np.__version__}
    (out / "budget.json").write_text(json.dumps({"meta": meta, "protocols": rows}, indent=1), encoding="utf-8")
    (out / "hashes.json").write_text(json.dumps({r.protocol.id: r.sha256 for r in results}, indent=1, sort_keys=True), encoding="utf-8")
    for r in results:
        np.savez(out / f"{r.protocol.id}.npz", **r.arrays)
    if render:
        from gallery.render import render_index, render_protocol

        for r in results:
            render_protocol(r, out / f"{r.protocol.id}.png")
        render_index(results, rows, meta, out / "index.html")


def record_measured(results: list[Result]) -> pathlib.Path:
    """Write `docs/gallery/measured.json`: per-protocol measured cost from this run and per-model
    micro-benchmarks, with the commit, platform and date. Machine-dependent by nature; the reference
    platform's numbers are in the CI artifact."""
    import datetime

    from gallery.bench import run as bench

    data = {
        "meta": {
            "commit": git_commit(),
            "platform": platform.platform(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "date": datetime.date.today().isoformat(),
        },
        "protocols": {
            r.protocol.id: {
                "points": r.points,
                "compile_s": round(r.compile_s, 4),
                "run_ms_per_point": round(1e3 * r.run_s / max(r.points, 1), 5),
                "model_evals_per_point": r.counts["model_evals"] / max(r.points, 1),
            }
            for r in results
        },
        "models": bench(),
    }
    MEASURED.write_text(json.dumps(data, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return MEASURED
