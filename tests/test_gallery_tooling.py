"""P0 for the gallery tooling (ADR 0008/0009): every protocol file loads and resolves, the canonical
hash is deterministic, and the cheapest protocols run with counts that match their declared budgets."""

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from gallery.runner import canonical_hash, load_protocols, resolve_runner, run_protocol


def test_every_protocol_loads_and_its_runner_resolves():
    protocols = load_protocols()
    assert len(protocols) >= 11
    for p in protocols:
        assert callable(resolve_runner(p.runner))
        assert p.budget["run"]["solver_calls"] == 0


def test_canonical_hash_is_deterministic_and_sensitive():
    a = {"x": np.arange(5.0), "y": np.ones((2, 3))}
    b = {"y": np.ones((2, 3)), "x": np.arange(5.0)}
    assert canonical_hash(a) == canonical_hash(b)
    c = {"x": np.arange(5.0), "y": np.ones((3, 2))}
    assert canonical_hash(a) != canonical_hash(c)
    d = {"x": np.arange(5.0), "y": np.ones((2, 3), dtype=np.float32)}
    assert canonical_hash(a) != canonical_hash(d)


def test_cheap_protocols_run_with_declared_counts():
    for p in load_protocols(["turn-on-single-gate", "pinch-off-barrier", "sensor-coulomb-peak"]):
        r = run_protocol(p, time_factor=1e9)  # counts only; wall time is CI's business
        assert r.points == p.budget["run"]["points"]
        assert r.counts["model_evals"] == p.budget["run"]["points"] * p.budget["run"]["model_evals_per_point"]
        assert r.counts["solver_calls"] == 0
        assert r.budget_ok, r.budget_notes
        assert r.signoff_status == "unsigned"
        assert run_protocol(p, time_factor=1e9).sha256 == r.sha256


def test_compare_arrays_tolerance_and_structure():
    from gallery.runner import compare_arrays

    a = {"x": np.linspace(0, 1, 100), "n": np.arange(4)}
    b = {"x": a["x"] * (1 + 1e-15), "n": np.arange(4)}
    ok, worst = compare_arrays(a, b, 1e-12)
    assert ok and worst < 1e-12
    c = {"x": a["x"] * (1 + 1e-9), "n": np.arange(4)}
    ok, worst = compare_arrays(a, c, 1e-12)
    assert not ok and 1e-10 < worst < 1e-8
    assert compare_arrays(a, {"x": a["x"]}, 1e-12) == (False, float("inf"))
    assert compare_arrays(a, {"x": a["x"], "n": np.arange(4) + 1}, 1e-12) == (False, float("inf"))
