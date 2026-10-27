"""phantom.contract — the consumer contract as a test any consumer can run (ADR 0003, 0024; docs/CONTRACT.md).

    python -m phantom.contract                    # the reference device (the gallery's recipe)
    python -m phantom.contract device.json ...    # and each device file given; exit status 0 or 1

`run(device)` checks the promises of `phantom.api` on a device: compile is deterministic in the seed;
the JSON round trip is text-identical and the reloaded device answers bit for bit; provenance names the
format and every model version; `at` has the documented shapes and finite values and passes the
sensor bias; a fresh session replays the device's noise, a noise seed gives a new realisation of the
same noise and replays it too; the clock only moves forward and every sample sits at its documented
time; with no noise a session in `mean` mode reads `at`'s equilibrium current; a sample's charge is
an integer state in `sampled` mode and ⟨N⟩ in `mean` mode; holding voltages moves the charge there
and a fast read afterwards still shows it; `knowable` holds exactly the fields a user may see; `truth`
carries the seed, the scenario, the dots, the charge lattice and the sensor peak with its wing;
`optimum` returns all five parts, finite; `spin_readout` returns one outcome per shot and takes its
time on the clock. phantom's own CI runs it (`tests/test_contract.py`); qdt's adapter runs it against
the pinned release and its committed device files.
"""

from __future__ import annotations

import itertools
import json
import math
import sys
from collections.abc import Callable
from dataclasses import dataclass, fields

import numpy as np

from phantom import api
from phantom.api import (
    DEVICE_FORMAT_VERSION,
    CompileParams,
    Device,
    compile_device,
    optimum,
    two_dot_layout,
)

REFERENCE_SEED = 2
REFERENCE_OPERATING_POINT_V = {"B_L": 0.15, "P1": 0.6, "B_M": 0.2, "P2": 0.6, "B_R": 0.15, "SG": 0.6, "R_L": 0.7, "R_R": 0.7}
"""The gallery's device: two_dot_layout("SiGe", 80 nm) at this programmed point, seed 2."""

BATCH_RTOL = 1e-12
"""A point evaluated alone and inside a batch agree within this, relative; identical calls agree bit for bit."""

KNOWABLE_FIELDS = ("material", "gate_names", "layout", "temperature_k", "programmed_point_v", "sensor_gate", "barrier_gate", "reservoir_gates")


@dataclass(frozen=True)
class Check:
    name: str
    ok: bool
    detail: str = ""


class ContractError(AssertionError):
    pass


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def reference_device(seed: int = REFERENCE_SEED) -> Device:
    return compile_device(two_dot_layout("SiGe", 80.0), CompileParams(operating_point_v=REFERENCE_OPERATING_POINT_V), seed=seed)


def _v0(device: Device) -> np.ndarray:
    return np.array(device.programmed_point_v, dtype=float)


def _sweep(device: Device, n: int = 7, half_v: float = 0.01) -> np.ndarray:
    """A short sweep of the first plunger-like gate around the programmed point, (n, gates)."""
    v = np.broadcast_to(_v0(device), (n, len(device.gate_names))).copy()
    j = device.gate_names.index("P1") if "P1" in device.gate_names else 0
    v[:, j] += np.linspace(-half_v, half_v, n)
    return v


def _arrays_equal(a: object, b: object) -> bool:
    return all(np.array_equal(getattr(a, f.name), getattr(b, f.name)) for f in fields(a))


# --- the checks ------------------------------------------------------------------------------------


def check_compile_is_deterministic(_device: Device) -> str:
    a, b = reference_device(), reference_device()
    _require(a.dumps() == b.dumps(), "two compiles of the reference recipe differ")
    return "reference recipe compiled twice: identical text"


def check_round_trip(device: Device) -> str:
    text = device.dumps()
    again = Device.loads(text)
    _require(again.dumps() == text, "the JSON round trip is not text-identical")
    v = _sweep(device)
    _require(_arrays_equal(again.at(v), device.at(v)), "the reloaded device's State differs")
    m1 = again.session().sense(v, 1e-3)
    m2 = device.session().sense(v, 1e-3)
    _require(_arrays_equal(m1, m2), "the reloaded device's session differs")
    return f"{len(text)} bytes, reload answers bit for bit"


def check_provenance(device: Device) -> str:
    p = device.provenance
    _require(int(p["format_version"]) == DEVICE_FORMAT_VERSION, "format version")
    _require(isinstance(p.get("phantom_version"), str), "phantom_version missing")
    mv = p.get("model_versions", {})
    _require(len(mv) > 0 and all(isinstance(x, int) for x in mv.values()), "compile-time model versions missing")
    runtime = api.model_versions()
    for name in ("device.session", "device.truth", "readout.measure", "readout.noise", "readout.master_equation", "spin.exchange", "spin.psb", "optimum.targets"):
        _require(isinstance(runtime.get(name), int), f"model_versions() lacks {name}")
    return f"phantom {p['phantom_version']}, {len(mv)} compile-time and {len(runtime)} runtime model versions"


def check_at_shapes(device: Device) -> str:
    v0, g = _v0(device), len(device.gate_names)
    for lead in ((), (3,), (2, 3)):
        v = np.broadcast_to(v0, lead + (g,))
        s = device.at(v)
        _require(s.v_programmed_v.shape == lead + (g,) and s.v_carried_v.shape == lead + (g,), f"voltages for {lead}")
        _require(s.n_mean.shape == lead + (2,) and s.n_ground.shape == lead + (2,), f"charge for {lead}")
        for name in ("sensor_potential_v", "sensor_current_a", "tunnel_coupling_ev", "exchange_ev"):
            x = np.asarray(getattr(s, name))
            _require(x.shape == lead, f"{name} has shape {x.shape} for leading axes {lead}")
            _require(bool(np.all(np.isfinite(x))), f"{name} not finite")
        for name in ("gamma_left_hz", "gamma_right_hz", "gamma_interdot_hz"):
            _require(math.isfinite(getattr(s, name)) and getattr(s, name) >= 0, f"{name}")
    _require(bool(np.all(device.at(v0).n_ground == np.rint(device.at(v0).n_ground))), "n_ground is not integer")
    bias = np.array([0.0, 50e-6, 200e-6])
    i = device.at(np.broadcast_to(v0, (3, g)), v_sd=bias).sensor_current_a
    _require(i.shape == (3,) and len(set(i.tolist())) == 3, "the sensor bias does not reach the current")
    return "(), (3,), (2, 3) leading axes; finite; bias per point"


def check_sessions_replay(device: Device) -> str:
    v = _sweep(device)

    def script(run: api.Session) -> list:
        out = [run.sense(v, 1e-3)]
        run.hold(v[0], 2e-3)
        out.append(run.sense(v[-1], 1e-4, samples=5, v_sd=100e-6))
        run.wait(1e-3)
        out.append(run.spin_readout(0.5, 0.0, 16, 1e-5, 1e-3))
        return out

    def same(a: list, b: list) -> bool:
        return all(_arrays_equal(x, y) for x, y in zip(a, b, strict=True))

    _require(same(script(device.session()), script(device.session())), "two fresh sessions differ")
    _require(same(script(device.session(noise_seed=None)), script(device.session())), "noise_seed=None is not the default")
    _require(same(script(device.session(noise_seed=11)), script(device.session(noise_seed=11))), "a noise seed does not replay")
    draw = device.scenario.get("noise")
    if draw is None:
        return "replays; the device has no noise draw, so seeds cannot differ"
    a = device.session(noise_seed=11, charge="mean").sense(v, 1e-3).current_a
    b = device.session(noise_seed=12, charge="mean").sense(v, 1e-3).current_a
    _require(not np.array_equal(a, b), "two noise seeds give the same realisation")
    m1 = device.session(noise_seed=11).sense(v, 1e-3)
    _require(_arrays_equal(m1, device.measure(v, 1e-3, noise_seed=11)), "measure is not a fresh session's sense")
    return "fresh sessions replay; noise seeds replay and differ from each other"


def check_clock(device: Device) -> str:
    run = device.session()
    v = _sweep(device, 4)
    _require(run.time_s == 0.0, "a fresh session's clock is not zero")
    m = run.sense(v, 2e-3)
    _require(np.allclose(m.times_s, 2e-3 * np.arange(1, 5), rtol=0, atol=1e-15) and run.time_s == m.times_s[-1], "sample times")
    run.hold(v[0], 0.5)
    _require(math.isclose(run.time_s, 0.508, rel_tol=1e-12), "hold does not advance the clock")
    m = run.sense(v[:2], 1e-3, samples=3)
    _require(m.times_s.shape == (2, 3) and bool(np.all(np.diff(m.times_s.reshape(-1)) > 0)), "samples are not in order")
    _require(m.times_s[0, 0] > 0.508, "time went backwards")
    t = run.time_s
    run.sense(np.zeros((0, len(device.gate_names))), 1e-3)
    _require(run.time_s == t, "an empty sense moved the clock")
    for bad in (lambda: run.sense(v, 0.0), lambda: run.sense(v, -1e-3), lambda: run.hold(v[0], -1.0)):
        try:
            bad()
        except ValueError:
            continue
        raise ContractError("the clock accepted a step backwards or a zero dwell")
    return "samples at t0 + (k+1)·dwell; hold advances; no step backwards"


def check_equilibrium(device: Device) -> str:
    v = _sweep(device)
    eq = device.at(v).sensor_current_a
    fast = device.session(noise=False, dynamics=False, charge="mean").sense(v, 1e-3).current_a
    _require(np.array_equal(fast, eq), "noise-free mean session without dynamics differs from at()")
    slow = device.session(noise=False, dynamics=True, charge="mean").sense(v, 10.0).current_a
    _require(np.allclose(slow, eq, rtol=1e-6, atol=1e-18), f"long dwell differs from at() by {np.max(np.abs(slow - eq)):.3g} A")
    return "exact without dynamics; within 1e-6 after a 10 s dwell"


def check_bias(device: Device) -> str:
    v = np.broadcast_to(_v0(device), (3, len(device.gate_names)))
    bias = np.array([10e-6, 50e-6, 300e-6])
    m = device.session(noise=False, dynamics=False, charge="mean").sense(v, 1e-3, v_sd=bias)
    _require(np.array_equal(m.current_a, device.at(v, v_sd=bias).sensor_current_a), "sense does not pass the bias")
    m = device.session(noise=False, dynamics=False, charge="mean").sense(v[0], 1e-3, samples=2, v_sd=bias[2])
    # one point alone and in a batch agree within 1e-12 relative, not bit for bit (numpy's kernels)
    single = device.at(v[0], v_sd=bias[2]).sensor_current_a
    _require(np.allclose(m.current_a, single, rtol=BATCH_RTOL, atol=0), "scalar bias with samples")
    return "per-point and scalar bias reach the sensor"


def check_sample_semantics(device: Device) -> str:
    v = _sweep(device)
    for dynamics in (True, False):
        n = device.session(dynamics=dynamics).sense(v, 1e-3, samples=3).n_present
        _require(n.shape == (len(v), 3, 2) and bool(np.all(n == np.rint(n))), f"sampled charge is not integer (dynamics={dynamics})")
    n = device.session(noise=False, dynamics=False, charge="mean").sense(v, 1e-3).n_present
    _require(np.array_equal(n, device.at(v).n_mean), "mean charge is not ⟨N⟩")
    return "sampled: integer states; mean: ⟨N⟩"


def _pair_with_different_charge(device: Device) -> tuple[np.ndarray, np.ndarray] | None:
    """Two programmed points one electron apart along a plunger, both away from a degeneracy."""
    v0 = _v0(device)
    lattice = np.array(device.truth().charge_lattice_e_per_v)
    for k, gate in enumerate(device.truth().lattice_gates):
        per_e = lattice[k, k]
        if not per_e > 0:
            continue
        j = device.gate_names.index(gate)
        for start, step in itertools.product((0.0, 0.5, -0.5, 0.25, -0.25), (1.0, -1.0)):
            a = v0.copy()
            a[j] += start / per_e
            b = a.copy()
            b[j] += step / per_e
            sa, sb = device.at(a), device.at(b)
            clean = np.allclose(sa.n_mean, sa.n_ground, atol=1e-3) and np.allclose(sb.n_mean, sb.n_ground, atol=1e-3)
            if clean and not np.array_equal(sa.n_ground, sb.n_ground):
                return a, b
    return None


def check_hold_moves_the_charge(device: Device) -> str:
    pair = _pair_with_different_charge(device)
    if pair is None:
        return "skipped: no clean pair of points one electron apart near the programmed point"
    a, b = pair
    s = device.at(a)
    slowest = min(r for r in (s.gamma_left_hz, s.gamma_right_hz) if r > 0)
    fastest = max(s.gamma_left_hz, s.gamma_right_hz, s.gamma_interdot_hz)
    settle = 200.0 / slowest
    if settle > 1e6:
        return f"skipped: the slowest lead ({slowest:.3g} Hz) does not settle within 1e6 s"
    quick = 1e-6 / fastest
    run = device.session(noise=False, dynamics=True, charge="mean")
    run.hold(a, settle)
    n = run.sense(b, quick).n_present
    _require(np.allclose(n, device.at(a).n_ground, atol=1e-3), f"a fast read after holding at A shows {n}, not A's charge")
    run.hold(b, settle)
    n = run.sense(b, quick).n_present
    _require(np.allclose(n, device.at(b).n_ground, atol=1e-3), f"after holding at B the charge is {n}, not B's")
    return f"held {settle:.3g} s at each point; fast reads of {quick:.3g} s show the held charge"


def check_knowable(device: Device) -> str:
    k = device.knowable()
    d = k.to_dict()
    _require(tuple(d) == KNOWABLE_FIELDS, f"knowable fields {tuple(d)}")
    _require(tuple(k.programmed_point_v) == tuple(device.programmed_point_v), "programmed point")
    _require(tuple(k.gate_names) == tuple(device.gate_names), "gate names")
    json.dumps(d)
    return ", ".join(KNOWABLE_FIELDS)


def check_truth(device: Device) -> str:
    t = device.truth()
    _require(t.seed == device.seed and t.scenario == device.scenario, "seed and scenario")
    _require(len(t.dots) == 2 and len(t.valley_splittings_uev) == 2 and len(t.g_factors) == 2, "dots, valleys, g-factors")
    lat = np.array(t.charge_lattice_e_per_v)
    _require(lat.shape == (2, 2) and bool(np.all(np.isfinite(lat))), "charge lattice")
    json.dumps(t.to_dict())
    s = t.sensor
    if s is None:
        return "the sensor gate does not move the sensor: no sensor truth"
    _require(all(math.isfinite(getattr(s, f.name)) for f in fields(s) if f.name != "gate"), "sensor truth not finite")
    _require(math.isclose(s.wing_v - s.peak_v, math.atanh(1 / math.sqrt(3)) * s.width_v, rel_tol=1e-9), "the wing is not the inflection")
    _require(math.isclose(s.lever_arm_ev_per_v, s.charging_energy_ev / s.period_v, rel_tol=1e-12), "lever arm")
    s2 = device.truth(v_sd=2 * s.v_sd_v).sensor
    _require(s2 is not None and math.isclose(s2.sensitivity_a_per_v, 2 * s.sensitivity_a_per_v, rel_tol=1e-12), "sensitivity is not linear in bias")
    near = device.truth(near_v=s.peak_v + s.period_v).sensor
    _require(near is not None and math.isclose(near.peak_v, s.peak_v + s.period_v, rel_tol=0, abs_tol=1e-9 * s.period_v), "the next peak")
    return f"sensor peak {s.peak_v:.4f} V, wing {s.wing_v:.4f} V, period {s.period_v:.4f} V"


def check_optimum(device: Device) -> str:
    o = optimum(device)
    d = o.to_dict()
    _require(set(d) == {"sensor", "lattice", "coupling", "readout", "coherence"}, f"parts {sorted(d)}")
    flat = np.concatenate([np.ravel(np.asarray(x, dtype=float)) for part in d.values() for k, x in part.items() if k != "plungers"])
    _require(bool(np.all(np.isfinite(flat))), "a part of the optimum is not finite")
    return "sensor, lattice, coupling, readout, coherence"


def check_spin_readout(device: Device) -> str:
    run = device.session()
    out = run.spin_readout(0.0, 0.0, 50, 1e-5, 1e-3)
    _require(out.prepared_triplet.shape == (50,) and out.in_02.shape == (50,), "one outcome per shot")
    _require(out.prepared_triplet.dtype == bool and out.in_02.dtype == bool, "outcomes are bool")
    _require(not out.prepared_triplet.any(), "p_triplet = 0 prepared a triplet")
    _require(math.isclose(run.time_s, 50 * 1e-5, rel_tol=1e-12), "readout time is not on the clock")
    return "shots, dtypes, clock"


def check_surface(_device: Device) -> str:
    missing = [n for n in api.__all__ if not hasattr(api, n)]
    _require(not missing, f"phantom.api lacks {missing}")
    return f"{len(api.__all__)} public names"


CHECKS: tuple[Callable[[Device], str], ...] = (
    check_surface,
    check_compile_is_deterministic,
    check_round_trip,
    check_provenance,
    check_at_shapes,
    check_knowable,
    check_truth,
    check_sessions_replay,
    check_clock,
    check_equilibrium,
    check_bias,
    check_sample_semantics,
    check_hold_moves_the_charge,
    check_optimum,
    check_spin_readout,
)


def run(device: Device | None = None) -> list[Check]:
    """Every check on `device` (the reference device when None); a failing check never stops the rest."""
    device = reference_device() if device is None else device
    results = []
    for check in CHECKS:
        name = check.__name__.removeprefix("check_")
        try:
            results.append(Check(name, True, check(device)))
        except Exception as exc:  # noqa: BLE001 - a contract failure of any kind is reported, not raised
            results.append(Check(name, False, f"{type(exc).__name__}: {exc}"))
    return results


def main(argv: list[str] | None = None) -> int:
    paths = sys.argv[1:] if argv is None else argv
    targets: list[tuple[str, Device]] = [("reference device", reference_device())]
    targets += [(p, Device.from_json(p)) for p in paths]
    ok = True
    print(f"phantom {api.__version__} contract")
    for label, device in targets:
        print(f"{label} (seed {device.seed}):")
        for c in run(device):
            ok &= c.ok
            print(f"  {'ok  ' if c.ok else 'FAIL'} {c.name}: {c.detail}")
    print("contract: OK" if ok else "contract: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
