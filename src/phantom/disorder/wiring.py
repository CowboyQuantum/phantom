"""Gate wiring faults: open gates, shorted gates, gate-to-gate leakage, and sensor-coupling asymmetry.

Catalog entries `gate-open-short` (#4), `gate-overlap-leakage` (#3) and `sensor-coupling-asymmetry`
(#11), all `defect`s with placeholder rates (`calibration: none`) until the owner writes them. These
faults need no new electrostatics: they transform the gate vector before the basis, add a current between
two gates, or scale the sensor's per-electron coupling.

- **Open gate**: the gate is disconnected and floats at a stuck potential; its programmed voltage has no
  effect. In the basis this is a column of the wiring matrix replaced by a constant.
- **Shorted gates**: two gates are electrically one; both follow the voltage programmed on the first.
  Their transition lines move together in every map, and the second gate's own sweep does nothing.
- **Gate-to-gate leakage**: a finite resistance between two gates. Its signature in `gate-leakage-iv`
  is a linear current between the pair that appears on both gates with opposite sign.
- **Sensor-coupling asymmetry**: the per-electron shift of the sensor potential differs between the
  dots (log-normal draws around the nominal), so one dot's lines are faint.

Severity multiplies the fault probabilities; severity 0 draws no faults and the wiring is the identity.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MODEL_VERSION = 1


@dataclass(frozen=True)
class WiringFaultParams:
    """Placeholder rates per device (`calibration: none`)."""

    p_open_per_gate: float = 0.03
    p_short_per_adjacent_pair: float = 0.03
    p_leak_per_adjacent_pair: float = 0.05
    leak_resistance_ohm_log10_mean: float = 9.0
    """log10 of the leak resistance in ohm: 1 GΩ typical when it happens."""
    leak_resistance_ohm_log10_sigma: float = 1.0
    stuck_potential_v: float = 0.0
    """Where an open gate floats."""
    sensor_beta_log_sigma: float = 0.35
    """Log-normal spread of each dot's per-electron sensor coupling around the nominal."""

    def __post_init__(self) -> None:
        for name in ("p_open_per_gate", "p_short_per_adjacent_pair", "p_leak_per_adjacent_pair"):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} must be a probability")


@dataclass(frozen=True)
class WiringDraw:
    gate_names: tuple[str, ...]
    open_gates: tuple[str, ...] = ()
    shorted_pairs: tuple[tuple[str, str], ...] = ()
    """(kept, follower): the follower takes the kept gate's voltage."""
    leaks: tuple[tuple[str, str, float], ...] = ()
    """(gate a, gate b, resistance in ohm)."""
    sensor_beta_factors: tuple[float, float] = (1.0, 1.0)
    stuck_potential_v: float = 0.0

    @property
    def is_clean(self) -> bool:
        return not self.open_gates and not self.shorted_pairs and not self.leaks and self.sensor_beta_factors == (1.0, 1.0)

    def to_dict(self) -> dict:
        return {
            "gate_names": list(self.gate_names),
            "open_gates": list(self.open_gates),
            "shorted_pairs": [list(p) for p in self.shorted_pairs],
            "leaks": [[a, b, r] for a, b, r in self.leaks],
            "sensor_beta_factors": list(self.sensor_beta_factors),
            "stuck_potential_v": self.stuck_potential_v,
        }

    @classmethod
    def from_dict(cls, d: dict) -> WiringDraw:
        return cls(
            tuple(d["gate_names"]),
            tuple(d.get("open_gates", [])),
            tuple((a, b) for a, b in d.get("shorted_pairs", [])),
            tuple((a, b, float(r)) for a, b, r in d.get("leaks", [])),
            (float(d["sensor_beta_factors"][0]), float(d["sensor_beta_factors"][1])),
            float(d.get("stuck_potential_v", 0.0)),
        )

    # --- effect on the gate vector ------------------------------------------------------------------

    def wiring_matrix(self) -> tuple[np.ndarray, np.ndarray]:
        """(M, c) such that the voltages the gates actually carry are M · V + c for programmed V."""
        n = len(self.gate_names)
        m = np.eye(n)
        c = np.zeros(n)
        for kept, follower in self.shorted_pairs:
            i, j = self.gate_names.index(kept), self.gate_names.index(follower)
            m[j, :] = m[i, :]
        for g in self.open_gates:
            j = self.gate_names.index(g)
            m[j, :] = 0.0
            c[j] = self.stuck_potential_v
        return m, c

    def apply(self, v: np.ndarray) -> np.ndarray:
        """Actual gate voltages (..., gates) for programmed `v` (..., gates)."""
        m, c = self.wiring_matrix()
        return np.asarray(v, dtype=float) @ m.T + c

    def leakage_current_a(self, v: np.ndarray) -> np.ndarray:
        """Current (..., gates) flowing *out of* each gate through the leaks, for actual voltages `v`."""
        v = np.asarray(v, dtype=float)
        out = np.zeros_like(v)
        for a, b, r in self.leaks:
            i, j = self.gate_names.index(a), self.gate_names.index(b)
            i_leak = (v[..., i] - v[..., j]) / r
            out[..., i] += i_leak
            out[..., j] -= i_leak
        return out


def draw_wiring_faults(
    params: WiringFaultParams,
    gate_names: tuple[str, ...],
    adjacent_pairs: tuple[tuple[str, str], ...],
    rng: np.random.Generator,
    severity: float = 1.0,
) -> WiringDraw:
    """Independent Bernoulli faults at severity × probability; a gate is open or in at most one short."""
    if severity <= 0:
        return WiringDraw(gate_names, stuck_potential_v=params.stuck_potential_v)
    s = min(severity, 1.0 / max(params.p_open_per_gate, params.p_short_per_adjacent_pair, params.p_leak_per_adjacent_pair, 1e-12))
    open_gates = tuple(g for g in gate_names if rng.random() < s * params.p_open_per_gate)
    used = set(open_gates)
    shorts = []
    for a, b in adjacent_pairs:
        if rng.random() < s * params.p_short_per_adjacent_pair and a not in used and b not in used:
            shorts.append((a, b))
            used.update((a, b))
    leaks = []
    for a, b in adjacent_pairs:
        if rng.random() < s * params.p_leak_per_adjacent_pair:
            r = 10 ** rng.normal(params.leak_resistance_ohm_log10_mean, params.leak_resistance_ohm_log10_sigma)
            leaks.append((a, b, float(r)))
    factors = np.exp(rng.normal(0.0, severity * params.sensor_beta_log_sigma, 2))
    return WiringDraw(gate_names, open_gates, tuple(shorts), tuple(leaks), (float(factors[0]), float(factors[1])), params.stuck_potential_v)
