"""The runtime surface of a compiled device: `State` at a gate point and a `Session` with a clock (plan §4, ADR 0024).

`state_at(device, v)` is a pure function of the programmed gate voltages (..., gates): the voltages
the gates carry after the wiring draw, the thermal and ground-state charge, the sensor potential and
its noise-free current at a source–drain bias, the interdot coupling and the exchange energy there,
the dots' lead and interdot rates there (#53), and the sensor's rates through its barriers and baths. Nothing
in it has a time or a random number.

A `Session` is one run on the device. It owns a clock that only moves forward, the device's noise
processes, and the charge on the dots. Its three verbs are all a consumer needs to drive the device:

- `sense(v, dwell_s, samples, v_sd)`: each gate point is held for `samples` dwells (one when
  `samples` is None) and read once at the end of each dwell. Sample k of a call is read at
  t₀ + (k + 1)·dwell, t₀ the clock when the call starts; the clock ends at the last sample.
- `hold(v, duration_s)`: the gates sit at `v` for `duration_s`; the charge evolves there, nothing is
  read. The gate noise is sampled once at the start of the hold and held through it (quasi-static).
- `wait(duration_s)`: `hold` at the last point the session was given (the programmed point if none).

**What one sample is.** The noise processes are sampled at the sample's time: the common-mode gate
noise is added to the carried voltages, the sensor's switchers shift its potential, and white readout
noise is added to the current with its density fixed at the device's 1 ms value (σ scales as
√(1 ms / dwell)). The charge that enters the sensor depends on two choices, fixed per session:

- `charge="sampled"` (default): one charge state per sample. With `dynamics=True` it is the state of
  a Markov jump trajectory of the rate equation (`readout.master_equation`): from the state present,
  the transition probabilities over the dwell are the column of exp(W·dwell), and the next state is
  drawn from them, so a sample is correlated with the last one when the rates are slow against the
  dwell and is a fresh Boltzmann draw when they are fast. The rates are those at the sample's
  voltages, gate noise included (#53): closing a barrier gate slows its lead until the charge
  latches along a sweep and telegraphs at a parked point. With `dynamics=False` every sample is an
  independent Boltzmann draw at its voltages (qdt-sim's behaviour). `Measurement.n_present` holds
  the integer state.
- `charge="mean"`: the ensemble average. With `dynamics=True` the distribution is propagated with
  exp(W·dwell) and ⟨N⟩ enters the sensor (latching and lag, as the Step 4 protocols); with
  `dynamics=False` it is the thermal ⟨N⟩ at the voltages, so a noise-free session reads exactly
  `State.sensor_current_a`. `n_present` holds ⟨N⟩.

Before its first call the device has sat at the first point long enough to be in equilibrium there.

**Seeds.** The device seed is the device's identity. `noise_seed=None` replays the device's own noise
realisation (`scenario.noise.seed`), so two fresh sessions given identical calls return identical
results. An integer `noise_seed` keeps every noise parameter of the draw (amplitudes, band, switcher
rates, white floor) and gives a new realisation: the noise draw's seed becomes a hash of
(device seed, noise seed), and the charge and spin samplers use streams keyed by both.

`spin_readout` is the Pauli-blockade readout of a prepared pair as a protocol step: it does not read
the session's charge or voltages (a spin state driven by voltages over time is #41). It takes
`shots × t_integration_s` on the clock, during which the charge relaxes at the last point.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Literal

import numpy as np
from scipy.linalg import expm

from phantom.constants import E_CHARGE_C, K_B_EV_PER_K
from phantom.readout import sensor_dynamics
from phantom.readout.master_equation import PROBABILITY_FLOOR, RateEquation
from phantom.readout.measure import SensedDoubleDot
from phantom.spin.exchange import exchange_j_ev
from phantom.spin.psb import psb_occupation

if TYPE_CHECKING:
    from phantom.device import Device

MODEL_VERSION = 4
"""Bumped whenever any output of `state_at` or `Session` changes for the same device, calls and seeds.
2: the sensor reads its leads at the carried voltages; `State.sensor_gamma_hz`; the switchers in the
sensor dot's units (ADR 0026). 3: the dots' rates follow the gates, in `State` and in the dynamics (#53).
4: with dynamics, the sensor dot's charge is a process on the clock that acts back on the dots (#54)."""

SESSION_STREAM = 8
"""Key of the session's random streams under the device seed: [seed, 8, (noise_seed,) 0] samples the
charge, [..., 1] the spin readout, [..., 2] the sensor dot's charge and counted electrons."""

WHITE_REFERENCE_DWELL_S = 1e-3
"""The dwell at which the noise draw states the white readout floor (`white_sigma_a_at_1ms`)."""

ChargeMode = Literal["sampled", "mean"]


@dataclass(frozen=True)
class State:
    """The device at programmed gate voltages, noise-free; arrays carry the leading axes of the input."""

    v_programmed_v: np.ndarray
    """(..., gates) as given, in V."""
    v_carried_v: np.ndarray
    """(..., gates) the voltages the gates carry after the wiring draw, in V."""
    n_mean: np.ndarray
    """(..., 2) thermal-average occupation of (dot 1, dot 2), electrons."""
    n_ground: np.ndarray
    """(..., 2) lowest-energy integer occupation, electrons."""
    sensor_potential_v: np.ndarray
    """(...) sensor potential φ_s with the thermal occupation, V-equivalent."""
    sensor_current_a: np.ndarray
    """(...) noise-free sensor current at `v_sd_v`, A."""
    v_sd_v: float | np.ndarray
    """Source–drain bias of the sensor, V (the sensor's default bias when none was given)."""
    tunnel_coupling_ev: np.ndarray
    """(...) interdot tunnel coupling t at the carried barrier-gate voltage, eV."""
    exchange_ev: np.ndarray
    """(...) J of the two-electron sector at the carried voltages, eV; meaningful in the (1,1) cell."""
    gamma_left_hz: np.ndarray
    """(...) dot 1 to the left lead at the carried voltages, Hz: the compiled rate moved by its gate
    slopes, capped at the attempt rate (#53)."""
    gamma_right_hz: np.ndarray
    """(...) dot 2 to the right lead at the carried voltages, Hz."""
    gamma_interdot_hz: np.ndarray
    """(...) electron transfer between the dots at the carried voltages, Hz."""
    sensor_gamma_hz: np.ndarray | None = None
    """(..., 2) the sensor's source and drain rates at the carried voltages, Hz, through its barriers
    and baths (ADR 0026); None for a typed sensor without leads."""


@dataclass(frozen=True)
class Measurement:
    """What `Session.sense` returns; leading axes are the input's, then `samples` when given."""

    times_s: np.ndarray
    """Time of each sample on the session clock, s."""
    current_a: np.ndarray
    """Sensor current with noise, A."""
    n_present: np.ndarray
    """(..., 2) the charge that entered the sensor: an integer state (`charge="sampled"`) or ⟨N⟩."""
    sensor_charge: np.ndarray | None = None
    """(...) the sensor dot's charge at the end of each dwell, electrons on the level numbering of φ_s
    (only its changes are physical): an integer (`sampled`) or its mean; None without dynamics or for a
    typed sensor without leads (ADR 0028)."""


@dataclass(frozen=True)
class SpinOutcome:
    """Per shot of `Session.spin_readout`."""

    prepared_triplet: np.ndarray
    """(shots,) bool: the pair was prepared in a triplet."""
    in_02: np.ndarray
    """(shots,) bool: the charge read (0,2) (a singlet signature) at the end of the integration."""


def _check_gates(device: Device, v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    if v.ndim == 0 or v.shape[-1] != len(device.gate_names):
        raise ValueError(f"gate voltages need a last axis of {len(device.gate_names)}: {device.gate_names}")
    return v


def state_at(device: Device, v_programmed: np.ndarray, v_sd: float | np.ndarray | None = None) -> State:
    """The noise-free `State` at programmed voltages (..., gates); `v_sd` broadcasts against (...)."""
    v = _check_gates(device, v_programmed)
    vc = device.effective_voltages(v)
    sensed = device.runtime_sensed(noise=False)
    ci, sensor = sensed.charge, sensed.sensor
    n_mean = ci.mean_occupation(vc)
    bias = sensor.bias_v if v_sd is None else v_sd
    hub = device.runtime_hubbard()
    gammas = device.runtime_rates().rates_hz(vc)
    return State(
        v_programmed_v=v,
        v_carried_v=vc,
        n_mean=n_mean,
        n_ground=ci.ground_state_occupation(vc),
        sensor_potential_v=sensor.potential(vc, n_mean),
        sensor_current_a=sensed.current_a(vc, bias),
        v_sd_v=bias,
        tunnel_coupling_ev=hub.barrier.tunnel_coupling_ev(vc[..., hub.gate_names.index(hub.barrier_gate)]),
        exchange_ev=exchange_j_ev(hub, vc),
        gamma_left_hz=gammas[..., 0],
        gamma_right_hz=gammas[..., 1],
        gamma_interdot_hz=gammas[..., 2],
        sensor_gamma_hz=sensor.rates_hz(vc),
    )


def noise_realisation_seed(device_seed: int, noise_seed: int) -> int:
    """The noise draw's seed for a session with `noise_seed` on the device with `device_seed`."""
    return int(np.random.SeedSequence([device_seed, SESSION_STREAM, noise_seed]).generate_state(1)[0])


class Session:
    """One run on a compiled device: a forward-only clock, the device's noise, and the dot charge."""

    def __init__(
        self,
        device: Device,
        noise: bool = True,
        dynamics: bool = True,
        charge: ChargeMode = "sampled",
        noise_seed: int | None = None,
    ) -> None:
        if charge not in ("sampled", "mean"):
            raise ValueError("charge must be 'sampled' or 'mean'")
        if noise_seed is not None and int(noise_seed) < 0:
            raise ValueError("noise_seed must be a non-negative integer")
        self.device = device
        self.noise = bool(noise)
        self.dynamics = bool(dynamics)
        self.charge: ChargeMode = charge
        self.noise_seed = None if noise_seed is None else int(noise_seed)
        draw = device.scenario_draw().noise
        ci, sensor = device.runtime_ci(), device.runtime_sensor()
        if self.noise and draw is not None:
            if self.noise_seed is not None:
                draw = replace(draw, seed=noise_realisation_seed(device.seed, self.noise_seed))
            self._sensed = SensedDoubleDot.from_noise_draw(draw, WHITE_REFERENCE_DWELL_S, device.switcher_scale(), charge=ci, sensor=sensor)
        else:
            self._sensed = SensedDoubleDot(charge=ci, sensor=sensor)
        self._rates: RateEquation = device.runtime_rates()
        key = [int(device.seed), SESSION_STREAM] + ([] if self.noise_seed is None else [self.noise_seed])
        self._charge_rng = np.random.default_rng(key + [0])
        self._spin_rng = np.random.default_rng(key + [1])
        self._sensor_rng = np.random.default_rng(key + [2])
        self._sensor_dynamics = self.dynamics and sensor.leads is not None
        self._m_states: np.ndarray | None = None
        self._m_p: np.ndarray | None = None
        self._bias_last = float(sensor.bias_v)
        self._t = 0.0
        self._v_last: np.ndarray | None = None
        self._states: np.ndarray | None = None
        self._p: np.ndarray | None = None

    # --- the clock --------------------------------------------------------------------------------

    @property
    def time_s(self) -> float:
        """The session clock, s: zero at the start, the time of the last sample or the end of the last hold."""
        return self._t

    @property
    def last_point_v(self) -> np.ndarray:
        """The programmed voltages the gates were last given (the device's programmed point if none)."""
        return (self._v_last if self._v_last is not None else np.array(self.device.programmed_point_v, dtype=float)).copy()

    # --- charge -----------------------------------------------------------------------------------

    def _start(self, v_carried: np.ndarray) -> None:
        """Equilibrium at the first point: the Boltzmann distribution (the rate equation's stationary
        state, by detailed balance), or one state drawn from it."""
        states, p = self._sensed.charge.occupation_probabilities(v_carried)
        keep = p > PROBABILITY_FLOOR
        states, p = states[keep], p[keep] / p[keep].sum()
        if self.charge == "sampled":
            i = int(np.searchsorted(np.cumsum(p), self._charge_rng.random() * p.sum(), side="right"))
            i = min(i, len(p) - 1)
            states, p = states[i : i + 1], np.ones(1)
        self._states, self._p = np.asarray(states, dtype=int), p

    def _propagate(self, v_carried: np.ndarray, dt_s: float, level_shift_ev: np.ndarray | None = None) -> np.ndarray:
        """Evolve the charge for `dt_s` at `v_carried` (gates,), each dot's levels shifted by
        `level_shift_ev` (2,) per electron; return what enters the sensor (2,)."""
        if self._states is None:
            self._start(v_carried)
        states, p = self._states, self._p
        keep = states[p > PROBABILITY_FLOOR]
        active = np.unique(np.concatenate([self._rates.window_states(v_carried), keep], axis=0), axis=0)
        index = {tuple(s): i for i, s in enumerate(active.tolist())}
        q = np.zeros(len(active))
        for s, prob in zip(states.tolist(), p, strict=True):
            if prob > PROBABILITY_FLOOR:
                q[index[tuple(s)]] += prob
        q /= q.sum()
        if dt_s > 0:
            _, w = self._rates.rate_matrix(v_carried, active, level_shift_ev)
            q = np.maximum(expm(w * dt_s) @ q, 0.0)
            q /= q.sum()
        if self.charge == "sampled":
            i = int(np.searchsorted(np.cumsum(q), self._charge_rng.random(), side="right"))
            i = min(i, len(q) - 1)
            self._states, self._p = active[i : i + 1], np.ones(1)
            return active[i].astype(float)
        self._states, self._p = active, q
        return q @ active

    # --- the sensor dot's charge (ADR 0028) -----------------------------------------------------------

    def _back_action_ev(self, v_carried: np.ndarray, delta_phi_v: float) -> np.ndarray | None:
        """(2,) shift of each dot's levels per electron from the sensor's charge, eV: −β_i·(M − M̄),
        M the sensor's charge (its mean in `mean` mode) and M̄ the mean of its stationary state at these
        voltages and the dots' present charge. The closed-form potential, which has no electrons, stands
        for the sensor at its mean; its switching and its lag are what act back (ADR 0028)."""
        if not self._sensor_dynamics or self._m_states is None or self._states is None:
            return None
        sensor = self._sensed.sensor
        n = self._p @ self._states
        phi = float(sensor.potential(v_carried, n, delta_phi_v))
        lo, hi = sensor_dynamics.window_states(sensor, phi, self._bias_last)
        c = sensor_dynamics.chain(sensor, phi, v_carried, self._bias_last, lo, hi)
        m_bar = float(sensor_dynamics.stationary(c) @ c.states)
        m = float(self._m_p @ self._m_states)
        return -sensor.beta_v_per_e * sensor.lever_arm_ev_per_v * (m - m_bar)

    def _sensor_step(self, phi: float, v_carried: np.ndarray, v_sd: float, dt_s: float) -> tuple[float, float]:
        """Evolve the sensor's charge for `dt_s` at potential `phi` (with the dots' charge and the
        switchers), carried voltages and bias; return (current over the dwell in A, charge after)."""
        sensor = self._sensed.sensor
        lo, hi = sensor_dynamics.window_states(sensor, phi, v_sd)
        rng = self._sensor_rng
        if self._m_states is None:
            c = sensor_dynamics.chain(sensor, phi, v_carried, v_sd, lo, hi)
            p = sensor_dynamics.stationary(c)
            if self.charge == "sampled":
                i = min(int(np.searchsorted(np.cumsum(p), rng.random(), side="right")), len(p) - 1)
                self._m_states, self._m_p = c.states[i : i + 1], np.ones(1)
            else:
                self._m_states, self._m_p = c.states, p
        held = self._m_states[self._m_p > PROBABILITY_FLOOR]
        lo, hi = min(lo, int(held.min())), max(hi, int(held.max()))
        c = sensor_dynamics.chain(sensor, phi, v_carried, v_sd, lo, hi)
        q = np.zeros(hi - lo + 1)
        for m, prob in zip(self._m_states.tolist(), self._m_p, strict=True):
            if prob > PROBABILITY_FLOOR:
                q[m - lo] += prob
        q /= q.sum()
        if self.charge == "mean":
            after, counted = sensor_dynamics.propagate(c, q, dt_s)
            self._m_states, self._m_p = c.states, after
            return E_CHARGE_C * counted / dt_s + c.background_a, float(after @ c.states)
        m0 = int(self._m_states[0])
        if sensor_dynamics.expected_jumps(c, dt_s) <= sensor_dynamics.EXACT_JUMPS:
            counted, m = sensor_dynamics.trajectory(sensor, phi, v_carried, v_sd, m0, dt_s, rng)
        else:
            after, mean = sensor_dynamics.propagate(c, q, dt_s)
            _, var = sensor_dynamics.counting_cumulants(c)
            counted = mean + np.sqrt(max(var, 0.0) * dt_s) * rng.standard_normal()
            i = min(int(np.searchsorted(np.cumsum(after), rng.random(), side="right")), len(after) - 1)
            m = int(c.states[i])
        self._m_states, self._m_p = np.array([m]), np.ones(1)
        return E_CHARGE_C * counted / dt_s + c.background_a, float(m)

    def _equilibrium(self, v_carried: np.ndarray) -> np.ndarray:
        """Charge without dynamics for points (n, gates): independent Boltzmann draws or ⟨N⟩."""
        ci = self._sensed.charge
        if self.charge == "mean":
            return ci.mean_occupation(v_carried)
        states, p = ci.occupation_probabilities(v_carried)
        u = self._charge_rng.random(len(v_carried))
        idx = (np.cumsum(p, axis=-1) <= u[:, None]).sum(axis=-1)
        idx = np.minimum(idx, p.shape[-1] - 1)
        return np.take_along_axis(states, idx[:, None, None], axis=-2)[:, 0, :].astype(float)

    # --- the verbs ----------------------------------------------------------------------------------

    def sense(
        self,
        v_programmed: np.ndarray,
        dwell_s: float,
        samples: int | None = None,
        v_sd: float | np.ndarray | None = None,
    ) -> Measurement:
        """Read the sensor at programmed voltages (..., gates), points in C order, each held for
        `samples` dwells (one when None) and read at the end of each; `v_sd` (V) broadcasts against
        (...), the sensor's default bias when None. Output arrays have shape (...) or (..., samples)."""
        v = _check_gates(self.device, v_programmed)
        if not dwell_s > 0:
            raise ValueError("dwell_s must be positive: the clock only moves forward")
        if samples is not None and int(samples) < 1:
            raise ValueError("samples must be at least 1")
        reps = 1 if samples is None else int(samples)
        lead = v.shape[:-1]
        out_shape = lead + (() if samples is None else (reps,))
        pts = np.repeat(v.reshape(-1, v.shape[-1]), reps, axis=0)
        n = len(pts)
        if n == 0:
            return Measurement(np.zeros(out_shape), np.zeros(out_shape), np.zeros(out_shape + (2,)))
        bias = self._sensed.sensor.bias_v if v_sd is None else v_sd
        bias = np.repeat(np.broadcast_to(np.asarray(bias, dtype=float), lead).reshape(-1), reps)
        times = self._t + dwell_s * np.arange(1, n + 1)
        g = self._sensed.gate_noise.trace(times)
        s = self._sensed.sensor_noise.trace(times)
        r = self._sensed.readout_noise.trace(times) * np.sqrt(WHITE_REFERENCE_DWELL_S / dwell_s)
        v_eff = self.device.effective_voltages(pts) + g[:, None]
        sensor = self._sensed.sensor
        m_sensor = None
        if self._sensor_dynamics:
            n_present, current, m_sensor = np.empty((n, 2)), np.empty(n), np.empty(n)
            for k in range(n):
                self._bias_last = float(bias[k])
                n_present[k] = self._propagate(v_eff[k], dwell_s, self._back_action_ev(v_eff[k], float(s[k])))
                phi = float(sensor.potential(v_eff[k], n_present[k], s[k]))
                current[k], m_sensor[k] = self._sensor_step(phi, v_eff[k], float(bias[k]), dwell_s)
            current = current + r
            m_sensor = m_sensor.reshape(out_shape)
        else:
            if self.dynamics:
                n_present = np.empty((n, 2))
                for k in range(n):
                    n_present[k] = self._propagate(v_eff[k], dwell_s)
            else:
                n_present = self._equilibrium(v_eff)
            current = sensor.current_a(sensor.potential(v_eff, n_present, s), bias, v_eff) + r
        self._t = float(times[-1])
        self._v_last = pts[-1].copy()
        return Measurement(times.reshape(out_shape), current.reshape(out_shape), n_present.reshape(out_shape + (2,)), m_sensor)

    def hold(self, v_programmed: np.ndarray, duration_s: float) -> None:
        """Keep the gates at programmed `v` (gates,) for `duration_s`: the clock advances, the charge
        evolves at those voltages (with the gate noise of the hold's start), nothing is read."""
        v = _check_gates(self.device, v_programmed)
        if v.ndim != 1:
            raise ValueError("hold takes one gate point (gates,)")
        if duration_s < 0:
            raise ValueError("duration_s must be non-negative: the clock only moves forward")
        if self.dynamics and duration_s > 0:
            t = np.array([self._t])
            g = float(self._sensed.gate_noise.trace(t)[0])
            v_eff = self.device.effective_voltages(v) + g
            if self._sensor_dynamics:
                s = float(self._sensed.sensor_noise.trace(t)[0])
                n = self._propagate(v_eff, float(duration_s), self._back_action_ev(v_eff, s))
                self._sensor_step(float(self._sensed.sensor.potential(v_eff, n, s)), v_eff, self._bias_last, float(duration_s))
            else:
                self._propagate(v_eff, float(duration_s))
        self._t += float(duration_s)
        self._v_last = v.copy()

    def wait(self, duration_s: float) -> None:
        """`hold` at the last point the session was given (the programmed point if none)."""
        self.hold(self.last_point_v, duration_s)

    def spin_readout(
        self,
        p_triplet: float,
        detuning_ev: float,
        shots: int,
        t_integration_s: float,
        t1_read_s: float,
    ) -> SpinOutcome:
        """Pauli-blockade readout of `shots` pairs, each prepared as a triplet with probability
        `p_triplet` and read at detuning ε (eV, from the (1,1)–(0,2) anticrossing) for
        `t_integration_s`: a triplet that relaxes (T₁ = `t1_read_s`) during the integration reads as
        a singlet; the (0,2) probability of each spin state follows `spin.psb.psb_occupation` with the
        receiving dot's valley splitting, the compiled coupling and the electron temperature."""
        if not 0.0 <= p_triplet <= 1.0:
            raise ValueError("p_triplet must be a probability")
        if shots < 0 or t_integration_s < 0 or not t1_read_s > 0:
            raise ValueError("need shots >= 0, t_integration_s >= 0 and t1_read_s > 0")
        dev = self.device
        kt = K_B_EV_PER_K * dev.temperature_k
        e_v = dev.dots[1].valley_splitting_uev * 1e-6
        p_s, p_t = psb_occupation(np.array([float(detuning_ev)]), e_v, dev.tunnel_coupling_ev, kt)
        prepared = self._spin_rng.random(shots) < p_triplet
        relaxed = self._spin_rng.random(shots) < -np.expm1(-t_integration_s / t1_read_s)
        p02 = np.where(prepared & ~relaxed, p_t[0], p_s[0])
        in_02 = self._spin_rng.random(shots) < p02
        self.wait(shots * t_integration_s)
        return SpinOutcome(prepared, in_02)


__all__ = ["MODEL_VERSION", "Measurement", "Session", "SpinOutcome", "State", "noise_realisation_seed", "state_at"]
