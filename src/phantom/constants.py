"""Physical constants (CODATA 2018), below every layer so any layer may import them (ADR 0004)."""

E_CHARGE_C = 1.602176634e-19
"""Elementary charge in coulomb (exact)."""
EPSILON_0_F_PER_M = 8.8541878128e-12
"""Vacuum permittivity in F/m."""
K_B_EV_PER_K = 8.617333262e-5
"""Boltzmann constant in eV/K."""
HBAR_EV_S = 6.582119569e-16
"""Reduced Planck constant in eV·s."""
HBAR2_OVER_2ME_EV_NM2 = 0.0380998212
"""ħ² / (2 m_e) in eV·nm²."""
E_UNIT_V_PER_E_AF = E_CHARGE_C / 1e-18
"""Volts per electron per attofarad: e / (1 aF). Converts aF·V to electrons and e²/aF to eV."""
