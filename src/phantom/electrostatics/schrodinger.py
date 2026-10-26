"""Two-dimensional effective-mass Schrödinger solve on a potential-energy window: dot levels and spacings.

Offline (compile tier, ADR 0001): a sparse finite-difference Hamiltonian on a rectangular window of
the grid, lowest eigenvalues by `scipy.sparse.linalg.eigsh`. Used to check the harmonic estimate of
`electrostatics.dots` and, later, by `compile`. Runtime layers never call it (rule 4).

    H = −(ħ² / 2m*) ∇² + U(x, y),  Dirichlet boundary at the window edge.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import diags, identity, kron
from scipy.sparse.linalg import eigsh

from phantom.electrostatics.dots import HBAR2_OVER_2ME_EV_NM2

MODEL_VERSION = 1


def dot_levels_ev(u_ev: np.ndarray, dx_nm: float, dy_nm: float, effective_mass_me: float, n_levels: int = 4) -> np.ndarray:
    """Lowest `n_levels` eigen-energies (eV, ascending) of an electron in the window `u_ev` (ny, nx)."""
    ny, nx = u_ev.shape
    t = HBAR2_OVER_2ME_EV_NM2 / effective_mass_me
    lap_x = diags([np.full(nx - 1, -1.0), np.full(nx, 2.0), np.full(nx - 1, -1.0)], [-1, 0, 1]) / dx_nm**2
    lap_y = diags([np.full(ny - 1, -1.0), np.full(ny, 2.0), np.full(ny - 1, -1.0)], [-1, 0, 1]) / dy_nm**2
    kinetic = t * (kron(identity(ny), lap_x) + kron(lap_y, identity(nx)))
    h = kinetic + diags(u_ev.reshape(-1))
    vals = eigsh(h.tocsc(), k=n_levels, sigma=float(u_ev.min()), which="LM", return_eigenvectors=False)
    return np.sort(vals)
