"""
3-D finite-volume conduction solver for the heat sink base, used to verify the
spreading-resistance correlation.

Rectangular plate W × L × tb, uniform heat flux over a centred source on one face,
uniform heat transfer coefficient h_eff on the opposite face (representing the
fins), adiabatic edges. Only a quarter of the plate is solved (symmetry).
"""

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def _solve_quarter(W, L, tb, k, a, b, h_eff, n_xy, n_z):
    """Temperature rise per watt in a quarter of the plate; returns the field and helpers."""
    nx, ny, nz = n_xy, n_xy, n_z
    dx, dy, dz = W / 2 / nx, L / 2 / ny, tb / nz
    xc = (np.arange(nx) + 0.5) * dx            # cell centres from the plate centre
    yc = (np.arange(ny) + 0.5) * dy
    idx = np.arange(nx * ny * nz).reshape(nx, ny, nz)

    rows, cols, vals = [], [], []
    rhs = np.zeros(nx * ny * nz)
    diag = np.zeros(nx * ny * nz)

    def link(c1, c2, G):
        rows.extend([c1, c2, c1, c2]); cols.extend([c2, c1, c1, c2]); vals.extend([-G, -G, 0, 0])
        diag[c1] += G; diag[c2] += G

    Gx, Gy, Gz = k * dy * dz / dx, k * dx * dz / dy, k * dx * dy / dz
    link_arrays = [(idx[:-1, :, :], idx[1:, :, :], Gx), (idx[:, :-1, :], idx[:, 1:, :], Gy),
                   (idx[:, :, :-1], idx[:, :, 1:], Gz)]
    for A, B, G in link_arrays:
        for c1, c2 in zip(A.ravel(), B.ravel()):
            link(c1, c2, G)

    # cooled face z = 0: conductance in series (half cell + film) to ambient at 0
    G_film = 1.0 / (dz / 2 / (k * dx * dy) + 1.0 / (h_eff * dx * dy))
    diag[idx[:, :, 0].ravel()] += G_film

    # heated face z = tb: source fraction of each cell (exact for partially covered cells)
    fx = np.clip((a / 2 - (xc - dx / 2)) / dx, 0, 1)
    fy = np.clip((b / 2 - (yc - dy / 2)) / dy, 0, 1)
    frac = np.outer(fx, fy)
    Q_quarter = 1.0
    q_flux = Q_quarter / (a / 2 * b / 2)
    rhs[idx[:, :, -1].ravel()] = (q_flux * frac * dx * dy).ravel()

    n = nx * ny * nz
    M = coo_matrix((vals + list(diag), (rows + list(range(n)), cols + list(range(n)))), shape=(n, n)).tocsr()
    T = spsolve(M, rhs).reshape(nx, ny, nz)

    return T / (4 * Q_quarter), xc, yc, dz, frac, q_flux / (4 * Q_quarter)


def base_plate_resistance(W: float, L: float, tb: float, k: float, a: float, b: float,
                          h_eff: float, n_xy: int = 60, n_z: int = 12) -> float:
    """Total resistance [K/W] from the mean source temperature to the ambient."""
    T, _, _, dz, frac, q = _solve_quarter(W, L, tb, k, a, b, h_eff, n_xy, n_z)
    T_top = T[:, :, -1] + q * frac * dz / 2 / k        # heated face (extrapolated half cell)
    return float((T_top * frac).sum() / frac.sum())


def base_temperature_field(W: float, L: float, tb: float, k: float, a: float, b: float,
                           h_eff: float, n_xy: int = 40, n_z: int = 8):
    """Temperature rise per watt [K/W] in the whole plate, as (x, y, z, T[x, y, z]).

    x across W and y along L are measured from the plate centre; z from the cooled face (z = 0)
    to the heated face (z = tb)."""
    Tq, xc, yc, dz, _, _ = _solve_quarter(W, L, tb, k, a, b, h_eff, n_xy, n_z)
    T = np.concatenate([Tq[::-1], Tq], axis=0)
    T = np.concatenate([T[:, ::-1], T], axis=1)
    x = np.concatenate([-xc[::-1], xc])
    y = np.concatenate([-yc[::-1], yc])
    z = (np.arange(n_z) + 0.5) * dz
    return x, y, z, T
