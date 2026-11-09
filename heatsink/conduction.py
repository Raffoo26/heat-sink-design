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


def base_plate_resistance(W: float, L: float, tb: float, k: float, a: float, b: float,
                          h_eff: float, n_xy: int = 60, n_z: int = 12) -> float:
    """Total resistance [K/W] from the mean source temperature to the ambient."""
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

    # surface temperature of the heated face (extrapolated half cell) weighted by the source
    T_top = T[:, :, -1] + q_flux * frac * dz / 2 / k
    T_mean_source = (T_top * frac).sum() / frac.sum()
    return T_mean_source / (4 * Q_quarter)
