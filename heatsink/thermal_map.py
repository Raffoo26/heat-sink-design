"""
Temperature field of a solved heat sink, for colouring 3-D pictures.

Base: 3-D conduction solution (heatsink.conduction) with the fins replaced by the
effective heat transfer coefficient of the model. Fins: 1-D fin solution from the
local root temperature, θ(s) = θ_root · cosh(m (H_c − s)) / cosh(m H_c).
This is a model-based picture, not a CFD result: h is uniform over the fins.
"""

import math

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .conduction import base_temperature_field
from .correlations import fin_m
from .model import HeatSinkResult


def temperature_function(r: HeatSinkResult, n_xy: int = 40, n_z: int = 8):
    """Return T(P) [°C] for points P (n × 3) in millimetres, in the STL coordinates of
    heatsink.export: x across the fins, y up from the bottom of the base, z along the fins."""
    hs, src = r.hs, r.source
    h_eff = 1 / (r.R_sink * hs.A_footprint)
    x, y, z, field = base_temperature_field(hs.W, hs.L, hs.tb, hs.k, src.a, src.b, h_eff, n_xy, n_z)
    T_amb = r.T_amb - 273.15
    base = RegularGridInterpolator((x, y, z), field * r.Q, bounds_error=False, fill_value=None)
    root = RegularGridInterpolator((x, y), field[:, :, 0] * r.Q, bounds_error=False, fill_value=None)
    m = fin_m(r.h, hs.k, hs.t)
    Hc = hs.H + (hs.t / 2 if r.mode == "natural" else 0.0)

    def T(P):
        P = np.asarray(P, float) / 1e3
        xc, height, yc = P[:, 0] - hs.W / 2, P[:, 1], P[:, 2] - hs.L / 2
        xc, yc = np.clip(xc, x[0], x[-1]), np.clip(yc, y[0], y[-1])
        out = np.empty(len(P))
        in_base = height <= hs.tb + 1e-9
        zc = np.clip(hs.tb - height[in_base], z[0], z[-1])
        out[in_base] = base(np.c_[xc[in_base], yc[in_base], zc])
        s = np.clip(height[~in_base] - hs.tb, 0, hs.H)
        theta_root = root(np.c_[xc[~in_base], yc[~in_base]])
        out[~in_base] = theta_root * np.cosh(m * (Hc - s)) / math.cosh(m * Hc)
        return T_amb + out

    return T
