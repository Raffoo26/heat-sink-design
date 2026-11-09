"""Plots: heat sink cross-section, design sweeps, fan operating point, resistance breakdown."""

import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

from .air import KELVIN
from .model import Fan, HeatSinkResult, forced_convection


def plot_profile(hs, source=None, ax=None, title=None):
    """Cross-section of the heat sink (perpendicular to the fins), dimensions in mm."""
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 3.6))
    mm = 1e3
    ax.add_patch(Rectangle((0, 0), hs.W * mm, hs.tb * mm, color="#8a9bb0"))
    for i in range(hs.N):
        x = i * hs.pitch * mm
        ax.add_patch(Rectangle((x, hs.tb * mm), hs.t * mm, hs.H * mm, color="#8a9bb0"))
    if source is not None:
        x0 = (hs.W - source.a) / 2 * mm
        ax.add_patch(Rectangle((x0, -0.12 * hs.H * mm), source.a * mm, 0.12 * hs.H * mm, color="#d62728",
                               label="heat source"))
        ax.legend(loc="upper right", fontsize=9)
    ax.set_xlim(-0.05 * hs.W * mm, 1.05 * hs.W * mm)
    ax.set_ylim(-0.2 * hs.H * mm, (hs.H + hs.tb) * mm * 1.12)
    ax.set_aspect("equal")
    ax.set_xlabel("mm")
    ax.set_title(title or f"{hs.N} fins, spacing {hs.s * mm:.1f} mm, L = {hs.L * mm:.0f} mm")
    return ax


def plot_sweep(points, ax=None, title=None, highlight=None):
    """Heat sink resistance vs number of fins, one curve per fin thickness."""
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))
    for t in sorted({p.t for p in points}):
        pts = [p for p in points if p.t == t and p.result is not None]
        ax.plot([p.N for p in pts], [p.R for p in pts], "o-", ms=3, label=f"t = {t * 1e3:g} mm")
    if highlight is not None:
        ax.plot(highlight.N, highlight.R, "*", ms=16, color="k", label=f"best: {highlight.N} fins")
        ax.set_ylim(0, 3 * highlight.R)
    ax.set_xlabel("Number of fins")
    ax.set_ylabel("Heat sink resistance [K/W]")
    ax.set_title(title or "Design sweep")
    ax.grid(alpha=0.3)
    ax.legend()
    return ax


def plot_operating_point(hs, fan: Fan, T_amb_C: float = 25.0, ax=None, title=None):
    """Fan curve and heat sink system curve; the intersection is the operating point."""
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 4.5))
    frontal = hs.W * hs.H
    q = np.linspace(fan.max_flow * 0.01, fan.max_flow, 60)
    T = T_amb_C + KELVIN
    dp = [forced_convection(hs, f / frontal, T + 20, T)[2]["total"] for f in q]
    ax.plot(np.array(fan.flow) * 1e3, fan.pressure, "C0-", lw=2, label=fan.name)
    ax.plot(q * 1e3, dp, "C3-", lw=2, label="heat sink (system curve)")
    diff = np.array([fan.dp(f) for f in q]) - np.array(dp)
    i = np.where(np.diff(np.sign(diff)))[0]
    if len(i):
        j = i[0]
        f0 = q[j] - diff[j] * (q[j + 1] - q[j]) / (diff[j + 1] - diff[j])
        ax.plot(f0 * 1e3, fan.dp(f0), "ko", ms=8)
        ax.annotate(f"  {f0 * 1e3:.1f} L/s, {fan.dp(f0):.0f} Pa", (f0 * 1e3, fan.dp(f0)))
    ax.set_xlabel("Air flow [L/s]")
    ax.set_ylabel("Static pressure [Pa]")
    ax.set_ylim(bottom=0, top=max(fan.pressure) * 1.1)
    ax.set_title(title or f"Operating point — {hs.N} fins")
    ax.grid(alpha=0.3)
    ax.legend()
    return ax


def plot_resistances(r: HeatSinkResult, ax=None, title=None):
    """Temperature rise contributed by each resistance (stacked from ambient to source)."""
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 3.2))
    parts = [("air (convection" + (" + radiation" if not math.isinf(r.R_rad) else "") + ")", r.R_sink),
             ("base conduction", r.R_base), ("spreading", r.R_spreading), ("interface", r.R_contact)]
    left = r.T_amb - KELVIN
    for i, (name, R) in enumerate(parts):
        dT = r.Q * R
        ax.barh(0, dT, left=left, height=0.5, color=f"C{i}", label=f"{name}: {dT:.1f} K")
        left += dT
    ax.set_yticks([])
    ax.set_xlabel("Temperature [°C]  (ambient → source)")
    ax.set_title(title or f"Source at {r.T_source - KELVIN:.1f} °C with {r.Q:.0f} W")
    ax.set_ylim(-0.4, 1.3)
    ax.legend(loc="upper center", ncol=2, fontsize=9, frameon=False)
    ax.set_xlim(r.T_amb - KELVIN - 2, r.T_source - KELVIN + 2)
    return ax


def plot_spacing_study(base, source, Q, T_amb_C, spacings, ax=None, title=None):
    """Natural convection: resistance vs fin spacing, with the Bar-Cohen & Rohsenow optimum."""
    from .model import analyze
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))
    R, s_act, s_opt = [], [], None
    for s in spacings:
        hs = type(base).from_spacing(base.W, base.L, base.H, base.t, s, base.tb, k=base.k,
                                     emissivity=base.emissivity)
        r = analyze(hs, source, Q, T_amb_C)
        R.append(r.R_heatsink)
        s_act.append(hs.s * 1e3)
        if r.R_heatsink <= min(R):
            s_opt = r.s_opt          # optimum evaluated at the best design's temperature
    ax.plot(s_act, R, "o-", ms=3)
    ax.axvline(s_opt * 1e3, color="k", ls="--", label=f"Bar-Cohen & Rohsenow optimum {s_opt * 1e3:.1f} mm")
    ax.set_xlabel("Fin spacing [mm]")
    ax.set_ylabel("Heat sink resistance [K/W]")
    ax.set_title(title or "Natural convection: effect of fin spacing")
    ax.grid(alpha=0.3)
    ax.legend()
    return ax
