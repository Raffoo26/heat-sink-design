"""Validation plots: forced convection benchmark, spreading resistance, natural convection optimum."""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from heatsink import PlateFinHeatSink, heat_rate_forced  # noqa: E402
from heatsink import correlations as c  # noqa: E402
from heatsink.air import air  # noqa: E402
from heatsink.conduction import base_plate_resistance  # noqa: E402

fig, axes = plt.subplots(1, 3, figsize=(17, 5))

# 1. Forced convection: Teertstra heat sink, Jouhara & Axcell benchmark
t, b, N = 1.25e-3, 2.18e-3, 15
hs = PlateFinHeatSink(W=N * t + (N - 1) * b, L=0.115, H=0.049, t=t, N=N, tb=0.005, k=200)
u_num = [1.88, 3.75, 5.63, 7.51, 9.39, 11.3, 13.1, 15.0]
Q_num = [89, 129, 151, 166, 178, 187, 196, 203]
u = np.linspace(1.5, 15.5, 60)
ax = axes[0]
ax.plot(u, [heat_rate_forced(hs, 50, 20, V_channel=v) for v in u], "C0-", lw=2, label="this model")
ax.plot(u_num, Q_num, "s", color="C1", label="numerical (Jouhara & Axcell)")
ax.plot([3.75, 7.51], [135, 172], "o", color="C3", ms=9, label="CFD (Jouhara & Axcell)")
Re_crit = 2300 * air(308).nu / (2 * b)
ax.axvspan(Re_crit, 16, color="0.9", label="Re > 2300 (transitional)")
ax.set_xlabel("Mean channel velocity [m/s]")
ax.set_ylabel("Heat removed [W]")
ax.set_title("Forced convection: Teertstra heat sink\nbase 50 °C, inlet 20 °C")
ax.grid(alpha=0.3)
ax.legend(fontsize=9)

# 2. Spreading resistance vs 3-D conduction
ax = axes[1]
cases = [(0.10, 0.008, 0.030, 0.30, 200), (0.10, 0.004, 0.020, 0.10, 200), (0.10, 0.010, 0.040, 0.50, 200),
         (0.06, 0.003, 0.010, 0.50, 390), (0.10, 0.006, 0.030, 0.05, 200), (0.08, 0.005, 0.015, 0.20, 200)]
fd, corr = [], []
for W, tb, a, R0, k in cases:
    A = W * W
    fd.append(base_plate_resistance(W, W, tb, k, a, a, 1 / (R0 * A), n_xy=50, n_z=10))
    corr.append(R0 + tb / (k * A) + c.spreading_resistance(a * a, A, tb, k, R0))
lim = [0, max(fd) * 1.1]
ax.plot(lim, lim, "k--", lw=1)
ax.plot(fd, corr, "o", ms=8)
for f, r in zip(fd, corr):
    ax.annotate(f"{(r / f - 1) * 100:+.1f} %", (f, r), xytext=(6, -10), textcoords="offset points", fontsize=9)
ax.set_xlabel("3-D finite-volume solution [K/W]")
ax.set_ylabel("Correlation (Lee et al.) [K/W]")
ax.set_title("Base: spreading + conduction + sink resistance")
ax.grid(alpha=0.3)

# 3. Natural convection optimum spacing vs Yazıcıoğlu & Yüncu experiments
ax = axes[2]
dT = np.array([50, 75, 100, 125])
exp = {0.25: [11.0, 10.9, 10.8, 10.7], 0.34: [11.9, 11.8, 11.7, 11.6]}   # H = 25 mm
for i, (L, s_exp) in enumerate(exp.items()):
    s_bcr = [c.optimum_spacing_natural(air(298 + d / 2), d, L) * 1e3 for d in dT]
    ax.plot(dT, s_exp, "s", color=f"C{i}", ms=8, label=f"experiment, L = {L * 1e3:.0f} mm")
    ax.plot(dT, s_bcr, "-", color=f"C{i}", label=f"Bar-Cohen & Rohsenow, L = {L * 1e3:.0f} mm")
ax.set_xlabel("Base – ambient temperature difference [K]")
ax.set_ylabel("Optimum fin spacing [mm]")
ax.set_ylim(0, 14)
ax.set_title("Natural convection: short fins (H = 25 mm)\nexperiments: Yazıcıoğlu & Yüncu (2007)")
ax.grid(alpha=0.3)
ax.legend(fontsize=9)

fig.tight_layout()
fig.savefig(ROOT / "images" / "validation.png", dpi=130)
print("saved images/validation.png")
