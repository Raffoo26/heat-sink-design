"""Case study 1: passive (natural convection) heat sink for a 40 W COB LED."""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from heatsink import HeatSource, PlateFinHeatSink, analyze  # noqa: E402
from heatsink.air import KELVIN  # noqa: E402
from heatsink.optimize import best_design, sweep_fins  # noqa: E402
from heatsink.plots import plot_profile, plot_resistances, plot_spacing_study, plot_sweep  # noqa: E402

Q, T_AMB, T_MAX = 40.0, 35.0, 80.0          # heat to dissipate [W], ambient, max LED case temperature [°C]
LED = HeatSource(0.028, 0.028, R_contact=0.10)   # 28 × 28 mm COB substrate, thermal paste
BASE = PlateFinHeatSink(W=0.15, L=0.15, H=0.04, t=2e-3, N=12, tb=0.006, emissivity=0.85)

print(f"Requirement: R_total ≤ ({T_MAX:.0f} − {T_AMB:.0f}) / {Q:.0f} = {(T_MAX - T_AMB) / Q:.3f} K/W\n")
points = sweep_fins(BASE, LED, Q, T_AMB, range(4, 31), [1.5e-3, 2e-3, 3e-3])
best = best_design(points)
r = best.result
print(r.summary())

bare = analyze(PlateFinHeatSink(**{**r.hs.__dict__, "emissivity": 0.05}), LED, Q, T_AMB)
print(f"\nSame heat sink, bare aluminium (ε = 0.05): LED at {bare.T_source - KELVIN:.1f} °C "
      f"instead of {r.T_source - KELVIN:.1f} °C")

fig, ax = plt.subplots(2, 2, figsize=(14, 9))
plot_sweep(points, ax=ax[0, 0], highlight=best, title="40 W LED, 150 × 150 mm, fins 40 mm, black anodised")
plot_spacing_study(PlateFinHeatSink(**{**r.hs.__dict__}), LED, Q, T_AMB, np.linspace(3e-3, 25e-3, 40), ax=ax[0, 1])
plot_profile(r.hs, LED, ax=ax[1, 0])
plot_resistances(r, ax=ax[1, 1])
fig.tight_layout()
fig.savefig(ROOT / "images" / "led_natural.png", dpi=130)
print("saved images/led_natural.png")
