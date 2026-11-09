"""Case study 2: ducted heat sink with a fan for a 125 W processor."""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from heatsink import Fan, HeatSource, PlateFinHeatSink, analyze  # noqa: E402
from heatsink.optimize import best_design, sweep_fins  # noqa: E402
from heatsink.plots import plot_operating_point, plot_profile, plot_resistances, plot_sweep  # noqa: E402

Q, T_AMB, T_MAX = 125.0, 35.0, 70.0              # CPU power [W], air inside the case, max case temperature
CPU = HeatSource(0.0375, 0.0375, R_contact=0.03)  # 37.5 mm heat spreader, thermal paste
FAN = Fan.quadratic(p_max=60, flow_max=0.025, name="fan: 60 Pa max, 25 L/s max")
BASE = PlateFinHeatSink(W=0.08, L=0.10, H=0.04, t=1e-3, N=20, tb=0.008)

print(f"Requirement: R_total ≤ ({T_MAX:.0f} − {T_AMB:.0f}) / {Q:.0f} = {(T_MAX - T_AMB) / Q:.3f} K/W\n")
points = sweep_fins(BASE, CPU, Q, T_AMB, range(6, 41), [0.8e-3, 1e-3, 1.5e-3], fan=FAN)
feasible = [p for p in points if p.result and p.result.hs.s >= 1.5e-3]   # extrusion / cleaning limit
best = best_design(feasible)
r = best.result
print(r.summary())

few = analyze(PlateFinHeatSink(**{**r.hs.__dict__, "N": 12}), CPU, Q, T_AMB, fan=FAN)
print(f"\nWith only 12 fins: {few.flow_rate * 1e3:.1f} L/s but R = {few.R_heatsink:.3f} K/W")

fig, ax = plt.subplots(2, 2, figsize=(14, 9))
plot_sweep(points, ax=ax[0, 0], highlight=best, title="125 W CPU, 80 × 100 mm, fins 40 mm, with fan")
plot_operating_point(r.hs, FAN, T_AMB, ax=ax[0, 1])
plot_profile(r.hs, CPU, ax=ax[1, 0])
plot_resistances(r, ax=ax[1, 1])
fig.tight_layout()
fig.savefig(ROOT / "images" / "cpu_fan.png", dpi=130)
print("saved images/cpu_fan.png")
