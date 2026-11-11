"""Banner for the README: the CPU heat sink coloured by the temperature computed by the model."""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from heatsink import Fan, HeatSource, PlateFinHeatSink, analyze  # noqa: E402
from heatsink.export import triangles  # noqa: E402
from heatsink.render import render  # noqa: E402
from heatsink.thermal_map import temperature_function  # noqa: E402

hs = PlateFinHeatSink(W=0.08, L=0.10, H=0.04, t=0.8e-3, N=32, tb=0.008)
r = analyze(hs, HeatSource(0.0375, 0.0375, R_contact=0.03), 125, 35, fan=Fan.quadratic(60, 0.025))
T = temperature_function(r)

tri = triangles(hs)                                                       # mm, STL coordinates
view = np.stack([tri[..., 0], -tri[..., 2], tri[..., 1]], axis=-1)        # fins up for the picture
gx, gz, gy = np.meshgrid(np.linspace(0, hs.W * 1e3, 41), np.linspace(0, hs.L * 1e3, 41),
                         [0, hs.tb * 1e3, (hs.tb + hs.H) * 1e3])
probe = T(np.c_[gx.ravel(), gy.ravel(), gz.ravel()])
vmin, vmax = np.floor(probe.min()), np.ceil(probe.max())
img = render(view, 1500, 1000, elev=26, azim=-38, positions=tri, scalar_fn=T, cmap="turbo",
             vmin=vmin, vmax=vmax, background=(0.051, 0.067, 0.09), alpha_background=True)

BG = "#0d1117"
fig = plt.figure(figsize=(16, 5.6), facecolor=BG)
ax_img = fig.add_axes([0.47, 0.02, 0.46, 0.96])
ax_img.imshow(img)
ax_img.axis("off")

fig.text(0.04, 0.78, "Heat Sink Design", color="white", fontsize=40, fontweight="bold")
fig.text(0.04, 0.665, "Plate-fin heat sinks for electronics, from equations to a printed part",
         color="#c9d1d9", fontsize=15)
steps = ["Analytical model", "Validation", "Optimisation", "CFD", "3D print"]
x = 0.04
for i, s in enumerate(steps):
    t = fig.text(x, 0.50, f" {s} ", color="white", fontsize=12.5, fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.45", fc="#1f6feb" if i < 3 else "#30363d", ec="none"))
    fig.canvas.draw()
    x += t.get_window_extent().width / fig.bbox.width + 0.018
fig.text(0.04, 0.30, f"125 W processor · {hs.N} fins · fan at {r.flow_rate * 1e3:.1f} L/s",
         color="#8b949e", fontsize=12.5)
fig.text(0.04, 0.23, f"Heat sink resistance {r.R_heatsink:.3f} K/W · CPU at {r.T_source - 273.15:.0f} °C",
         color="#8b949e", fontsize=12.5)
fig.text(0.04, 0.08, "Colours: temperature from the model (base: 3-D conduction, fins: fin equation)",
         color="#6e7681", fontsize=10)

cax = fig.add_axes([0.945, 0.2, 0.012, 0.6])
cb = fig.colorbar(plt.cm.ScalarMappable(Normalize(vmin, vmax), "turbo"), cax=cax)
cb.set_label("°C", color="#c9d1d9")
cb.ax.tick_params(colors="#c9d1d9")
cb.outline.set_visible(False)

fig.savefig(ROOT / "images" / "banner.png", dpi=110, facecolor=BG)
print(f"saved images/banner.png  (T from {probe.min():.1f} to {probe.max():.1f} °C)")
