"""Export the two optimised heat sinks as STL files (mm) for 3D printing or CAD, and draw a 3-D view."""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from heatsink import PlateFinHeatSink  # noqa: E402
from heatsink.export import triangles, write_box_stl, write_stl  # noqa: E402
from heatsink.render import render  # noqa: E402




DESIGNS = {
    "led_natural_heatsink": ("Case 1 — 40 W LED, natural convection",
                             PlateFinHeatSink(W=0.15, L=0.15, H=0.04, t=1.5e-3, N=16, tb=0.006)),
    "cpu_fan_heatsink": ("Case 2 — 125 W CPU with fan",
                         PlateFinHeatSink(W=0.08, L=0.10, H=0.04, t=0.8e-3, N=32, tb=0.008)),
}
SOURCES = {"led_source": (28, 28, 3), "cpu_source": (37.5, 37.5, 3)}   # heat sources [mm], other colour

for name, (a, b, c) in SOURCES.items():
    path = write_box_stl(a, b, c, ROOT / "cad" / f"{name}.stl", name=name)
    print(f"{path.relative_to(ROOT)}: {a:g} × {b:g} × {c:g} mm")

fig = plt.figure(figsize=(15, 6.6))
for i, (name, (title, hs)) in enumerate(DESIGNS.items(), start=1):
    path = write_stl(hs, ROOT / "cad" / f"{name}.stl", name=name)
    print(f"{path.relative_to(ROOT)}: {hs.describe()}, {hs.mass() * 1e3:.0f} g in aluminium")
    tri = triangles(hs)
    tri = np.stack([tri[..., 0], -tri[..., 2], tri[..., 1]], axis=-1)   # rotate: fins up (z), length along −y
    ax = fig.add_subplot(1, 2, i)
    ax.imshow(render(tri, 1300, 950, elev=30, azim=-35))
    ax.axis("off")
    ax.set_title(f"{title}\n{hs.W * 1e3:.0f} × {hs.L * 1e3:.0f} mm, {hs.N} fins {hs.H * 1e3:.0f} × {hs.t * 1e3:g} mm, "
                 f"spacing {hs.s * 1e3:.2f} mm")
fig.tight_layout(pad=1.5)
fig.savefig(ROOT / "images" / "heatsinks_3d.png", dpi=130)
print("saved images/heatsinks_3d.png")
