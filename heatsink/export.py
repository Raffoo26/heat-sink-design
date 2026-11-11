"""
Export a plate-fin heat sink as a watertight STL solid (for 3D printing, CAD or CFD meshing).

The cross-section (base + fins) is extruded along the fin length L. The two end caps are
split into rectangles that share all their vertices, so the mesh has no T-junctions:
every edge belongs to exactly two triangles.
"""

import struct
from pathlib import Path

import numpy as np

from .geometry import PlateFinHeatSink


def _profile(hs: PlateFinHeatSink, scale: float):
    """Outline of the cross-section (counter-clockwise, x across the fins, y up) and cap rectangles."""
    W, H, t, tb, p = (v * scale for v in (hs.W, hs.H, hs.t, hs.tb, hs.pitch))
    fins = [(i * p, i * p + t) for i in range(hs.N)]
    fins[-1] = (W - t, W)                      # last fin flush with the edge (removes rounding)
    xs = sorted({0.0, W, *[x for f in fins for x in f]})

    top = tb + H
    outline = [(x, 0.0) for x in xs] + [(W, tb)]   # bottom split like the caps: no T-junctions
    for i in range(hs.N - 1, -1, -1):          # along the top, from the right edge to the left edge
        xl, xr = fins[i]
        outline += [(xr, top), (xl, top)]
        if i > 0:
            outline += [(xl, tb), (fins[i - 1][1], tb)]
    outline.append((0.0, tb))
    rects = [(xs[i], 0.0, xs[i + 1], tb) for i in range(len(xs) - 1)]
    rects += [(xl, tb, xr, tb + H) for xl, xr in fins]
    return outline, rects


def triangles(hs: PlateFinHeatSink, scale: float = 1e3) -> np.ndarray:
    """Triangles (n × 3 × 3) of the solid, in mm by default (scale = 1000)."""
    outline, rects = _profile(hs, scale)
    L = hs.L * scale
    tris = []
    # end caps: z = 0 (normal −z) and z = L (normal +z)
    for x0, y0, x1, y1 in rects:
        a, b, c, d = (x0, y0), (x1, y0), (x1, y1), (x0, y1)
        for z, flip in ((0.0, True), (L, False)):
            q = [(*a, z), (*b, z), (*c, z), (*d, z)]
            t1, t2 = [q[0], q[1], q[2]], [q[0], q[2], q[3]]
            tris += [t1[::-1], t2[::-1]] if flip else [t1, t2]
    # side walls: each outline edge extruded along z (outline is counter-clockwise → outward normals)
    n = len(outline)
    for i in range(n):
        (x0, y0), (x1, y1) = outline[i], outline[(i + 1) % n]
        p0, p1, p2, p3 = (x0, y0, 0.0), (x1, y1, 0.0), (x1, y1, L), (x0, y0, L)
        tris += [[p0, p1, p2], [p0, p2, p3]]
    return np.array(tris, dtype=float)


def write_stl(hs: PlateFinHeatSink, path, scale: float = 1e3, name: str = "heatsink") -> Path:
    """Write a binary STL (millimetres by default) and return its path."""
    tri = triangles(hs, scale)
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        f.write(name.encode()[:80].ljust(80, b" "))
        f.write(struct.pack("<I", len(tri)))
        for nrm, t in zip(normals, tri):
            f.write(struct.pack("<12fH", *nrm, *t.ravel(), 0))
    return path


def mesh_volume(tri: np.ndarray) -> float:
    """Signed volume of a closed triangle mesh (divergence theorem)."""
    return float(np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() / 6)


def box_triangles(a: float, b: float, c: float) -> np.ndarray:
    """Closed box a × b × c (same units as given), corner at the origin."""
    v = np.array([[x, y, z] for z in (0, c) for y in (0, b) for x in (0, a)], float)
    faces = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5)]
    tris = []
    for p, q, r, s in faces:
        tris += [[v[p], v[q], v[r]], [v[p], v[r], v[s]]]
    return np.array(tris)


def write_box_stl(a: float, b: float, c: float, path, name: str = "box") -> Path:
    """Binary STL of a box with sides a, b, c [mm] (e.g. the heat source, to print in another colour)."""
    tri = box_triangles(a, b, c)
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        f.write(name.encode()[:80].ljust(80, b" "))
        f.write(struct.pack("<I", len(tri)))
        for nrm, t in zip(normals, tri):
            f.write(struct.pack("<12fH", *nrm, *t.ravel(), 0))
    return path
