"""Small z-buffer renderer for triangle meshes (numpy only): shaded pictures of the heat sinks."""

import numpy as np


def render(tri: np.ndarray, width: int = 1200, height: int = 900, elev: float = 28, azim: float = -50,
           color=(0.62, 0.69, 0.78), background=(1.0, 1.0, 1.0), light=(-0.6, -0.4, 0.7),
           positions: np.ndarray | None = None, scalar_fn=None, cmap="inferno", vmin=None, vmax=None,
           alpha_background: bool = False) -> np.ndarray:
    """Render triangles (n × 3 × 3, z up) to an RGB(A) image array with Lambert shading.

    With scalar_fn, each pixel is coloured by scalar_fn(p) through the colormap, where p are the
    pixel's coordinates interpolated from `positions` (same shape as tri; default: tri itself),
    e.g. a temperature field defined in the heat sink's own coordinates."""
    e, a = np.radians(elev), np.radians(azim)
    view = np.array([np.cos(e) * np.sin(a), -np.cos(e) * np.cos(a), np.sin(e)])   # towards the viewer
    right = np.cross([0, 0, 1], view)
    right /= np.linalg.norm(right)
    up = np.cross(view, right)
    pts = tri.reshape(-1, 3) - tri.reshape(-1, 3).mean(axis=0)
    u, v, depth = pts @ right, pts @ up, pts @ view
    span = max(np.ptp(u) / width, np.ptp(v) / height) * 1.08
    px = ((u - u.mean()) / span + width / 2).reshape(-1, 3)
    py = (height / 2 - (v - v.mean()) / span).reshape(-1, 3)
    pz = depth.reshape(-1, 3)

    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    nrm /= np.linalg.norm(nrm, axis=1)[:, None]
    L = np.array(light, float) / np.linalg.norm(light)
    lum = 0.35 + 0.65 * np.clip(np.abs(nrm @ L), 0, 1)
    facing = nrm @ view > 0

    img = np.ones((height, width, 3)) * background
    zbuf = np.full((height, width), -np.inf)
    pos = tri if positions is None else positions
    if scalar_fn is not None:
        import matplotlib
        cm = matplotlib.colormaps[cmap]
        sval = np.full((height, width), np.nan)
        shade = np.zeros((height, width))
    for k in np.where(facing)[0]:
        x, y, z = px[k], py[k], pz[k]
        x0, x1 = int(max(np.floor(x.min()), 0)), int(min(np.ceil(x.max()), width - 1))
        y0, y1 = int(max(np.floor(y.min()), 0)), int(min(np.ceil(y.max()), height - 1))
        if x1 < x0 or y1 < y0:
            continue
        gx, gy = np.meshgrid(np.arange(x0, x1 + 1) + 0.5, np.arange(y0, y1 + 1) + 0.5)
        d = (y[1] - y[2]) * (x[0] - x[2]) + (x[2] - x[1]) * (y[0] - y[2])
        if abs(d) < 1e-12:
            continue
        w0 = ((y[1] - y[2]) * (gx - x[2]) + (x[2] - x[1]) * (gy - y[2])) / d
        w1 = ((y[2] - y[0]) * (gx - x[2]) + (x[0] - x[2]) * (gy - y[2])) / d
        w2 = 1 - w0 - w1
        inside = (w0 >= -1e-6) & (w1 >= -1e-6) & (w2 >= -1e-6)
        zz = w0 * z[0] + w1 * z[1] + w2 * z[2]
        sub = zbuf[y0:y1 + 1, x0:x1 + 1]
        win = inside & (zz > sub)
        sub[win] = zz[win]
        if scalar_fn is None:
            img[y0:y1 + 1, x0:x1 + 1][win] = np.array(color) * lum[k]
        else:
            P = (w0[win][:, None] * pos[k, 0] + w1[win][:, None] * pos[k, 1] + w2[win][:, None] * pos[k, 2])
            sval[y0:y1 + 1, x0:x1 + 1][win] = scalar_fn(P)
            shade[y0:y1 + 1, x0:x1 + 1][win] = lum[k]
    if scalar_fn is not None:
        drawn = ~np.isnan(sval)
        lo = np.nanmin(sval) if vmin is None else vmin
        hi = np.nanmax(sval) if vmax is None else vmax
        rgb = cm(np.clip((sval[drawn] - lo) / (hi - lo), 0, 1))[:, :3]
        img[drawn] = rgb * (0.6 + 0.4 * shade[drawn])[:, None]
    if alpha_background:
        alpha = np.isfinite(zbuf).astype(float)
        return np.dstack([img, alpha])
    return img
