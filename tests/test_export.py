import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from heatsink import PlateFinHeatSink  # noqa: E402
from heatsink.export import mesh_volume, triangles, write_stl  # noqa: E402

DESIGNS = [PlateFinHeatSink(W=0.15, L=0.15, H=0.04, t=1.5e-3, N=16, tb=0.006),
           PlateFinHeatSink(W=0.08, L=0.10, H=0.04, t=0.8e-3, N=32, tb=0.008),
           PlateFinHeatSink(W=0.05, L=0.03, H=0.01, t=2e-3, N=2, tb=0.003)]


@pytest.mark.parametrize("hs", DESIGNS)
def test_volume_matches_geometry(hs):
    assert mesh_volume(triangles(hs)) == pytest.approx(hs.volume * 1e9, rel=1e-9)   # mm³


@pytest.mark.parametrize("hs", DESIGNS)
def test_watertight_and_consistently_oriented(hs):
    tri = np.round(triangles(hs), 9)
    edges = Counter()
    for t in tri:
        for a, b in ((0, 1), (1, 2), (2, 0)):
            edges[(tuple(t[a]), tuple(t[b]))] += 1
    for (a, b), n in edges.items():
        assert n == 1 and edges[(b, a)] == 1      # each directed edge once, its reverse once


def test_binary_stl_file(tmp_path):
    hs = DESIGNS[1]
    path = write_stl(hs, tmp_path / "hs.stl")
    data = path.read_bytes()
    n = int.from_bytes(data[80:84], "little")
    assert n == len(triangles(hs)) and len(data) == 84 + 50 * n


def test_box_is_closed_with_outward_normals():
    from heatsink.export import box_triangles
    tri = box_triangles(28, 28, 3)
    assert mesh_volume(tri) == pytest.approx(28 * 28 * 3)
