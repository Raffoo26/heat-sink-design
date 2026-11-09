"""
Validation of the analytical model against published results and exact solutions.

Forced convection — Jouhara & Axcell, "Modelling and simulation techniques for forced
convection heat transfer in heat sinks with rectangular fins", Simulation Modelling
Practice and Theory (2009): heat sink tested by Teertstra et al. (1999), L = 115 mm,
15 fins, H = 49 mm, t = 1.25 mm, spacing 2.18 mm, k = 200 W/(m·K), base at 50 °C,
inlet air at 20 °C. Their CFD gives 135 W at 3.75 m/s and 172 W at 7.51 m/s in the
channels; their numerical integration gives the other values.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from heatsink import HeatSource, PlateFinHeatSink, analyze, heat_rate_forced  # noqa: E402
from heatsink import correlations as c  # noqa: E402
from heatsink.air import air  # noqa: E402
from heatsink.conduction import base_plate_resistance  # noqa: E402

T, B, N = 1.25e-3, 2.18e-3, 15
TEERTSTRA = PlateFinHeatSink(W=N * T + (N - 1) * B, L=0.115, H=0.049, t=T, N=N, tb=0.005, k=200)


@pytest.mark.parametrize("u_channel, Q_cfd", [(3.75, 135), (7.51, 172)])
def test_forced_vs_published_cfd(u_channel, Q_cfd):
    assert heat_rate_forced(TEERTSTRA, 50, 20, V_channel=u_channel) == pytest.approx(Q_cfd, rel=0.06)


@pytest.mark.parametrize("u_channel, Q_num", [(1.88, 89), (3.75, 129), (5.63, 151), (7.51, 166)])
def test_forced_vs_published_numerical(u_channel, Q_num):
    # laminar range (Re_Dh < 2300); above it the flow becomes transitional
    assert heat_rate_forced(TEERTSTRA, 50, 20, V_channel=u_channel) == pytest.approx(Q_num, rel=0.06)


@pytest.mark.parametrize("h, t, H", [(10, 2e-3, 0.04), (50, 1e-3, 0.03), (150, 0.8e-3, 0.05), (5, 3e-3, 0.025)])
def test_fin_efficiency_corrected_length_vs_exact(h, t, H):
    k = 200
    assert c.fin_efficiency(h, k, t, H) == pytest.approx(c.fin_efficiency_exact(h, k, t, H), rel=0.005)


@pytest.mark.parametrize("W, tb, a, R0, k", [
    (0.10, 0.008, 0.030, 0.30, 200), (0.10, 0.004, 0.020, 0.10, 200),
    (0.06, 0.003, 0.010, 0.50, 390), (0.10, 0.006, 0.030, 0.05, 200)])
def test_spreading_vs_3d_conduction(W, tb, a, R0, k):
    """Lee et al. correlation vs a 3-D finite-volume solution of the base plate."""
    A = W * W
    fd = base_plate_resistance(W, W, tb, k, a, a, h_eff=1 / (R0 * A), n_xy=40, n_z=8)
    corr = R0 + tb / (k * A) + c.spreading_resistance(a * a, A, tb, k, R0)
    assert corr == pytest.approx(fd, rel=0.03)


def test_bar_cohen_optimum_matches_model_optimum():
    """The best spacing found by sweeping the full model is close to the BCR formula."""
    src = HeatSource(0.03, 0.03)
    best = None
    for n in range(6, 40):
        hs = PlateFinHeatSink(W=0.15, L=0.15, H=0.05, t=1.5e-3, N=n, tb=0.006)
        r = analyze(hs, src, 30, 25)
        if best is None or r.R_heatsink < best.R_heatsink:
            best = r
    assert best.hs.s == pytest.approx(best.s_opt, rel=0.15)


def test_teertstra_limits():
    a = air(300)
    # long channel → fully developed limit Re* Pr / 2 ; short channel → developing (boundary layer)
    assert c.nu_channel_forced(1.0, a.Pr, 0.002, 1.0) == pytest.approx(1.0 * 0.002 / 1.0 * a.Pr / 2, rel=0.01)
    Re_star = 5000 * 0.002 / 0.01
    dev = 0.664 * Re_star ** 0.5 * a.Pr ** (1 / 3) * (1 + 3.65 / Re_star ** 0.5) ** 0.5
    assert c.nu_channel_forced(5000, a.Pr, 0.002, 0.01) == pytest.approx(dev, rel=0.01)
