import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from heatsink import Fan, HeatSource, PlateFinHeatSink, analyze  # noqa: E402
from heatsink.air import KELVIN  # noqa: E402
from heatsink.optimize import best_design, sweep_fins  # noqa: E402

HS = PlateFinHeatSink(W=0.08, L=0.10, H=0.04, t=1e-3, N=25, tb=0.008)
SRC = HeatSource(0.0375, 0.0375, R_contact=0.03)
FAN = Fan.quadratic(60, 0.025)


def test_energy_balance_forced():
    r = analyze(HS, SRC, 100, 30, velocity=3.0)
    m_cp = r.flow_rate * 1.15 * 1007
    assert (r.T_air_out - r.T_amb) * m_cp == pytest.approx(100, rel=0.05)
    assert r.T_source > r.T_base > r.T_air_out > r.T_amb


def test_resistance_does_not_depend_on_power_in_forced_flow():
    r1, r2 = analyze(HS, SRC, 50, 30, velocity=3), analyze(HS, SRC, 150, 30, velocity=3)
    assert r1.R_heatsink == pytest.approx(r2.R_heatsink, rel=0.03)


def test_natural_convection_improves_with_temperature_difference():
    hs = PlateFinHeatSink(W=0.15, L=0.15, H=0.04, t=2e-3, N=15, tb=0.006)
    src = HeatSource(0.03, 0.03)
    assert analyze(hs, src, 60, 25).R_heatsink < analyze(hs, src, 10, 25).R_heatsink


def test_radiation_helps():
    bare = PlateFinHeatSink(W=0.15, L=0.15, H=0.04, t=2e-3, N=15, tb=0.006, emissivity=0.05)
    black = PlateFinHeatSink(W=0.15, L=0.15, H=0.04, t=2e-3, N=15, tb=0.006, emissivity=0.85)
    src = HeatSource(0.03, 0.03)
    rb, rk = analyze(bare, src, 40, 25), analyze(black, src, 40, 25)
    assert rk.R_heatsink < rb.R_heatsink and rk.Q_rad > 4 * rb.Q_rad


def test_fan_operating_point_balances_pressures():
    r = analyze(HS, SRC, 125, 35, fan=FAN)
    assert FAN.dp(r.flow_rate) == pytest.approx(r.pressure_drop, rel=0.02)


def test_more_fins_choke_the_fan():
    flows = [analyze(PlateFinHeatSink(W=0.08, L=0.10, H=0.04, t=1e-3, N=n, tb=0.008), SRC, 100, 30,
                     fan=FAN).flow_rate for n in (10, 20, 30)]
    assert flows[0] > flows[1] > flows[2]


def test_fan_sweep_has_an_interior_optimum():
    pts = sweep_fins(HS, SRC, 125, 35, range(6, 45), fan=FAN)
    best = best_design(pts)
    assert 6 < best.N < 44


def test_thicker_base_reduces_spreading():
    thin = analyze(PlateFinHeatSink(W=0.1, L=0.1, H=0.03, t=1e-3, N=20, tb=0.003), SRC, 100, 30, velocity=3)
    thick = analyze(PlateFinHeatSink(W=0.1, L=0.1, H=0.03, t=1e-3, N=20, tb=0.010), SRC, 100, 30, velocity=3)
    assert thick.R_spreading < thin.R_spreading


def test_temperatures_in_kelvin_and_summary():
    r = analyze(HS, SRC, 100, 30, velocity=3)
    assert r.T_amb == pytest.approx(30 + KELVIN)
    assert "HEAT SINK RESISTANCE" in r.summary()


@pytest.mark.parametrize("kwargs", [dict(N=1), dict(N=100), dict(t=-1e-3), dict(emissivity=2)])
def test_invalid_geometry(kwargs):
    base = dict(W=0.08, L=0.10, H=0.04, t=1e-3, N=25, tb=0.008)
    base.update(kwargs)
    with pytest.raises(ValueError):
        PlateFinHeatSink(**base)


def test_invalid_inputs():
    with pytest.raises(ValueError):
        analyze(HS, SRC, -5, 30)
    with pytest.raises(ValueError):
        analyze(HS, HeatSource(0.2, 0.2), 50, 30)
    with pytest.raises(ValueError):
        analyze(HS, SRC, 50, 30, velocity=2, fan=FAN)
    with pytest.raises(ValueError):
        Fan((0, 1), (10, 20))
