"""Plate-fin heat sink design: analytical model, optimisation and CFD comparison."""
from .air import air
from .geometry import MATERIALS, HeatSource, PlateFinHeatSink
from .model import Fan, HeatSinkResult, analyze, heat_rate_forced
