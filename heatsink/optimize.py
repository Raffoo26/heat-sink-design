"""Design studies: sweep the number and thickness of fins and pick the best heat sink."""

from dataclasses import dataclass, replace

from .geometry import HeatSource, PlateFinHeatSink
from .model import Fan, HeatSinkResult, analyze


@dataclass
class DesignPoint:
    N: int
    t: float
    result: HeatSinkResult | None
    error: str = ""

    @property
    def R(self) -> float:
        return self.result.R_heatsink if self.result else float("inf")


def sweep_fins(base: PlateFinHeatSink, source: HeatSource, Q: float, T_amb_C: float,
               N_values, t_values=None, velocity: float | None = None,
               fan: Fan | None = None) -> list[DesignPoint]:
    """Analyse every combination of fin count and thickness (other dimensions from `base`)."""
    points = []
    for t in (t_values or [base.t]):
        for N in N_values:
            try:
                hs = replace(base, N=int(N), t=t)
                points.append(DesignPoint(int(N), t, analyze(hs, source, Q, T_amb_C, velocity, fan)))
            except ValueError as e:
                points.append(DesignPoint(int(N), t, None, str(e)))
    return points


def best_design(points: list[DesignPoint]) -> DesignPoint:
    """Design with the lowest heat sink resistance."""
    valid = [p for p in points if p.result is not None]
    if not valid:
        raise ValueError("No feasible design in the sweep.")
    return min(valid, key=lambda p: p.R)
