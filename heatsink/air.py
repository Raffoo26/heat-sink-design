"""Air properties at atmospheric pressure (CoolProp, real-gas reference data)."""

from dataclasses import dataclass
from functools import lru_cache

from CoolProp.CoolProp import PropsSI

KELVIN = 273.15
G = 9.80665          # gravitational acceleration [m/s²]
SIGMA = 5.670374e-8  # Stefan–Boltzmann constant [W/(m²·K⁴)]


@dataclass(frozen=True)
class AirProperties:
    T: float      # temperature [K]
    rho: float    # density [kg/m³]
    mu: float     # dynamic viscosity [Pa·s]
    k: float      # thermal conductivity [W/(m·K)]
    cp: float     # specific heat [J/(kg·K)]

    @property
    def nu(self) -> float:
        """Kinematic viscosity [m²/s]."""
        return self.mu / self.rho

    @property
    def alpha(self) -> float:
        """Thermal diffusivity [m²/s]."""
        return self.k / (self.rho * self.cp)

    @property
    def Pr(self) -> float:
        return self.mu * self.cp / self.k

    @property
    def beta(self) -> float:
        """Volumetric expansion coefficient of an ideal gas, 1/T [1/K]."""
        return 1.0 / self.T


@lru_cache(maxsize=4096)
def _props(T_rounded: float, p: float) -> AirProperties:
    args = ("T", T_rounded, "P", p, "Air")
    return AirProperties(T=T_rounded, rho=PropsSI("D", *args), mu=PropsSI("V", *args),
                         k=PropsSI("L", *args), cp=PropsSI("C", *args))


def air(T_K: float, p_Pa: float = 101325.0) -> AirProperties:
    """Air properties at temperature T [K] and pressure p [Pa]."""
    if not 200 < T_K < 600:
        raise ValueError(f"Air temperature {T_K - KELVIN:.0f} °C is outside the model range (−73…327 °C).")
    return _props(round(T_K, 2), p_Pa)
