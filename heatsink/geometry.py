"""
Geometry of a plate-fin heat sink and of the heat source attached to it.

                 ┌┐   ┌┐   ┌┐   ┌┐      ▲
                 ││   ││   ││   ││      │ H  fin height
                 ││ s ││   ││   ││      │
                 ││<->││   ││   ││      ▼
              ┌──┴┴───┴┴───┴┴───┴┴──┐   ▲ tb base thickness
              └─────────────────────┘   ▼
              <──────── W ──────────>     L = fin length (into the page):
                                          flow direction (forced convection)
                                          or vertical direction (natural convection)
"""

from dataclasses import dataclass

MATERIALS = {               # thermal conductivity [W/(m·K)]
    "aluminium 6063": 209.0,
    "aluminium 6061": 167.0,
    "aluminium 1050": 229.0,
    "copper": 391.0,
}


@dataclass(frozen=True)
class PlateFinHeatSink:
    W: float                  # base width, across the fins [m]
    L: float                  # fin length, along the channels [m]
    H: float                  # fin height above the base [m]
    t: float                  # fin thickness [m]
    N: int                    # number of fins
    tb: float                 # base thickness [m]
    k: float = MATERIALS["aluminium 6063"]   # conductivity of fins and base [W/(m·K)]
    emissivity: float = 0.0   # surface emissivity (≈0.05 bare, ≈0.85 black anodised)

    def __post_init__(self):
        if min(self.W, self.L, self.H, self.t, self.tb) <= 0:
            raise ValueError("All dimensions must be positive.")
        if self.N < 2:
            raise ValueError("A plate-fin heat sink needs at least 2 fins.")
        if self.N * self.t >= self.W:
            raise ValueError(f"{self.N} fins {self.t * 1e3:g} mm thick do not fit in a "
                             f"{self.W * 1e3:g} mm wide base.")
        if self.k <= 0:
            raise ValueError("Thermal conductivity must be positive.")
        if not 0 <= self.emissivity <= 1:
            raise ValueError("Emissivity must be between 0 and 1.")

    @classmethod
    def from_spacing(cls, W, L, H, t, s, tb, **kw) -> "PlateFinHeatSink":
        """Build the heat sink with as many fins as fit in W for a target spacing s."""
        N = max(2, int(round((W + s) / (s + t))))
        return cls(W=W, L=L, H=H, t=t, N=N, tb=tb, **kw)

    @property
    def s(self) -> float:
        """Clear spacing between fins (channel width) [m]."""
        return (self.W - self.N * self.t) / (self.N - 1)

    @property
    def pitch(self) -> float:
        return self.s + self.t

    @property
    def n_channels(self) -> int:
        return self.N - 1

    @property
    def A_fins(self) -> float:
        """Wetted area of all fins: both faces and the tip [m²]."""
        return self.N * (2 * self.H + self.t) * self.L

    @property
    def A_base_exposed(self) -> float:
        """Base area between the fins [m²]."""
        return self.n_channels * self.s * self.L

    @property
    def A_total(self) -> float:
        return self.A_fins + self.A_base_exposed

    @property
    def A_footprint(self) -> float:
        return self.W * self.L

    @property
    def open_area_ratio(self) -> float:
        """σ = free flow area / frontal area of the finned region."""
        return self.s / self.pitch

    @property
    def volume(self) -> float:
        """Metal volume [m³]."""
        return self.W * self.L * self.tb + self.N * self.t * self.H * self.L

    def mass(self, density: float = 2700.0) -> float:
        """Mass [kg] (default density: aluminium)."""
        return self.volume * density

    def describe(self) -> str:
        return (f"{self.W * 1e3:.0f} × {self.L * 1e3:.0f} mm base, {self.N} fins "
                f"{self.H * 1e3:.0f} × {self.t * 1e3:.1f} mm, spacing {self.s * 1e3:.2f} mm, "
                f"base {self.tb * 1e3:.1f} mm, k = {self.k:.0f} W/(m·K)")


@dataclass(frozen=True)
class HeatSource:
    """Rectangular heat source (chip, LED, power module) centred on the base."""
    a: float          # side along W [m]
    b: float          # side along L [m]
    R_contact: float = 0.0   # interface (thermal paste / pad) resistance [K/W]

    @property
    def area(self) -> float:
        return self.a * self.b
