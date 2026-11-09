"""
Thermal model of a plate-fin heat sink: resistance network

    source ──R_contact── base (under the source) ──R_spreading── ──R_base(1-D)── fins ──R_conv── air
                                                                                  └──R_rad─── surroundings

Natural convection: vertical fins, the air rises along L.
Forced convection: air ducted along L through the channels at a given approach velocity,
or at the operating point of a fan.
"""

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq

from . import correlations as c
from .air import KELVIN, SIGMA, air
from .geometry import HeatSource, PlateFinHeatSink


@dataclass(frozen=True)
class Fan:
    """Fan curve given by points (volumetric flow [m³/s], static pressure [Pa])."""
    flow: tuple
    pressure: tuple
    name: str = "fan"

    def __post_init__(self):
        if len(self.flow) != len(self.pressure) or len(self.flow) < 2:
            raise ValueError("The fan curve needs at least two (flow, pressure) points.")
        if any(np.diff(self.flow) <= 0) or any(np.diff(self.pressure) > 0):
            raise ValueError("Fan curve flows must increase and pressures must not increase.")

    @classmethod
    def quadratic(cls, p_max: float, flow_max: float, name: str = "fan") -> "Fan":
        """Simple curve Δp = p_max (1 − (V̇/V̇_max)²)."""
        q = np.linspace(0, flow_max, 11)
        return cls(tuple(q), tuple(p_max * (1 - (q / flow_max) ** 2)), name)

    def dp(self, flow: float) -> float:
        return float(np.interp(flow, self.flow, self.pressure, right=-1e9))

    @property
    def max_flow(self) -> float:
        return self.flow[-1]


@dataclass
class HeatSinkResult:
    hs: PlateFinHeatSink
    source: HeatSource
    mode: str                    # "natural" or "forced"
    Q: float                     # heat dissipated [W]
    T_amb: float                 # ambient / inlet air temperature [K]
    h: float                     # mean convective coefficient [W/(m²·K)]
    eta_fin: float
    eta_overall: float
    R_conv: float                # fins + base to air [K/W]
    R_rad: float                 # radiation [K/W] (inf if emissivity = 0)
    R_base: float                # 1-D conduction through the base
    R_spreading: float
    R_contact: float
    Q_rad: float                 # heat removed by radiation [W]
    velocity: float = 0.0        # approach velocity [m/s] (forced)
    channel_velocity: float = 0.0
    flow_rate: float = 0.0       # volumetric flow [m³/s] (forced)
    pressure_drop: float = 0.0   # [Pa]
    Re: float = 0.0              # channel Reynolds number (based on D_h)
    T_air_out: float = 0.0       # mean outlet air temperature [K] (forced)
    s_opt: float = 0.0           # Bar-Cohen & Rohsenow optimum spacing [m] (natural)
    warnings: list = field(default_factory=list)

    @property
    def R_sink(self) -> float:
        """Resistance from the base to the ambient: convection ∥ radiation [K/W]."""
        return 1 / (1 / self.R_conv + (0 if math.isinf(self.R_rad) else 1 / self.R_rad))

    @property
    def R_heatsink(self) -> float:
        """Heat sink resistance as in datasheets: mean base under the source → ambient [K/W]."""
        return self.R_sink + self.R_base + self.R_spreading

    @property
    def R_total(self) -> float:
        return self.R_heatsink + self.R_contact

    @property
    def T_base(self) -> float:
        """Mean temperature of the finned side of the base [K]."""
        return self.T_amb + self.Q * self.R_sink

    @property
    def T_source(self) -> float:
        """Mean temperature of the source case [K]."""
        return self.T_amb + self.Q * self.R_total

    def summary(self) -> str:
        hs = self.hs
        lines = [f"Heat sink: {hs.describe()}",
                 f"Mode: {self.mode} convection, Q = {self.Q:.1f} W, ambient {self.T_amb - KELVIN:.1f} °C"]
        if self.mode == "forced":
            lines.append(f"Air: approach {self.velocity:.2f} m/s, channels {self.channel_velocity:.2f} m/s, "
                         f"{self.flow_rate * 1e3:.2f} L/s, Δp = {self.pressure_drop:.1f} Pa, Re = {self.Re:.0f}, "
                         f"outlet {self.T_air_out - KELVIN:.1f} °C")
        else:
            lines.append(f"Optimum spacing (Bar-Cohen & Rohsenow): {self.s_opt * 1e3:.2f} mm "
                         f"(actual {hs.s * 1e3:.2f} mm)")
        lines += [f"h = {self.h:.1f} W/(m²·K), fin efficiency {self.eta_fin * 100:.1f} %, "
                  f"overall {self.eta_overall * 100:.1f} %",
                  f"Radiation: {self.Q_rad:.1f} W ({self.Q_rad / self.Q * 100:.0f} %)",
                  "Resistances [K/W]: " + ", ".join(f"{n} {v:.3f}" for n, v in [
                      ("convection", self.R_conv), ("radiation", self.R_rad), ("sink", self.R_sink),
                      ("base", self.R_base), ("spreading", self.R_spreading), ("contact", self.R_contact)]),
                  f"HEAT SINK RESISTANCE = {self.R_heatsink:.3f} K/W  (total with contact {self.R_total:.3f} K/W)",
                  f"Base {self.T_base - KELVIN:.1f} °C, source {self.T_source - KELVIN:.1f} °C"]
        lines += [f"WARNING: {w}" for w in self.warnings]
        return "\n".join(lines)


# ------------------------------------------------------------------ building blocks
def _envelope_area(hs: PlateFinHeatSink) -> float:
    """Area of the smallest box around the fins that radiates to the surroundings
    (front, two sides, top and bottom; the back of the base is on the device)."""
    Ht = hs.H + hs.tb
    return hs.W * hs.L + 2 * Ht * hs.L + 2 * hs.W * Ht


def _radiation_resistance(hs: PlateFinHeatSink, T_s: float, T_amb: float) -> float:
    if hs.emissivity == 0:
        return math.inf
    h_rad = hs.emissivity * SIGMA * (T_s ** 2 + T_amb ** 2) * (T_s + T_amb)
    return 1 / (h_rad * _envelope_area(hs))


def natural_convection_h(hs: PlateFinHeatSink, T_s: float, T_amb: float):
    """Mean h of the finned surface in natural convection, and the optimum spacing."""
    a = air((T_s + T_amb) / 2)
    dT = max(T_s - T_amb, 1e-3)
    Ra_s = c.rayleigh(a, dT, hs.s)
    h_ch = c.nu_channel_natural(Ra_s, hs.s, hs.L) * a.k / hs.s
    h_out = c.nu_vertical_plate(c.rayleigh(a, dT, hs.L), a.Pr) * a.k / hs.L
    A_in = hs.n_channels * (2 * hs.H + hs.s) * hs.L          # channel walls and floors
    A_out = (2 * hs.H + hs.N * hs.t) * hs.L                   # two outer faces and fin tips
    h = (h_ch * A_in + h_out * A_out) / (A_in + A_out)
    return h, c.optimum_spacing_natural(a, dT, hs.L)


def forced_convection(hs: PlateFinHeatSink, V_approach: float, T_s: float, T_in: float):
    """h (inlet-temperature based), channel velocity, pressure drop data (fully ducted)."""
    a = air((T_s + T_in) / 2)
    V_ch = V_approach / hs.open_area_ratio
    Re_b = V_ch * hs.s / a.nu
    h = c.nu_channel_forced(Re_b, a.Pr, hs.s, hs.L) * a.k / hs.s
    dp = c.channel_pressure_drop(a, V_ch, hs.s, hs.H, hs.L, hs.open_area_ratio)
    return h, V_ch, dp, a


def heat_rate_forced(hs: PlateFinHeatSink, T_base_C: float, T_in_C: float,
                     V_channel: float | None = None, V_approach: float | None = None) -> float:
    """Heat removed [W] by a ducted heat sink with an isothermal base (benchmark use)."""
    if (V_channel is None) == (V_approach is None):
        raise ValueError("Give either V_channel or V_approach.")
    V_app = V_approach if V_approach is not None else V_channel * hs.open_area_ratio
    T_s, T_in = T_base_C + KELVIN, T_in_C + KELVIN
    h, *_ = forced_convection(hs, V_app, T_s, T_in)
    A_f = 2 * hs.H * hs.L * hs.N
    eta = c.fin_efficiency(h, hs.k, hs.t, hs.H, tip_convection=False)
    eta_o = c.overall_efficiency(eta, A_f, A_f + hs.A_base_exposed)
    return eta_o * h * (A_f + hs.A_base_exposed) * (T_s - T_in)


# ------------------------------------------------------------------ main entry point
def analyze(hs: PlateFinHeatSink, source: HeatSource, Q: float, T_amb_C: float = 25.0,
            velocity: float | None = None, fan: Fan | None = None) -> HeatSinkResult:
    """
    Solve the heat sink for a dissipated power Q [W].

    velocity : approach air velocity [m/s] over the finned frontal area W·H (forced, ducted)
    fan      : Fan curve; the operating point is found from the heat sink pressure drop
    Neither  : natural convection with vertical fins.
    """
    if Q <= 0:
        raise ValueError("The dissipated power must be positive.")
    if velocity is not None and fan is not None:
        raise ValueError("Give either an air velocity or a fan, not both.")
    if source.a > hs.W or source.b > hs.L:
        raise ValueError("The heat source is larger than the heat sink base.")
    forced = velocity is not None or fan is not None
    T_amb = T_amb_C + KELVIN
    frontal = hs.W * hs.H

    if fan is not None:                     # operating point: fan Δp = heat sink Δp
        T_guess = T_amb + 20
        def mismatch(flow):
            *_, dp, _ = forced_convection(hs, flow / frontal, T_guess, T_amb)
            return fan.dp(flow) - dp["total"]
        if mismatch(1e-6 * fan.max_flow) <= 0:
            raise ValueError("The fan cannot push air through this heat sink.")
        velocity = brentq(mismatch, 1e-6 * fan.max_flow, fan.max_flow, xtol=1e-9) / frontal

    if forced and velocity <= 0:
        raise ValueError("The air velocity must be positive.")

    T_s = T_amb + 10.0
    for _ in range(100):
        if forced:
            h, V_ch, dp, a = forced_convection(hs, velocity, T_s, T_amb)
            A_f, tip = 2 * hs.H * hs.L * hs.N, False
            A_t = A_f + hs.A_base_exposed
        else:
            h, s_opt = natural_convection_h(hs, T_s, T_amb)
            A_f, A_t, tip = hs.A_fins, hs.A_total, True
        eta = c.fin_efficiency(h, hs.k, hs.t, hs.H, tip_convection=tip)
        eta_o = c.overall_efficiency(eta, A_f, A_t)
        R_conv = 1 / (eta_o * h * A_t)
        R_rad = _radiation_resistance(hs, T_s, T_amb)
        R_sink = 1 / (1 / R_conv + (0 if math.isinf(R_rad) else 1 / R_rad))
        T_new = T_amb + Q * R_sink
        if abs(T_new - T_s) < 1e-4:
            break
        T_s = 0.5 * (T_s + T_new)
    T_s = T_new

    R_base = hs.tb / (hs.k * hs.A_footprint)
    R_sp = c.spreading_resistance(source.area, hs.A_footprint, hs.tb, hs.k, R_sink)
    Q_rad = 0.0 if math.isinf(R_rad) else (T_s - T_amb) / R_rad
    res = HeatSinkResult(hs=hs, source=source, mode="forced" if forced else "natural", Q=Q, T_amb=T_amb,
                         h=h, eta_fin=eta, eta_overall=eta_o, R_conv=R_conv, R_rad=R_rad, R_base=R_base,
                         R_spreading=R_sp, R_contact=source.R_contact, Q_rad=Q_rad)

    w = res.warnings
    if forced:
        res.velocity, res.channel_velocity, res.flow_rate = velocity, V_ch, velocity * frontal
        res.pressure_drop, res.Re = dp["total"], dp["Re_Dh"]
        m_dot = a.rho * res.flow_rate
        res.T_air_out = T_amb + (Q - Q_rad) / (m_dot * a.cp)
        if res.Re > 2300:
            w.append(f"channel Reynolds number {res.Re:.0f} > 2300: the laminar model may be inaccurate.")
        if res.T_air_out > T_s - 1:
            w.append("the air leaves almost at the base temperature: the channels are too narrow or the flow too low.")
    else:
        res.s_opt = s_opt
        if hs.s < 0.5 * s_opt:
            w.append(f"fins far too close for natural convection (spacing {hs.s * 1e3:.1f} mm, "
                     f"optimum ≈ {s_opt * 1e3:.1f} mm): the channels choke.")
    if eta < 0.6:
        w.append(f"low fin efficiency ({eta * 100:.0f} %): the fins are too tall or too thin.")
    return res
