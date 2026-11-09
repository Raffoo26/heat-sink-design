"""
Heat transfer and pressure drop correlations for plate-fin heat sinks.

Fins
    Straight rectangular fin, efficiency with the adiabatic-tip solution
    η = tanh(mL)/(mL), m = √(2h / k t), optionally with the corrected length
    L_c = L + t/2 to account for convection from the tip (Incropera & DeWitt).

Natural convection (vertical channels)
    Bar-Cohen & Rohsenow (1984), symmetric isothermal parallel plates:
        Nu_s = [576 / El² + 2.873 / El^½]^(−½),   El = Ra_s · s / L
    Optimum spacing  s_opt = 2.714 · L / Ra_L^¼.
    Outer surfaces: Churchill & Chu (1975) correlation for a vertical plate.

Forced convection (ducted, no tip clearance)
    Teertstra, Yovanovich & Culham (1999), simultaneously developing laminar flow
    between parallel plates (h referred to the inlet air temperature):
        Nu_b = [(Re_b* Pr / 2)^(−3) + (0.664 √Re_b* Pr^⅓ √(1 + 3.65/√Re_b*))^(−3)]^(−⅓)
        Re_b* = Re_b · s / L,   Re_b = V_ch s / ν,   Nu_b = h s / k_air

Pressure drop (Kays & London; Muzychka & Yovanovich)
    Δp = (K_c + 4 f_app L / D_h + K_e) · ½ ρ V_ch²
    K_c = 0.42 (1 − σ²),  K_e = (1 − σ²)²,  σ = s / (s + t)
    f_app Re = √[(3.44 / √L*)² + (f Re)_fd²],  L* = L / (D_h Re)
    (f Re)_fd for rectangular ducts from Shah & London (1978).

Spreading resistance
    Lee, Song, Au & Moran (1995), centred source on a plate cooled by an
    effective heat transfer coefficient (average source temperature).
"""

import math

from .air import G, AirProperties


# ---------------------------------------------------------------- fins
def fin_m(h: float, k: float, t: float) -> float:
    return math.sqrt(2.0 * h / (k * t))


def fin_efficiency(h: float, k: float, t: float, H: float, tip_convection: bool = True) -> float:
    """Efficiency of a straight rectangular fin of height H and thickness t.

    tip_convection=True uses the corrected length H + t/2 (fin tip exposed to the air);
    False assumes an adiabatic tip (e.g. fins touching a shroud)."""
    Hc = H + t / 2 if tip_convection else H
    mL = fin_m(h, k, t) * Hc
    return 1.0 if mL < 1e-9 else math.tanh(mL) / mL


def fin_efficiency_exact(h: float, k: float, t: float, H: float) -> float:
    """Exact efficiency with a convective tip (per unit fin length), for checking."""
    m = fin_m(h, k, t)
    r = h / (m * k)
    q = math.sqrt(2 * h * k * t) * (math.sinh(m * H) + r * math.cosh(m * H)) / \
        (math.cosh(m * H) + r * math.sinh(m * H))                  # per unit length and θ_b
    return q / (h * (2 * H + t))


def overall_efficiency(eta_fin: float, A_fins: float, A_total: float) -> float:
    """η_o = 1 − (A_f / A_t)(1 − η_f)."""
    return 1.0 - A_fins / A_total * (1.0 - eta_fin)


# ---------------------------------------------------------------- natural convection
def rayleigh(air: AirProperties, dT: float, length: float) -> float:
    return G * air.beta * dT * length ** 3 / (air.nu * air.alpha)


def nu_channel_natural(Ra_s: float, s: float, L: float) -> float:
    """Bar-Cohen & Rohsenow, symmetric isothermal vertical plates."""
    El = Ra_s * s / L
    return (576.0 / El ** 2 + 2.873 / math.sqrt(El)) ** -0.5


def optimum_spacing_natural(air: AirProperties, dT: float, L: float) -> float:
    """Bar-Cohen & Rohsenow optimum spacing for isothermal vertical plates [m]."""
    return 2.714 * L / rayleigh(air, dT, L) ** 0.25


def nu_vertical_plate(Ra_L: float, Pr: float) -> float:
    """Churchill & Chu, average Nusselt number of a vertical isothermal plate."""
    return (0.825 + 0.387 * Ra_L ** (1 / 6) / (1 + (0.492 / Pr) ** (9 / 16)) ** (8 / 27)) ** 2


# ---------------------------------------------------------------- forced convection
def nu_channel_forced(Re_b: float, Pr: float, s: float, L: float) -> float:
    """Teertstra, Yovanovich & Culham composite model (h based on the inlet temperature)."""
    Re_star = Re_b * s / L
    fd = Re_star * Pr / 2
    dev = 0.664 * math.sqrt(Re_star) * Pr ** (1 / 3) * math.sqrt(1 + 3.65 / math.sqrt(Re_star))
    return (fd ** -3 + dev ** -3) ** (-1 / 3)


def fRe_rectangular(alpha: float) -> float:
    """Fully developed Fanning f·Re of a rectangular duct, aspect ratio α ≤ 1 (Shah & London)."""
    a = min(alpha, 1 / alpha)
    return 24 * (1 - 1.3553 * a + 1.9467 * a ** 2 - 1.7012 * a ** 3 + 0.9564 * a ** 4 - 0.2537 * a ** 5)


def channel_pressure_drop(air: AirProperties, V_ch: float, s: float, H: float, L: float,
                          sigma: float) -> dict:
    """Pressure drop through the fin channels (fully ducted heat sink)."""
    Dh = 2 * s * H / (s + H)
    Re = V_ch * Dh / air.nu
    L_star = L / (Dh * Re)
    f_app = math.sqrt((3.44 / math.sqrt(L_star)) ** 2 + fRe_rectangular(s / H) ** 2) / Re
    q_dyn = 0.5 * air.rho * V_ch ** 2
    Kc, Ke = 0.42 * (1 - sigma ** 2), (1 - sigma ** 2) ** 2
    friction = 4 * f_app * L / Dh * q_dyn
    return {"total": (Kc + Ke) * q_dyn + friction, "friction": friction,
            "contraction": Kc * q_dyn, "expansion": Ke * q_dyn, "Re_Dh": Re, "Dh": Dh, "f_app": f_app}


# ---------------------------------------------------------------- spreading resistance
def spreading_resistance(A_source: float, A_plate: float, t_base: float, k: float,
                         R_sink: float) -> float:
    """Average spreading resistance [K/W] of a centred source (Lee et al. 1995).

    R_sink is the resistance from the far side of the base to the air (fins + convection)."""
    if A_source >= A_plate:
        return 0.0
    a, b = math.sqrt(A_source / math.pi), math.sqrt(A_plate / math.pi)
    eps, tau = a / b, t_base / b
    Bi = 1.0 / (math.pi * k * b * R_sink)
    lam = math.pi + 1.0 / (math.sqrt(math.pi) * eps)
    th = math.tanh(lam * tau)
    phi = (th + lam / Bi) / (1 + lam / Bi * th)
    psi = 0.5 * (1 - eps) ** 1.5 * phi
    return psi / (math.sqrt(math.pi) * k * a)
