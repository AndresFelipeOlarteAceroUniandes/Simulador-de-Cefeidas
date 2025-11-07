
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from typing import Iterable

# Constantes físicas
SOLAR_RADIUS_M = 6.957e8           # m
DAY_S = 86400.0                    # s
SIGMA = 5.670374419e-8             # W m^-2 K^-4

@dataclass
class CepheidParams:
    """Parámetros mínimos para una Cefeida (MVP)."""
    P: float           # Periodo [días]
    R0: float          # Radio medio [R_sun]
    T0: float          # Temperatura efectiva media [K]
    v_gamma: float = 0.0        # Velocidad sistémica [km/s]
    p_factor: float = 1.36      # Factor de proyección (1.3–1.5 típico)
    # Armónicos (fracciones respecto al valor medio)
    AR: np.ndarray = field(default_factory=lambda: np.array([0.10], dtype=float))
    AT: np.ndarray = field(default_factory=lambda: np.array([0.03], dtype=float))
    phi_R: np.ndarray = field(default_factory=lambda: np.array([0.0], dtype=float))
    phi_T: np.ndarray = field(default_factory=lambda: np.array([0.6], dtype=float))
    # Fotometría
    mV_mean: float = 5.0      # magnitud promedio asumida para banda V
    A_V: float = 0.0          # extinción total en V [mag]

def _ensure_array(x: Iterable | float) -> np.ndarray:
    return np.asarray(x, dtype=float)

def angular_freq(P_days: float) -> float:
    """Frecuencia angular [rad/s] para un periodo en días."""
    return 2.0 * np.pi / (P_days * DAY_S)

def phase_to_time_s(phase: np.ndarray | float, P_days: float) -> np.ndarray:
    phase = _ensure_array(phase)
    return phase * P_days * DAY_S

def _frac_series(A: np.ndarray, phi: np.ndarray, w: float, t: np.ndarray) -> np.ndarray:
    s = np.zeros_like(t, dtype=float)
    for Ak, ph in zip(A, phi):
        s += Ak * np.sin(w * t + ph)
    return s

def radius_Rsun(params: CepheidParams, t_s: np.ndarray | float) -> np.ndarray:
    """Radio instantáneo en unidades de R_sun."""
    t = _ensure_array(t_s)
    w = angular_freq(params.P)
    f = _frac_series(params.AR, params.phi_R, w, t)
    return params.R0 * (1.0 + f)

def _radius_m(params: CepheidParams, t_s: np.ndarray | float) -> np.ndarray:
    return radius_Rsun(params, t_s) * SOLAR_RADIUS_M

def teff(params: CepheidParams, t_s: np.ndarray | float) -> np.ndarray:
    """Temperatura efectiva instantánea [K]."""
    t = _ensure_array(t_s)
    w = angular_freq(params.P)
    f = _frac_series(params.AT, params.phi_T, w, t)
    return params.T0 * (1.0 + f)

def luminosity_W(params: CepheidParams, t_s: np.ndarray | float) -> np.ndarray:
    """Luminosidad instantánea [W] via Stefan–Boltzmann."""
    R = _radius_m(params, t_s)
    T = teff(params, t_s)
    return 4.0 * np.pi * SIGMA * R**2 * T**4

def vrad_kms(params: CepheidParams, t_s: np.ndarray | float) -> np.ndarray:
    """Velocidad radial [km/s] = v_gamma + p * dR/dt proyectada."""
    t = _ensure_array(t_s)
    w = angular_freq(params.P)
    # dR/dt analítica a partir de los armónicos
    dRdt = np.zeros_like(t, dtype=float)
    for Ak, ph in zip(params.AR, params.phi_R):
        dRdt += params.R0 * SOLAR_RADIUS_M * (Ak * w * np.cos(w * t + ph))  # m/s
    return params.v_gamma + params.p_factor * (dRdt / 1_000.0)  # km/s
