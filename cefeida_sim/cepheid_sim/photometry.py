from __future__ import annotations
import numpy as np
from typing import Optional, Tuple

PI = np.pi
# Constantes físicas
H  = 6.62607015e-34     # J*s
C  = 2.99792458e8       # m/s
KB = 1.380649e-23       # J/K

# Atajos Planck
C1 = 2.0 * H * C**2     #W*m^2
C2 = (H * C) / KB       # m*K

# ------------------------------
# 1) Filtro top-hat para banda V
# ------------------------------
def top_hat_filter(center_nm: float = 550.0, width_nm: float = 88.0, n: int = 400):
    """Ventana top-hat simple para banda V (aprox)."""
    lam0 = center_nm * 1e-9
    half = (width_nm * 1e-9) / 2.0
    lam = np.linspace(lam0 - half, lam0 + half, n)
    tx  = np.ones_like(lam)
    return lam, tx

# ------------------------------------------
# 2) Ley de Planck B_lambda (W m^-3 sr^-1)
# ------------------------------------------
def planck_lambda_w_m3_sr(T: np.ndarray | float, lam_m: np.ndarray) -> np.ndarray:
    """
    B_lambda(T, λ) en W m^-3 sr^-1. Broadcasting en T y λ.
    """
    T   = np.asarray(T, dtype=float)
    lam = np.asarray(lam_m, dtype=float)

    # x = (hc)/(λ k_B T)
    x = C2 / (lam[None, ...] * T[..., None])  # (..., nlam)
    x = np.clip(x, 1e-6, 100.0)               # estabilidad numérica

    # B_lambda = (2 h c^2)/λ^5 * 1/(exp(x)-1)
    B = C1 / (lam[None, ...]**5 * (np.exp(x) - 1.0))
    return B  # W m^-3 sr^-1

# -----------------------------------------------------
# 3) Modelos de flujo relativo en banda V a partir de R,T
# -----------------------------------------------------
def flux_V_rel(R_rel: np.ndarray, T_rel: np.ndarray, beta: float = 1.6) -> np.ndarray:
    """
    Modelo rápido: F_rel ∝ R_rel^2 * T_rel^beta. Útil para calibrar amplitud de V.
    """
    R_rel = np.asarray(R_rel, dtype=float)
    T_rel = np.asarray(T_rel, dtype=float)
    F_rel = R_rel**2 * T_rel**beta
    return np.clip(F_rel, 1e-12, 1e12)

def flux_band_rel_from_RT(
    R_rel: np.ndarray,
    T_K:   np.ndarray,
    center_nm: float = 550.0,
    width_nm:  float = 88.0,
    n:         int   = 400
) -> np.ndarray:
    """
    Modelo BB: integra B_lambda en una top-hat centrada en V y multiplica por R_rel^2.
    Devuelve flujo *relativo* (normalizado por su media).
    """
    R_rel = np.asarray(R_rel, dtype=float)
    T_K   = np.asarray(T_K,   dtype=float)

    lam, tx = top_hat_filter(center_nm, width_nm, n)
    B = planck_lambda_w_m3_sr(T_K, lam)    # shape (..., nlam)
    F_lambda_surface = PI * B              # integra hemisferio: F_λ ≈ π B_λ

    # Integramos en λ
    F_band_surface = np.trapz(F_lambda_surface * tx[None, :], lam, axis=-1)  # shape (...)

    # Factor geométrico relativo
    F_rel = (R_rel**2) * F_band_surface
    # Normaliza a flujo relativo
    F_rel = F_rel / float(np.mean(F_rel))
    return np.clip(F_rel, 1e-12, 1e12)

# ---------------------------------------
# 4) Magnitudes desde flujo relativo
# ---------------------------------------
def magnitude_from_flux_rel(F_rel: np.ndarray, mV_mean: float, A_V: float = 0.0) -> np.ndarray:
    """
    m_V = mV_mean + A_V + [ -2.5 log10(F_rel) - < -2.5 log10(F_rel) > ].
    Siempre centrada a <m_V>=mV_mean + A_V.
    """
    F_rel = np.asarray(F_rel, dtype=float)
    dm = -2.5 * np.log10(F_rel)
    return mV_mean + A_V + (dm - float(np.mean(dm)))
