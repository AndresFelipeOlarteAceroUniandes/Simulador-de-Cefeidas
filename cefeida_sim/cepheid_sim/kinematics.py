# cepheid_sim/kinematics.py
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

R_SUN_KM = 695_700.0  # km

@dataclass(frozen=True)
class VRadOptions:
    mode: str = "sin"
    v_amp: float = 25.0
    skew: float = 0.35
    frac_rise: float = 0.35

def _phase_skew(phi: np.ndarray, skew: float) -> np.ndarray:
    phi_ = (1.0 - skew) * phi + skew * (phi ** 2)
    return np.mod(phi_, 1.0)

def _triangular(phi: np.ndarray, frac_rise: float) -> np.ndarray:
    fr = np.clip(float(frac_rise), 1e-3, 1 - 1e-3)
    out = np.empty_like(phi, dtype=float)
    # subida
    mask_up = (phi < fr)
    out[mask_up] = (phi[mask_up] / fr)  # 0..1
    # bajada
    mask_dn = ~mask_up
    out[mask_dn] = 1.0 - (phi[mask_dn] - fr) / (1.0 - fr)  # 1..0
    out = 2.0 * out - 1.0  # a [-1, 1]
    return out

def vrad_profile(phi: np.ndarray, v_gamma: float, opts: VRadOptions) -> np.ndarray:
    phi = np.mod(phi, 1.0).astype(float)

    if opts.mode == "sin":
        base = np.sin(2.0 * np.pi * phi)  # [-1, 1]
    elif opts.mode == "skew":
        phis = _phase_skew(phi, float(opts.skew))
        base = np.sin(2.0 * np.pi * phis)  # [-1, 1] pero con asimetría de fase
    elif opts.mode == "tri":
        base = _triangular(phi, float(opts.frac_rise))  # [-1, 1]
    else:
        raise ValueError(f"VRadOptions.mode inválido: {opts.mode}")

    return v_gamma + opts.v_amp * base  # [km/s]

def integrate_radius(phi: np.ndarray,
                     v_puls_kms: np.ndarray,
                     P_days: float,
                     R0_Rsun: float) -> np.ndarray:
    """
    Integra la velocidad pulsacional (ya con p-factor aplicado) para obtener R(ϕ) en R_sun.
    - phi: fase 0..1 (malla monótona, típicamente equiespaciada)
    - v_puls_kms: velocidad pulsacional en km/s (centra en 0 aprox.)
    - P_days: período en días
    - R0_Rsun: radio medio para recentrar el ciclo
    """
    phi = np.asarray(phi, dtype=float)
    v  = np.asarray(v_puls_kms, dtype=float)

    # dt por paso de fase (asumiendo malla uniforme):
    dphi = np.diff(np.r_[phi, phi[0] + 1.0])  # cierra el ciclo
    dt   = dphi * (P_days * 86400.0)          # [s]

    # dR [km] = v [km/s] * dt [s]
    dR_km = v * dt
    R_km  = np.cumsum(dR_km)
    # cerrar ciclo: eliminar deriva lineal (resta la recta que une extremo con inicio)
    drift = np.linspace(0.0, R_km[-1], num=R_km.size)
    R_km_corr = R_km - drift

    # escala a R_sun y centra en R0
    R_Rsun = R0_Rsun + (R_km_corr / R_SUN_KM)
    return R_Rsun

