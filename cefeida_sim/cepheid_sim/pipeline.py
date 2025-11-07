from __future__ import annotations
import numpy as np
from typing import Dict

from .kinematics import VRadOptions, vrad_profile, integrate_radius
from .thermal import temperature_from_radius
from .photometry import (
    flux_V_rel,
    flux_band_rel_from_RT,
    magnitude_from_flux_rel,
)

def _rescale_T_ptp(T: np.ndarray, T0: float, target_dT_pp: float | None) -> np.ndarray:
    """Reescala T(ϕ) para que tenga un pico-a-pico deseado (si se pide)."""
    if target_dT_pp is None:
        return T
    T = np.asarray(T, dtype=float)
    # Centra alrededor de T0 y reescala
    T_centered = T - float(np.mean(T))
    ptp = float(np.max(T_centered) - np.min(T_centered))
    if ptp <= 0:
        return T  # nada que escalar
    scale = target_dT_pp / ptp
    return T0 + T_centered * scale

def compute_coupled_curves(
    *,
    # Parámetros físicos/base
    P_days: float,
    R0_Rsun: float,
    T0_K: float,
    v_gamma_kms: float,
    p_factor: float,
    mV_mean: float,
    A_V: float = 0.0,

    # Temperatura acoplada
    alpha_T: float = 0.5,       # T ∝ R^{-alpha_T}
    T_lag_cycles: float = 0.08, # desfase de T respecto a R (en ciclos)
    target_dT_pp: float | None = None,  # amplitud T pico-a-pico (K), opcional

    # Fotometría
    photometry_mode: str = "powerlaw",  # "powerlaw" o "bandpass"
    beta_V: float = 1.6,                # solo para "powerlaw"
    band_center_nm: float = 550.0,      # solo para "bandpass"
    band_width_nm: float = 88.0,
    band_nlam: int = 400,

    # Mallado y kinemática
    phase_grid: np.ndarray | None = None,
    n_frames: int = 100,
    vopts: VRadOptions | None = None,
) -> Dict[str, np.ndarray]:
    """
    Devuelve:
      - phi, vr, vp, R, T, mV   (en la malla 'phi')
      - phi_seq, vr_seq, vp_seq, R_seq, T_seq, mV_seq (en frames)
    """
    if phase_grid is None:
        phi = np.linspace(0.0, 1.0, 600, endpoint=False)
    else:
        phi = np.mod(np.asarray(phase_grid, dtype=float), 1.0)

    phi_seq = np.linspace(0.0, 1.0, int(n_frames), endpoint=False)

    # --- 1) v_rad observada ---
    if vopts is None:
        raise ValueError("Se requiere vopts (VRadOptions) para definir el perfil de v_rad.")
    vr     = vrad_profile(phi,     v_gamma_kms, vopts)
    vr_seq = vrad_profile(phi_seq, v_gamma_kms, vopts)

    # --- 2) velocidad pulsacional (aplica p-factor) ---
    vp     = p_factor * (vr     - v_gamma_kms)
    vp_seq = p_factor * (vr_seq - v_gamma_kms)

    # --- 3) integra a R(ϕ) (coherente con p-factor) ---
    R     = integrate_radius(phi,     vp,     P_days, R0_Rsun)
    R_seq = integrate_radius(phi_seq, vp_seq, P_days, R0_Rsun)

    # Normaliza para fotometría relativa
    R_rel     = R     / float(np.mean(R))
    R_rel_seq = R_seq / float(np.mean(R_seq))

    # --- 4) T(ϕ) acoplada a R con desfase y exponente ---
    T     = temperature_from_radius(R_rel,     T0_K, alpha=alpha_T, lag_cycles=T_lag_cycles)
    T_seq = temperature_from_radius(R_rel_seq, T0_K, alpha=alpha_T, lag_cycles=T_lag_cycles)

    # Reescala amplitud de T si se pide (para igualar, p.ej., 5500–6800 K)
    T     = _rescale_T_ptp(T,     T0_K, target_dT_pp)
    T_seq = _rescale_T_ptp(T_seq, T0_K, target_dT_pp)

    # --- 5) Fotometría: flujo relativo y magnitud ---
    if photometry_mode.lower() == "bandpass":
        F_rel     = flux_band_rel_from_RT(R_rel,     T,     band_center_nm, band_width_nm, band_nlam)
        F_rel_seq = flux_band_rel_from_RT(R_rel_seq, T_seq, band_center_nm, band_width_nm, band_nlam)
    else:
        # Ruta rápida por potencia
        F_rel     = flux_V_rel(R_rel,     T / T0_K, beta=beta_V)
        F_rel_seq = flux_V_rel(R_rel_seq, T_seq / T0_K, beta=beta_V)

    mV     = magnitude_from_flux_rel(F_rel,     mV_mean, A_V=A_V)
    mV_seq = magnitude_from_flux_rel(F_rel_seq, mV_mean, A_V=A_V)

    return dict(
        # malla
        phi=phi, vr=vr, vp=vp, R=R, T=T, mV=mV,
        # frames
        phi_seq=phi_seq, vr_seq=vr_seq, vp_seq=vp_seq, R_seq=R_seq, T_seq=T_seq, mV_seq=mV_seq,
        # útiles para depuración/diagnóstico
        R_rel=R_rel, R_rel_seq=R_rel_seq, F_rel=F_rel, F_rel_seq=F_rel_seq
    )
