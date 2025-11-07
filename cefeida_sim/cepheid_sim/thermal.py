from __future__ import annotations
import numpy as np
from typing import Optional, Tuple

def _phase_shift_fractional(y: np.ndarray, shift_cycles: float) -> np.ndarray:
    """
    Desplaza periódicamente 'y' una fracción de ciclo (fase en [0,1)).
    Usa interpolación lineal y extensión de un punto para evitar bordes.
    """
    y = np.asarray(y, dtype=float)
    n = y.size
    if n == 0:
        return y
    x = np.arange(n, dtype=float)
    x_new = (x - shift_cycles * n) % n
    x_ext = np.concatenate([x[-1:] - n, x, x[:1] + n])
    y_ext = np.concatenate([y[-1:],    y, y[:1]])
    return np.interp(x_new, x_ext, y_ext)

def temperature_from_radius(
    R_rel: np.ndarray,          # R(ϕ) / R0 (adimensional)
    T0: float,                  # temperatura media (K)
    alpha: float = 0.5,         # T/T0 ≈ (R_rel)^(-alpha)
    lag_cycles: float = 0.08,   # desfase de T respecto a R (en ciclos)
    dT_pp: Optional[float] = None,                 # amplitud pico-a-pico deseada (K)
    T_bounds: Optional[Tuple[float, float]] = None # (T_min, T_max) opcional
) -> np.ndarray:
    """
    Calcula T(ϕ) acoplada a R(ϕ) con desfase y control explícito de amplitud.

    Pasos:
      1) T_rel_raw = (R_rel)^(-alpha)
      2) Aplica desfase en fase (lag_cycles)
      3) Convierte a K y recentra para que <T> = T0
      4) Si T_bounds: reescala para que min/max sean T_min/T_max
         Si dT_pp: reescala para que el pico-a-pico sea dT_pp

    Prioridad: T_bounds (si se da) tiene prioridad sobre dT_pp.
    """
    R_rel = np.asarray(R_rel, dtype=float)
    # 1) ley térmica acoplada (evita extremos numéricos)
    T_rel = np.power(np.clip(R_rel, 1e-6, 1e6), -alpha)

    # 2) desfase de temperatura respecto a radio
    T_rel = _phase_shift_fractional(T_rel, lag_cycles)

    # 3) pasa a Kelvin y recentra (asegura <T> = T0)
    T_raw = T0 * T_rel
    T = T0 + (T_raw - float(np.mean(T_raw)))  # ahora el promedio es T0

    # 4) normalización de amplitud
    if T_bounds is not None:
        T_min, T_max = T_bounds
        if T_max <= T_min:
            # caso degenerado: cae a constante T0
            return np.full_like(T, T0, dtype=float)

        ptp = float(np.max(T) - np.min(T))
        if ptp > 0:
            # reescala desde (T - <T>) para mantener el centro en T0
            scale = (T_max - T_min) / ptp
            T = T0 + (T - float(np.mean(T))) * scale
            # ajusta el centro exactamente al punto medio (por si numéricamente se movió)
            T = T - (float(np.mean(T)) - 0.5 * (T_min + T_max))
        else:
            T = np.full_like(T, 0.5 * (T_min + T_max), dtype=float)
        return T

    if dT_pp is not None and dT_pp > 0:
        ptp = float(np.max(T) - np.min(T))
        if ptp > 0:
            scale = dT_pp / ptp
            T = T0 + (T - float(np.mean(T))) * scale
        else:
            # si no hay variación, crea una pequeña sinusoidal alrededor de T0
            phi = np.linspace(0.0, 1.0, len(T), endpoint=False)
            T = T0 + 0.5 * dT_pp * np.sin(2 * np.pi * phi)
        return T

    # si no se pidió control de amplitud, devuelve la señal re-centrada
    return T
