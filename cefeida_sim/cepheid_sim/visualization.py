
from __future__ import annotations
import numpy as np
import plotly.graph_objects as go

def kelvin_to_rgb(T: float) -> tuple[int, int, int]:
    """Aproximación simple (Tanner Helland). T en Kelvin."""
    T = np.clip(T, 1000.0, 40000.0) / 100.0
    # Red
    if T <= 66:
        red = 255
    else:
        red = T - 60.0
        red = 329.698727446 * (red ** -0.1332047592)
        red = np.clip(red, 0, 255)
    # Green
    if T <= 66:
        green = 99.4708025861 * np.log(T) - 161.1195681661
    else:
        green = T - 60.0
        green = 288.1221695283 * (green ** -0.0755148492)
    green = np.clip(green, 0, 255)
    # Blue
    if T >= 66:
        blue = 255
    elif T <= 19:
        blue = 0
    else:
        blue = 138.5177312231 * np.log(T - 10.0) - 305.0447927307
        blue = np.clip(blue, 0, 255)
    return int(red), int(green), int(blue)

def sphere_figure(radius_rel: float, T_color: float):
    """Crea una esfera 3D con radio relativo y color desde Teff."""
    # Mallado de esfera
    u = np.linspace(0, 2*np.pi, 60)
    v = np.linspace(0, np.pi, 30)
    x = radius_rel * np.outer(np.cos(u), np.sin(v))
    y = radius_rel * np.outer(np.sin(u), np.sin(v))
    z = radius_rel * np.outer(np.ones_like(u), np.cos(v))

    r,g,b = kelvin_to_rgb(T_color)
    color = f"rgb({r},{g},{b})"

    surface = go.Surface(
        x=x, y=y, z=z,
        showscale=False,
        colorscale=[[0, color],[1, color]],
        cmin=0, cmax=1,
        opacity=1.0,
    )
    fig = go.Figure(data=[surface])
    fig.update_scenes(
        xaxis_visible=False, yaxis_visible=False, zaxis_visible=False,
        aspectmode="data"
    )
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))
    return fig
