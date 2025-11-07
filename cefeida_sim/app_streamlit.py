import numpy as np
import streamlit as st
from cepheid_sim.kinematics import VRadOptions
from cepheid_sim.pipeline import compute_coupled_curves
import plotly.graph_objects as go
from plotly.subplots import make_subplots

st.set_page_config(page_title="Simulador 3D de Cefeidas", layout="wide")

st.title("Simulador 3D de Estrellas Cefeidas — MVP")

def make_combined_anim(phase_grid, mV, vr, phi_seq, mV_seq, vr_seq, R_seq, T_seq,
                       scale_R=6.0, boost_K=1500, frame_ms=33):

    def kelvin_to_rgb(T):
        T = np.clip(T, 1000.0, 40000.0)/100.0
        if T <= 66: red = 255
        else:
            red = 329.698727446 * ((T-60.0)**-0.1332047592)
            red = np.clip(red, 0, 255)
        if T <= 66:
            green = 99.4708025861*np.log(T) - 161.1195681661
        else:
            green = 288.1221695283 * ((T-60.0)**-0.0755148492)
        green = np.clip(green, 0, 255)
        if T >= 66: blue = 255
        elif T <= 19: blue = 0
        else:
            blue = 138.5177312231*np.log(T-10.0) - 305.0447927307
            blue = np.clip(blue, 0, 255)
        return int(red), int(green), int(blue)

    # Esfera base
    u = np.linspace(0, 2*np.pi, 60)
    v = np.linspace(0, np.pi, 30)
    Xu = np.outer(np.cos(u), np.sin(v))
    Yu = np.outer(np.sin(u), np.sin(v))
    Zu = np.outer(np.ones_like(u), np.cos(v))

    R_mean = float(np.mean(R_seq))
    R_rel_seq = R_seq / R_mean
    R_vis_seq = 1.0 + (R_rel_seq - 1.0) * scale_R
    rmax = float(np.max(R_vis_seq))
    lim = rmax * 1.1
    T_mean = float(np.mean(T_seq))

    R_rel0 = float(R_seq[0] / R_mean)
    R_vis0 = 1.0 + (R_rel0 - 1.0) * scale_R
    dT0 = float(T_seq[0] - T_mean)
    T_vis0 = float(np.clip(T_seq[0] + np.sign(dT0)*min(abs(dT0)+boost_K, 4000.0), 1000.0, 40000.0))
    r0, g0, b0 = kelvin_to_rgb(T_vis0)
    color0 = f"rgb({r0},{g0},{b0})"

    y1_min, y1_max = float(np.min(mV)), float(np.max(mV))
    y2_min, y2_max = float(np.min(vr)), float(np.max(vr))

    fig = make_subplots(
        rows=1, cols=3,
        specs=[[{"type": "xy"}, {"type": "xy"}, {"type": "scene"}]],
        column_widths=[0.33, 0.33, 0.34],
        horizontal_spacing=0.06,
        subplot_titles=("Curva de Luz (V)", "Velocidad radial", "Esfera 3D")
    )

    fig.update_yaxes(autorange="reversed", row=1, col=1)

    # 0: mV (línea)
    fig.add_trace(go.Scatter(x=phase_grid, y=mV, mode="lines", name="mV"), row=1, col=1)
    # 1: vline mV
    fig.add_trace(go.Scatter(x=[phi_seq[0], phi_seq[0]], y=[y1_min, y1_max], mode="lines",
                             line=dict(dash="dash"), showlegend=False, name="ϕ"), row=1, col=1)
    # 2: marker mV
    fig.add_trace(go.Scatter(x=[phi_seq[0]], y=[float(mV_seq[0])], mode="markers",
                             marker=dict(size=10), name="fase mV"), row=1, col=1)

    # 3: vr (línea)
    fig.add_trace(go.Scatter(x=phase_grid, y=vr, mode="lines", name="vrad [km/s]"), row=1, col=2)
    # 4: vline vr
    fig.add_trace(go.Scatter(x=[phi_seq[0], phi_seq[0]], y=[y2_min, y2_max], mode="lines",
                             line=dict(dash="dash"), showlegend=False, name="ϕ"), row=1, col=2)
    # 5: marker vr
    fig.add_trace(go.Scatter(x=[phi_seq[0]], y=[float(vr_seq[0])], mode="markers",
                             marker=dict(size=10), name="fase vr"), row=1, col=2)

    # 6: esfera
    fig.add_trace(go.Surface(
        x=R_vis0*Xu, y=R_vis0*Yu, z=R_vis0*Zu,
        showscale=False, cmin=0, cmax=1,
        colorscale=[[0, color0], [1, color0]], opacity=1.0, name="esfera"
    ), row=1, col=3)

    frames = []
    for i, (phi, mv, vv, Ri, Ti) in enumerate(zip(phi_seq, mV_seq, vr_seq, R_seq, T_seq)):
        R_rel = float(Ri / R_mean)
        R_vis = 1.0 + (R_rel - 1.0) * scale_R
        dT = float(Ti - T_mean)
        T_vis = float(np.clip(Ti + np.sign(dT)*min(abs(dT)+boost_K, 4000.0), 1000.0, 40000.0))
        r, g, b = kelvin_to_rgb(T_vis)
        col = f"rgb({r},{g},{b})"

        frames.append(go.Frame(
            name=str(i),
            data=[
                go.Scatter(x=phase_grid, y=mV),                # 0
                go.Scatter(x=[phi, phi], y=[y1_min, y1_max]),  # 1
                go.Scatter(x=[phi], y=[float(mv)]),            # 2
                go.Scatter(x=phase_grid, y=vr),                # 3
                go.Scatter(x=[phi, phi], y=[y2_min, y2_max]),  # 4
                go.Scatter(x=[phi], y=[float(vv)]),            # 5
                go.Surface(                                     # 6
                    x=R_vis*Xu, y=R_vis*Yu, z=R_vis*Zu,
                    showscale=False, cmin=0, cmax=1,
                    colorscale=[[0, col], [1, col]], opacity=1.0
                )
            ]
        ))

    fig.frames = frames

    fig.update_layout(
        margin=dict(l=10, r=10, t=40, b=10),
        xaxis_title="fase", xaxis2_title="fase",
        yaxis_title="mV", yaxis2_title="vrad [km/s]",
        scene=dict(
            xaxis=dict(range=[-lim, lim], autorange=False, visible=False),
            yaxis=dict(range=[-lim, lim], autorange=False, visible=False),
            zaxis=dict(range=[-lim, lim], autorange=False, visible=False),
            aspectmode="cube",
            camera=dict(eye=dict(x=1.6, y=1.6, z=1.2))
        ),
        updatemenus=[dict(
            type="buttons", showactive=False,
            buttons=[
                dict(label="▶ Play", method="animate",
                     args=[None, {"frame": {"duration": frame_ms, "redraw": True},
                                  "fromcurrent": True, "transition": {"duration": 0}}]),
                dict(label="⏸ Pause", method="animate",
                     args=[[None], {"mode": "immediate",
                                    "frame": {"duration": frame_ms, "redraw": True},
                                    "transition": {"duration": 0}}])
            ],
            x=0.01, y=1.12, xanchor="left", yanchor="top"
        )]
    )

    return fig


# -------------------- Sidebar --------------------
with st.sidebar:
    st.header("Parámetros")
    P = st.number_input("Periodo P [días]", value=10.15, min_value=0.5, step=0.05)
    R0 = st.number_input("Radio medio R₀ [R☉]", value=60.0, min_value=1.0, step=1.0)
    T0 = st.number_input("Temperatura media T₀ [K]", value=5800.0, min_value=3000.0, step=50.0)
    v_gamma = st.number_input("Velocidad sistémica v_γ [km/s]", value=0.0, step=0.1)
    p_factor = st.slider("Factor de proyección p", min_value=1.20, max_value=1.60, value=1.36, step=0.01)
    mV_mean = st.number_input("Magnitud V media", value=3.9, step=0.1)
    A_V = st.number_input("Extinción A_V [mag]", value=0.2, step=0.01)

    perfil_v = st.selectbox(
        "Perfil de v_rad",
        ["Sinusoidal", "Asimétrica (skew)", "Dos tramos (triangular)"],
        index=0
    )
    v_amp = st.slider("Amplitud v_rad [km/s]", 0.1, 60.0, 25.0, 0.1)

    if "Asimétrica" in perfil_v:
        skew = st.slider("Asimetría (skew)", 0.0, 3.0, 0.35, 0.01)
    else:
        skew = 0.35

    if "tramos" in perfil_v:
        frac_rise = st.slider("Fracción de fase en ascenso", 0.1, 0.9, 0.35, 0.01)
    else:
        frac_rise = 0.35

    # **Controles de acople fotométrico**
    beta_V = st.slider("β fotométrico (banda V)", 1.0, 2.5, 1.6, 0.1)
    phiT_lag = st.slider("Desfase T↔R [ciclos]", 0.0, 0.25, 0.08, 0.005)
    
    
    # --- Acople térmico ---
    alpha_T = st.slider("Exponente térmico α (T ∝ R^{-α})", 0.2, 1.0, 0.55, 0.01)
    dT_pp   = st.number_input("Amplitud T pico a pico [K]", value=1200.0, min_value=0.0, step=50.0)
    
    # --- Modo fotométrico ---
    phot_mode_ui = st.radio(
        "Modo fotometría",
        ["R²·T^β", "Planck + filtro V"],
        index=0
    )
    photometry_mode = "bandpass" if phot_mode_ui.startswith("Planck") else "powerlaw"
    
    # Parámetros del modo bandpass
    if photometry_mode == "bandpass":
        band_center_nm = st.number_input("Centro banda [nm]", 400.0, 800.0, 550.0, 1.0)
        band_width_nm  = st.number_input("Ancho banda [nm]", 10.0, 200.0, 88.0, 1.0)
        band_nlam      = st.slider("Puntos espectrales", 100, 1000, 400, 50)
    else:
        band_center_nm = 550.0
        band_width_nm  = 88.0
        band_nlam      = 400

    n_frames = st.slider("Cuadros por ciclo (animación)", 30, 180, 100, 5)
    scale_R = st.slider("Escala visual del radio", 1.0, 10.0, 6.0, 0.5)
    boost_K = st.slider("Realce de color (K)", 0, 5000, 1500, 50)
    fps = st.select_slider("FPS animación", options=[15, 24, 30, 60], value=30)

# -------------------- Cálculo acoplado --------------------
phase_grid = np.linspace(0, 1, 600, endpoint=False)

mode = "sin" if "Sinusoidal" in perfil_v else ("skew" if "Asimétrica" in perfil_v else "tri")
vopts = VRadOptions(mode=mode, v_amp=v_amp, skew=skew, frac_rise=frac_rise)

out = compute_coupled_curves(
    P_days=P,
    R0_Rsun=R0,
    T0_K=T0,
    v_gamma_kms=v_gamma,
    p_factor=p_factor,
    mV_mean=mV_mean,
    A_V=A_V,  # pásalo aquí; ya no lo sumes luego a mano

    # térmico
    alpha_T=alpha_T,
    T_lag_cycles=phiT_lag,
    target_dT_pp=dT_pp,

    # fotometría
    beta_V=beta_V,  # usado solo en "powerlaw"

    # malla y cinemática
    phase_grid=phase_grid,
    n_frames=n_frames,
    vopts=vopts,
)

# Curvas en malla y secuencias (ya acopladas entre sí)
vr       = out["vr"]
mV_curve = out["mV"]
phi_seq  = out["phi_seq"]
vr_seq   = out["vr_seq"]
R_seq    = out["R_seq"]
T_seq    = out["T_seq"]
mV_seq   = out["mV_seq"]

# -------------------- Render --------------------
fig_all = make_combined_anim(
    phase_grid, mV_curve, vr,
    phi_seq, mV_seq, vr_seq, R_seq, T_seq,
    scale_R=scale_R, boost_K=boost_K,
    frame_ms=int(1000 / fps)
)

st.plotly_chart(fig_all, use_container_width=True, config={"displayModeBar": False})

st.markdown("---")
st.markdown(
    f"**Resumen:** ⟨mV⟩={mV_curve.mean():.3f}, Δm={mV_curve.max()-mV_curve.min():.3f} | "
    f"⟨vr⟩={vr.mean():.2f} km/s, Δv={vr.max()-vr.min():.2f} km/s | "
    f"⟨R⟩≈{np.mean(R_seq):.1f} R☉, ΔR≈{np.max(R_seq)-np.min(R_seq):.1f} R☉, "
    f"⟨T⟩≈{np.mean(T_seq):.0f} K"
)