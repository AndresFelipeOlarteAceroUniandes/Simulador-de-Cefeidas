# Simulador 3D de Estrellas Cefeidas (MVP)

Este repositorio contiene un simulador funcional de estrellas cefeidas. 
Incluye:
- Módulo físico (radio, temperatura, luminosidad, velocidad radial).
- Fotometría simplificada en banda V (magnitud relativa con extinción).
- Visualización 3D con Plotly (esfera que cambia de tamaño y color).
- App interactiva con **Streamlit**.

## Instalación

```bash``````````````````````````````````````
python -m venv .venv

# Windows
.venv\Scripts\activate

# macOS/Linux
source .venv/bin/activate

pip install -r requirements.txt
`````````````````````````````````````````````

## Ejecutar la app (Streamlit)

```bash``````````````````````````````````````
streamlit run app_streamlit.py
```

## Ejecutar el demo por consola

```bash
python -m examples.run_demo
`````````````````````````````````````````````


## Estructura

```
cepheid_sim/
  cepheid_sim/
    __init__.py
    physics.py
    photometry.py
    visualization.py
  examples/
    run_demo.py
  app_streamlit.py
  requirements.txt
  README.md
```

## Qué hace cada módulo

- `physics.py`: define `CepheidParams` y funciones para `radius`, `teff`, `luminosity`, `vrad`.
- `photometry.py`: Planck simplificado y magnitud instrumental (normalizada al promedio y con A_V).
- `visualization.py`: figura 3D (Plotly) con esfera cuyo radio y color dependen de la fase.

## Pendientes y mejoras sugeridas

- Integrar **relación Período–Luminosidad** para fijar `M_V` y el cero fotométrico físico.
- Añadir **curvas de transmisión reales** (banda V, B, R, I) y flujo integrado más realista.
- Validar **fase de luz vs. velocidad radial** y coherencia con observaciones.
- Permitir múltiples armónicos ajustables (AR_k, AT_k, φ_k) desde la UI.
- Exportar curvas simuladas a CSV.
