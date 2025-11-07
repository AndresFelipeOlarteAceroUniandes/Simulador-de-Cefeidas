# 🌟 Simulador 3D de Estrellas Cefeidas — MVP

Este repositorio contiene una implementación interactiva en **Python + Streamlit** de un simulador físico-minimalista de **estrellas Cefeidas clásicas**, que modela el acoplamiento entre velocidad radial, variaciones del radio, temperatura efectiva y respuesta fotométrica en banda V.  
El proyecto permite visualizar curvas de luz, velocidad radial y una representación 3D dinámica de la estrella.

---

## 🚀 Características principales

- Simulación acoplada:  
  \(v_{\text{rad}}(\phi) \rightarrow v_p(\phi) \rightarrow R(\phi) \rightarrow T(\phi) \rightarrow m_V(\phi)\)
- Módulos modulares: `kinematics`, `thermal`, `photometry`, `pipeline`
- Curvas de luz en banda V (modelo de cuerpo negro o ley de potencia)
- Interfaz interactiva en **Streamlit**
- Soporte para exportación de resultados
- Ideal para **docencia** y **visualización de física estelar**

---

## 🧩 Estructura del repositorio

cepheid_sim/
│
├── app_streamlit.py # Interfaz principal (Streamlit)
├── cepheid_sim/
│    ├── init.py
│    ├── kinematics.py # Perfil de velocidad radial e integración de radio
│    ├── thermal.py # Cálculo de temperatura efectiva
│    ├── photometry.py # Fotometría y magnitudes
│    ├── pipeline.py # Flujo acoplado de simulación
│    ├─ physics.py                 
│    ├─ visualization.py  
├── requirements.txt # Dependencias del entorno
└── README.md

## Ubicar directorio 
```bash``````````````````````````````````````

# Windows
cd c:\ubicacion_del_directorio

# macOS/Linux
cd c:/ubicacion_del_directorio
`````````````````````````````````````````````



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
`````````````````````````````````````````````

