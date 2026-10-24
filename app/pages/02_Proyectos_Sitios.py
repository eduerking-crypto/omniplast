import sys
from pathlib import Path

_here = Path(__file__).resolve()
for _p in _here.parents:
    if (_p / "app" / "ui" / "bootstrap.py").exists():
        for _s in (str(_p), str(_p / "app")):
            if _s not in sys.path:
                sys.path.insert(0, _s)
        break

import pandas as pd
import streamlit as st
from ui.bootstrap import (chart, exportar_dfs, footer, page_header, resumen_box,
                          section_title, tabla)

from omniplast.core import indices
from omniplast.db import dao
from omniplast.reporting import figures

page_header("Proyectos, estaciones y campañas",
            "Registro maestro georreferenciado · estructura lista para sus datos reales")

sites = dao.sites_df()
samples = dao.samples_df()
camps = dao.campaigns_df()
mp = dao.mp_df()

t1, t2, t3 = st.tabs(["Estaciones", "Campañas", "Cargar datos propios"])

with t1:
    section_title("Estaciones del proyecto activo")
    tabla(sites)
    section_title("Carga de MP por estación y matriz (todas las campañas)")
    if not mp.empty:
        piv = mp.pivot_table(index="sitio", columns="matriz", values="abundancia",
                             aggfunc="median").round(0)
        tabla(piv.reset_index())
        chart(figures.heatmap(piv.T, title="Mediana de MP: sitio × matriz"))

with t2:
    section_title("Campañas")
    tabla(camps)
    if not mp.empty:
        serie = mp.groupby(["campana", "matriz"], observed=True)["abundancia"].median().reset_index()
        chart(figures.bars(serie, x="campana", y="abundancia", color="matriz",
                           title="Mediana de MP por campaña y matriz", ytitle="abundancia"))

with t3:
    section_title("Plantillas oficiales de la base de datos")
    resumen_box("Descargue la <b>PLANTILLA</b> (estructura vacía con validaciones "
                "desplegables, catálogos y guía LEEME) o el <b>EJEMPLO LLENO</b> "
                "(dataset demo de la Ciénaga de La Virgen, DATOS SINTÉTICOS) para ver "
                "cómo se ve un proyecto completo y reproducir cualquier análisis de inmediato.")
    from omniplast.reporting import plantillas as _plt

    @st.cache_data(ttl=600, show_spinner="Generando plantilla…")
    def _bytes_plantilla():
        return _plt.generar_plantilla_bytes()

    @st.cache_data(ttl=600, show_spinner="Generando ejemplo…")
    def _bytes_ejemplo():
        return _plt.generar_ejemplo_bytes()

    cA, cB = st.columns(2, gap="large")
    with cA:
        st.markdown('<div class="omp-card tall"><h4>① PLANTILLA_OMNIPLAST.xlsx</h4>'
                    '<p>Estructura vacía: 9 hojas (Sitios → Microbioma_16S), listas '
                    'desplegables de polímeros/morfologías/colores/genes, catálogo de '
                    'parámetros con unidades y guía de análisis en LEEME.</p></div>',
                    unsafe_allow_html=True)
        st.download_button("Descargar plantilla vacía", _bytes_plantilla(),
                           file_name="PLANTILLA_OMNIPLAST.xlsx",
                           mime="application/vnd.openxmlformats-officedocument."
                                "spreadsheetml.sheet")
    with cB:
        st.markdown('<div class="omp-card tall"><h4>② EJEMPLO_OMNIPLAST.xlsx</h4>'
                    '<p>El mismo formato lleno con 10 estaciones, 4 campañas, 144 muestras, '
                    '1.248 registros de MP, parámetros ambientales, qPCR del panel de '
                    'genes y microbioma 16S: ábralo en la página 03 o ejecute cualquier '
                    'agente para ver el análisis completo en acción.</p></div>',
                    unsafe_allow_html=True)
        st.download_button("Descargar ejemplo lleno (demo)", _bytes_ejemplo(),
                           file_name="EJEMPLO_OMNIPLAST.xlsx",
                           mime="application/vnd.openxmlformats-officedocument."
                                "spreadsheetml.sheet")

    section_title("Carga de sus datos propios")
    resumen_box("Cargue sus CSV/Excel con la estructura de las plantillas (sitios, "
                "muestras, microplasticos, parametros_ambientales, biota, qpcr, "
                "microbioma_16S). El sistema validará las columnas mínimas.",
                "warn")
    up = st.file_uploader("Archivo de estaciones (sitios.csv)", type=["csv", "xlsx"])
    if up is not None:
        try:
            df = pd.read_csv(up) if up.name.endswith(".csv") else pd.read_excel(up)
            requeridas = {"codigo", "nombre", "lat", "lon"}
            faltan = requeridas - set(df.columns.str.lower())
            if faltan:
                resumen_box(f"Faltan columnas requeridas: {faltan}", "error")
            else:
                st.session_state["upload_sitios"] = df
                tabla(df.head(30))
                st.success(f"{len(df)} estaciones listas para validar. Use la API "
                           "(/api/sitios) o contacte al administrador del proyecto para "
                           "el volcado formal a la base.")
        except Exception as e:
            resumen_box(f"Error leyendo el archivo: {e}", "error")
    st.markdown("**Plantilla de columnas mínimas:**")
    st.code("codigo, nombre, lat, lon, zona, descripcion", language="text")
    exportar_dfs("plantillas_estaciones", {"sitios_ejemplo": sites})

footer()
