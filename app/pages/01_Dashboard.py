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
from ui.bootstrap import (chart, footer, kpi_row, page_header, resumen_box,
                          section_title, tabla)
from ui.components import mapa_estaciones

from omniplast.core import indices
from omniplast.db import dao
from omniplast.reporting import figures
from omniplast.utils import fmt_num

page_header("Dashboard del ecosistema",
           "Panorama integrado de carga de microplásticos, calidad del agua e índices compuestos · DATOS DEMO")

f1, f2, f3, _ = st.columns([2, 2, 2, 5])
mp = dao.mp_df()
if mp.empty:
    resumen_box("La base de datos está vacía. Ejecute `python scripts/generate_demo.py` "
                "o cargue datos desde la página Proyectos.", "warn")
    st.stop()
matrices = sorted(mp["matriz"].unique())
sitios = sorted(mp["sitio"].unique())
temporadas = sorted(mp["temporada"].dropna().unique())
matriz_f = f1.selectbox("Matriz", matrices, index=0)
sitio_f = f2.selectbox("Estación", ["Todas"] + sitios)
temp_f = f3.selectbox("Temporada", ["Todas"] + temporadas)

d = mp[(mp["matriz"] == matriz_f)]
if sitio_f != "Todas":
    d = d[d["sitio"] == sitio_f]
if temp_f != "Todas":
    d = d[d["temporada"] == temp_f]

per_sample = d.groupby(["sample_id", "sitio", "zona", "campana", "temporada"],
                       observed=True)["abundancia"].sum().reset_index()

mediana = per_sample["abundancia"].median()
tier = indices.mp_tier(mediana, matriz_f)
color_tier = {"Baja": "#3E9B4F", "Moderada": "#00A896", "Alta": "#E9B44C",
              "Muy alta": "#FF6B6B", "Extremadamente alta": "#FF6B6B"}.get(tier, "#00A896")

agua_env = dao.env_wide()
agua_env = agua_env[agua_env["matriz"] == "agua"]
wqi_vals = []
for sitio, g in agua_env.groupby("sitio"):
    w = indices.wqi_nsf(g.select_dtypes("number").mean(numeric_only=True).to_dict())
    if w["wqi"] == w["wqi"]:
        wqi_vals.append(w["wqi"])

kpi_row([
    (fmt_num(mediana, 0) + " p/m³" if matriz_f == "agua" else fmt_num(mediana, 0) + " p/kg",
     "Mediana de MP", f"Clase: {tier}", color_tier),
    (fmt_num(per_sample["abundancia"].max(), 0), "Máximo registrado",
     per_sample.loc[per_sample["abundancia"].idxmax(), "sitio"], "#FF6B6B"),
    (fmt_num(pd.Series(wqi_vals).mean(), 1) if wqi_vals else "—", "WQI-NSF medio",
     "calidad fisicoquímica", "#0A2540"),
    (fmt_num(len(per_sample), 0), "Muestras activas", f"filtro: {matriz_f}", "#05BFDB"),
])

section_title("Mapa de estaciones")
mapa_fig = mapa_estaciones(per_sample, matriz_f)
if mapa_fig is not None:
    from streamlit_folium import st_folium
    st_folium(mapa_fig, height=420, returned_objects=[])

section_title("Distribución espacial y temporal")
c1, c2 = st.columns(2, gap="large")
with c1:
    por_sitio = per_sample.groupby("sitio", observed=True)["abundancia"].median().sort_values()
    chart(figures.bars(por_sitio.reset_index(), "sitio", "abundancia",
                       title=f"Mediana de MP por estación ({matriz_f})",
                       ytitle="partículas"))
with c2:
    comp = d.groupby("polimero", observed=True)["abundancia"].sum().sort_values(ascending=False).head(8)
    chart(figures.pie(list(comp.index), list(comp.values), title="Composición de polímeros (top 8)"))

c3, c4 = st.columns(2, gap="large")
with c3:
    chart(figures.box(per_sample, x="zona", y="abundancia", title="MP por zona de impacto",
                      ytitle="abundancia"))
with c4:
    serie = per_sample.groupby("campana", observed=True)["abundancia"].median().reset_index()
    serie = serie.sort_values("campana")
    chart(figures.line(serie, x="campana", y="abundancia", title="Evolución por campaña",
                       ytitle="mediana"))

section_title("Radiografía por estación")
rows = []
for sitio, g in per_sample.groupby("sitio", observed=True):
    zona = g["zona"].iloc[0]
    agua_env_site = agua_env[agua_env["sitio"] == sitio]
    wqi = indices.wqi_nsf(agua_env_site.select_dtypes("number").mean(numeric_only=True).to_dict()) if len(agua_env_site) else {"wqi": float("nan"), "clase": "—"}
    phi = indices.polymer_hazard_index(d[d["sitio"] == sitio])
    rows.append({"sitio": sitio, "zona": zona,
                 "mediana_mp": float(g["abundancia"].median()),
                 "clase_carga": indices.mp_tier(g["abundancia"].median(), matriz_f),
                 "phi": round(phi["phi"], 1) if phi["phi"] == phi["phi"] else None,
                 "clase_phi": phi["clase"],
                 "wqi": round(wqi["wqi"], 1) if wqi["wqi"] == wqi["wqi"] else None,
                 "clase_wqi": wqi["clase"]})
resumen = pd.DataFrame(rows).sort_values("mediana_mp", ascending=False)
tabla(resumen)

resumen_box("<b>Lectura rápida:</b> contraste la columna <i>clase_carga</i> contra "
            "<i>wqi</i> y <i>phi</i>: estaciones con carga alta + WQI bajo + PHI alto "
            "concentran la intervención (ver recomendaciones del agente AG-100). "
            "Los umbrales son referenciales y editables en config.py.")
footer()
