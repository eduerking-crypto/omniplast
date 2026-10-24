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
from ui.bootstrap import chart, footer, page_header, resumen_box, section_title, tabla

from omniplast.core.statistics import bayes, timeseries
from omniplast.db import dao
from omniplast.reporting import figures

page_header("Estadística V — Series temporales",
            "Mann-Kendall · pendiente de Sen · MK estacional · STL · Holt · ACF")

mp = dao.mp_df()
if mp.empty:
    st.stop()

d = mp[mp["matriz"] == "agua"].groupby(["sitio", "fecha", "temporada", "campana"],
                                       observed=True)["abundancia"].sum().reset_index()

c1, c2 = st.columns([2, 2])
sitio = c1.selectbox("Estación", ["Todo el sistema"] + sorted(d["sitio"].unique()))
if sitio != "Todo el sistema":
    d = d[d["sitio"] == sitio]
serie = d.groupby("fecha")["abundancia"].mean().sort_index().reset_index()
resumen_box(f"Serie con <b>{len(serie)}</b> puntos temporales "
            f"({serie['fecha'].min():%Y-%m} a {serie['fecha'].max():%Y-%m}).")

if len(serie) >= 6:
    t1, t2 = st.tabs(["Tendencia", "Descomposición y pronóstico"])

    with t1:
        mk = timeseries.mann_kendall(serie["abundancia"])
        resumen_box(f"<b>Mann-Kendall:</b> {mk.get('interpretacion', mk.get('error', ''))}")
        sen = timeseries.sens_slope(serie["abundancia"])
        resumen_box(f"<b>Pendiente de Sen:</b> {sen.get('interpretacion', sen.get('error', ''))}")
        mks = timeseries.seasonal_mann_kendall(d.assign(fecha2=d["fecha"]),
                                                "fecha2", "abundancia")
        if "por_estacion" in mks:
            section_title("Mann-Kendall estacional (por mes de muestreo)")
            tabla(mks["por_estacion"].round(3))
        chart(figures.line(serie.assign(fecha=serie["fecha"].dt.strftime("%Y-%m")),
                           x="fecha", y="abundancia",
                           title="Serie temporal de MP (media por campaña)",
                           ytitle="part/m³"))

    with t2:
        stl = timeseries.stl_decompose(d, "fecha", "abundancia", period=2)
        if "componentes" in stl:
            comp = stl["componentes"].assign(
                fecha=lambda x: pd.to_datetime(x["fecha"]).dt.strftime("%Y-%m"))
            chart(figures.line(comp, x="fecha", y="observado",
                               title="Observado vs tendencia",
                               color=pd.Series(["observado"] * len(comp))))
            tabla(comp.head(12))
        else:
            resumen_box(f"STL: {stl.get('error', '')} — con campañas semestrales el "
                        "periodo estacional debe ≥2; con 4 campañas use Mann-Kendall "
                        "estacional.", "warn")
        acf = timeseries.acf_pacf(serie["abundancia"])
        section_title("Autocorrelación")
        chart(figures.bars(acf["acf"], "lag", "acf", title="ACF (IC Bartlett ±"
                         f"{acf['ic_bartlett']:.2f})", ytitle="ACF"))
        holt = timeseries.holt_forecast(serie["abundancia"], periods=4)
        if "pronostico" in holt:
            section_title("Pronóstico Holt amortiguado (4 pasos)")
            fut_idx = pd.DataFrame({"paso": range(1, 5), "pronostico": holt["pronostico"]})
            chart(figures.line(fut_idx, x="paso", y="pronostico",
                               title="Proyección de la media de MP", ytitle="part/m³"))
else:
    resumen_box("La serie requiere ≥6 puntos temporales. Añada campañas o reduzca el "
                "filtro de estación.", "warn")

footer()
