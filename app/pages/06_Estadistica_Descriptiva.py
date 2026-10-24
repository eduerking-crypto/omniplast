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
from ui.bootstrap import (chart, footer, page_header, render_result, resumen_box,
                          section_title, tabla)

from omniplast.core.statistics import bayes, descriptive, inferential
from omniplast.db import dao
from omniplast.reporting import figures

page_header("Estadística I — Descriptiva y distribuciones",
            "Diagnóstico de normalidad, homocedasticidad, transformaciones y atípicos")

mp = dao.mp_df()
if mp.empty:
    st.stop()


def merged():
    mpm = mp.groupby("sample_id")["abundancia"].sum().rename("mp_total")
    return dao.env_wide().merge(mpm, on="sample_id", how="left").dropna(subset=["mp_total"])

c1, c2, c3 = st.columns([3, 2, 2])
d = merged()
variables = ["mp_total"] + [c for c in d.select_dtypes("number").columns
                            if c not in ("mp_total", "sample_id", "lat", "lon", "anio",
                                         "profundidad_m") and d[c].notna().sum() >= 10]
variable = c1.selectbox("Variable respuesta", variables, index=0)
grupo = c2.selectbox("Factor de agrupación", ["zona", "temporada", "sitio", "matriz"])
sub = d[[variable, grupo]].dropna()

t_desc, t_norm, t_trans, t_out = st.tabs(["Descriptiva", "Normalidad y supuestos",
                                           "Transformaciones", "Atípicos"])

with t_desc:
    resumen_box("Resumen por grupo (media ± DE, IC95, mediana [IQR]):")
    tab = descriptive.describe_by_group(sub, variable, grupo)
    tabla(tab)
    chart(figures.box(sub, x=grupo, y=variable, title=f"{variable} por {grupo}",
                      color=grupo))
    boot = bayes.bootstrap_stat if hasattr(bayes, "bootstrap_stat") else None
    from omniplast.core.statistics.multilevel import bootstrap_stat
    b = bootstrap_stat(sub[variable], stat="median", n_boot=3000)
    resumen_box(f"<b>Mediana global con IC bootstrap:</b> {b.get('interpretacion', '')}")

with t_norm:
    nt = descriptive.normality_tests(sub[variable])
    tabla(nt)
    lev = descriptive.homogeneity_tests(sub, variable, grupo)
    tabla(lev)
    if not nt.empty and nt.iloc[0]["p"] < .05:
        resumen_box("Distribución <b>no normal</b> detectada → prefiera pruebas no "
                    "paramétricas (Kruskal-Wallis) o GLM para conteos. Use la página de "
                    "Comparaciones o el asesor PLASTI.", "warn")
    else:
        resumen_box("Sin evidencia fuerte de no normalidad; puede usar enfoques "
                    "paramétricos con verificación de homocedasticidad.")

with t_trans:
    st_ = descriptive.suggest_transformation(sub[variable])
    tabla(pd.DataFrame([st_["p_valores"]]).T.rename(columns={0: "p Shapiro"}))
    resumen_box(f"<b>Mejor transformación sugerida: {st_['mejor_transformacion']}</b> — "
                "para conteos ambientales la regla práctica es log(x+1) o sqrt(x); "
                "Box-Cox requiere x>0 en todos los casos.")

with t_out:
    flagged = descriptive.iqr_outliers(sub, variable, grupo)
    out = flagged[flagged["outlier"]] if "outlier" in flagged else pd.DataFrame()
    resumen_box(f"<b>{len(out)}</b> observaciones atípicas por criterio IQR×1.5 (Q1/Q3 "
                "por grupo).")
    tabla(flagged)
    if len(out):
        chart(figures.scatter(flagged, x=grupo, y=variable, color="outlier",
                             title="Atípicos marcados"))

footer()
