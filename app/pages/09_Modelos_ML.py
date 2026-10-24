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

from omniplast.core.statistics import ml, multilevel
from omniplast.db import dao
from omniplast.reporting import figures

page_header("Estadística IV — Modelos, ML y potencia",
            "OLS · GLM Poisson/NB · mixtos · GEE · cuantílica · RF · boosting · MCMC · potencia")

mp = dao.mp_df()
if mp.empty:
    st.stop()

mpm = mp.groupby("sample_id")["abundancia"].sum().rename("mp_total")
d = dao.env_wide().merge(mpm, on="sample_id", how="left").dropna(subset=["mp_total"])
numvars = [c for c in d.select_dtypes("number").columns
            if c not in ("mp_total", "sample_id", "lat", "lon", "anio", "profundidad_m")
            and d[c].notna().sum() >= 10]

t_reg, t_glm, t_mixtos, t_ml, t_bayes = st.tabs(
    ["Regresiones", "GLM de conteos", "Mixtos y GEE", "Machine Learning", "Bayes y potencia"])

with t_reg:
    sel = st.multiselect("Predictores", numvars, default=numvars[:4], key="ols")
    if len(sel) >= 2:
        r = multilevel.ols_regression(d.dropna(subset=sel + ["mp_total"]), "mp_total", sel)
        resumen_box(r.get("interpretacion", r.get("error", "")))
        if "coeficientes" in r:
            tabla(r["coeficientes"].round(4))
            tabla(r["vif"].to_frame("VIF"))
        qr = multilevel.quantile_regression(d.dropna(subset=sel[:1] + ["mp_total"]),
                                            "mp_total", sel[0])
        if "tabla" in qr:
            section_title(f"Regresión cuantílica ({sel[0]})")
            tabla(qr["tabla"])

with t_glm:
    sel2 = st.multiselect("Predictores del GLM", numvars, default=numvars[:2], key="glm")
    familia = st.radio("Familia", ["Binomial Negativa (recomendada)", "Poisson"], horizontal=True)
    fam = "negativebinomial" if "Negativa" in familia else "poisson"
    if len(sel2) >= 1:
        r = multilevel.glm_counts(d.dropna(subset=sel2 + ["mp_total"]), "mp_total",
                                 sel2, family=fam)
        resumen_box(r.get("interpretacion", r.get("error", "")))
        if "coeficientes" in r:
            tabla(r["coeficientes"].round(4))
            tabla(r["razones_de_cambio"].to_frame("exp(β) × unidad").round(3))

with t_mixtos:
    sel3 = st.selectbox("Predictor del modelo mixto", numvars, index=0)
    r = multilevel.mixedlm(d, "mp_total", [sel3], "sitio")
    resumen_box(r.get("interpretacion", r.get("error", "")))
    if "efectos_fijos" in r:
        tabla(r["efectos_fijos"].round(4))
    g = multilevel.gee_population(d, "mp_total", [sel3], "sitio", family="gaussian")
    section_title("GEE (SE robustas por sitio)")
    if "coeficientes" in g:
        tabla(g["coeficientes"].round(4))
    med = multilevel.mediation_analysis(d, sel3, "dbo_mg_l", "mp_total", n_boot=1500)
    section_title("Mediación (bootstrap)")
    resumen_box(med.get("interpretacion", med.get("error", "")))

with t_ml:
    sel4 = st.multiselect("Features para ML", numvars, default=numvars[:5], key="ml")
    if len(sel4) >= 3:
        dd = d.copy()
        dd["impacto"] = (dd["zona"] == "impacto").astype(int)
        rfc = ml.random_forest_classifier(dd, sel4, "impacto")
        resumen_box(f"<b>RF clasificador (zona impacto):</b> "
                    f"{rfc.get('interpretacion', rfc.get('error', ''))}")
        if "importancias" in rfc:
            chart(figures.bars(rfc["importancias"].head(8).rename("importancia")
                               .rename_axis("variable").reset_index(),
                               "variable", "importancia", horizontal=True,
                               title="Importancia por permutación"))
        rfr = ml.random_forest_regressor(d, sel4, "mp_total")
        section_title("RF regresor de MP")
        resumen_box(rfr.get("interpretacion", rfr.get("error", "")))
        if "matriz_confusion" in rfc:
            tabla(rfc["matriz_confusion"])
        scan = ml.silhouette_scan(d.dropna(subset=sel4), sel4)
        section_title("k-means (estaciones)")
        chart(figures.bars(scan, "k", "silueta", title="Selección de k por silueta"))

with t_bayes:
    sub = d.dropna(subset=["mp_total"])
    zonas = sorted(sub["zona"].unique())
    if len(zonas) >= 2:
        from omniplast.core.statistics import bayes
        a = sub[sub["zona"] == zonas[0]]["mp_total"].values
        b = sub[sub["zona"] == zonas[1]]["mp_total"].values
        r = bayes.mcmc_two_means(a, b, n_iter=8000)
        resumen_box(f"<b>MCMC (diferencia {zonas[0]} − {zonas[1]}):</b> "
                    f"{r.get('interpretacion', r.get('error', ''))}")
        post = bayes.normal_posterior(sub["mp_total"].values)
        section_title("Posterior de la media (conjugado)")
        resumen_box(post.get("interpretacion", post.get("error", "")))
    section_title("Análisis de potencia a priori")
    col1, col2 = st.columns([2, 3])
    d_esperado = col1.slider("d de Cohen esperado", 0.2, 1.5, 0.8, 0.05)
    potencia = col2.slider("Potencia objetivo", 0.6, 0.95, 0.8, 0.05)
    pw = bayes.power_analysis(effect_d=d_esperado, power=potencia)
    resumen_box(pw.get("interpretacion", pw.get("error", "")))

footer()
