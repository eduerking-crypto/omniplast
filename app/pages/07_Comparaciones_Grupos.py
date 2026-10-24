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

from omniplast.core.statistics import inferential
from omniplast.db import dao
from omniplast.reporting import figures

page_header("Estadística II — Comparaciones de grupos",
            "t/Welch · ANOVA+Tukey · Kruskal+Dunn · Games-Howell · pareadas · χ² · correlaciones")

mp = dao.mp_df()
if mp.empty:
    st.stop()

mpm = mp.groupby("sample_id")["abundancia"].sum().rename("mp_total")
d = dao.env_wide().merge(mpm, on="sample_id", how="left").dropna(subset=["mp_total"])

c1, c2 = st.columns([3, 2])
numvars = [c for c in d.select_dtypes("number").columns
            if c not in ("sample_id", "lat", "lon", "anio", "profundidad_m")
            and d[c].notna().sum() >= 10]
variable = c1.selectbox("Variable", ["mp_total"] + numvars, index=0)
grupo = c2.selectbox("Factor", ["zona", "temporada", "sitio"])
sub = d[[variable, grupo]].dropna()

t_t, t_anova, t_kw, t_pareadas, t_corr = st.tabs(
    ["t de Welch", "ANOVA + Tukey", "Kruskal-Wallis + Dunn", "Pareadas",
     "Correlaciones"])

with t_t:
    zonas = sorted(sub[grupo].unique())
    if len(zonas) >= 2:
        cc1, cc2 = st.columns(2)
        g1 = cc1.selectbox("Grupo 1", zonas, index=0)
        g2 = cc2.selectbox("Grupo 2", zonas, index=len(zonas) - 1)
        r = inferential.ttest_two(sub, variable, grupo, g1, g2)
        resumen_box(r.get("interpretacion", r.get("error", "")))
        chart(figures.box(sub[sub[grupo].isin([g1, g2])], x=grupo, y=variable,
                          title=f"{variable}: {g1} vs {g2}"))

with t_anova:
    r = inferential.anova_one(sub, variable, grupo)
    resumen_box(r.get("interpretacion", r.get("error", "")))
    if "error" not in r:
        tk = inferential.tukey_hsd(sub, variable, grupo)
        if len(tk):
            tabla(tk)
        gh = inferential.games_howell(sub, variable, grupo)
        if len(gh):
            section_title("Games-Howell (robusto a heterocedasticidad)")
            tabla(gh)

with t_kw:
    r = inferential.kruskal_wallis(sub, variable, grupo)
    resumen_box(r.get("interpretacion", r.get("error", "")))
    dunn = inferential.dunn_posthoc(sub, variable, grupo)
    if len(dunn):
        tabla(dunn)

with t_pareadas:
    per = mp[mp["matriz"] == "agua"].groupby(["sitio", "temporada"],
                                             observed=True)["abundancia"].sum().reset_index()
    piv = per.pivot_table(index="sitio", columns="temporada", values="abundancia").dropna()
    if piv.shape[1] >= 2 and len(piv) >= 3:
        cols = sorted(piv.columns)
        melted = piv.reset_index().melt(id_vars="sitio", value_vars=cols,
                                        var_name="temporada", value_name="abundancia")
        w = inferential.wilcoxon_paired(melted, "abundancia", "temporada", cols[0], cols[1])
        resumen_box(w.get("interpretacion", w.get("error", "")))
        tabla(piv)
    else:
        resumen_box("Se requieren ≥3 sitios medidos en ambas temporadas.", "warn")

with t_corr:
    sel = st.multiselect("Variables para correlación", numvars,
                         default=numvars[:6])
    if len(sel) >= 3:
        r_mat, p_mat = inferential.correlation_matrix(d[sel + ["mp_total"]].dropna(),
                                                      method="spearman")
        chart(figures.heatmap(r_mat, title="Spearman — incluida MP total"))
        pares = []
        for i in range(len(r_mat.columns)):
            for j in range(i + 1, len(r_mat.columns)):
                pares.append({"par": f"{r_mat.columns[i]} × {r_mat.columns[j]}",
                              "r": r_mat.iloc[i, j], "p": p_mat.iloc[i, j]})
        tabla(pd.DataFrame(pares).sort_values("p").head(20))
    elif len(sel) == 2:
        r = inferential.correlation_test(d[sel[0]], d[sel[1]], method="spearman")
        resumen_box(r.get("interpretacion", ""))
        chart(figures.scatter(d.dropna(subset=sel), x=sel[0], y=sel[1], trendline=True,
                              title="Dispersión con ajuste lineal"))
    else:
        resumen_box("Seleccione ≥2 variables.", "warn")

footer()
