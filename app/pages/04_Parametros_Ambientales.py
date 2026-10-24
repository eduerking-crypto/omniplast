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

page_header("Parámetros ambientales",
            "Fisicoquímica, nutrientes y microbiología del agua y sedimento")

env = dao.env_wide()
if env.empty:
    st.stop()

agua = env[env["matriz"] == "agua"]
sed = env[env["matriz"] == "sedimento"]
f1, f2 = st.columns([2, 2])
vista = f1.selectbox("Vista", ["Agua", "Sedimento", "Todas"])
df_v = {"Agua": agua, "Sedimento": sed, "Todas": env}[vista]
sitio_f = f2.selectbox("Estación", ["Todas"] + sorted(df_v["sitio"].unique()))
if sitio_f != "Todas":
    df_v = df_v[df_v["sitio"] == sitio_f]

t_datos, t_wqi, t_tsi, t_corr = st.tabs(["Datos", "WQI-NSF", "Estado trófico (TSI)",
                                          "Correlaciones"])

with t_datos:
    tabla(df_v)
    exportar_dfs("parametros_ambientales", {vista: df_v})

with t_wqi:
    section_title("Índice de calidad del agua NSF por estación")
    rows = []
    for sitio, g in df_v.groupby("sitio"):
        w = indices.wqi_nsf(g.select_dtypes("number").mean(numeric_only=True).to_dict())
        rows.append({"sitio": sitio, "wqi": w["wqi"], "clase": w["clase"],
                     "n_parametros": w.get("n_parametros", 0)})
    tab = pd.DataFrame(rows).dropna().sort_values("wqi")
    chart(figures.bars(tab, "sitio", "wqi", title="WQI-NSF por estación (0-100)",
                       ytitle="WQI"))
    tabla(tab)
    resumen_box("Fórmula: WQI = Σ(qᵢ·wᵢ)/Σwᵢ con curvas de valoración NSF aproximadas "
                "(Brown et al. 1970). Clases: >90 excelente · 71-90 buena · 51-70 media · "
                "26-50 mala · <26 muy mala.")
    tabla(tab.set_index("sitio")[["wqi"]].T)

with t_tsi:
    section_title("Estado trófico TSI (Carlson 1977)")
    rows = []
    for sitio, g in df_v.groupby("sitio"):
        chla = g["clorofila_a_ug_l"].median() if "clorofila_a_ug_l" in g else None
        tp = g["fosfato_mg_l"].median() * 3.066 * 1000 if "fosfato_mg_l" in g else None
        t = indices.tsi_carlson(chla=chla, tp_ug=tp)
        rows.append({"sitio": sitio, "tsi": t["tsi"], "clase": t["clase"],
                     "por_variable": str(t.get("por_variable", ""))})
    tab = pd.DataFrame(rows).dropna(subset=["tsi"]).sort_values("tsi")
    chart(figures.bars(tab, "sitio", "tsi", title="TSI por estación", ytitle="TSI"))
    resumen_box("TSI(clorofila) = 30.6 + 9.81·ln(chla); TSI(PT) = 14.42·ln(PT_µg/L) + 4.15. "
                "Oligotrófico <30 · Mesotrófico 40-50 · Eutrófico 60-70 · Hipereutrófico >70.")

with t_corr:
    num_cols = [c for c in df_v.select_dtypes("number").columns
                if c not in ("sample_id", "lat", "lon", "profundidad_m", "anio")
                and df_v[c].notna().sum() >= 10]
    if len(num_cols) >= 3:
        from omniplast.core.statistics import inferential
        r_mat, p_mat = inferential.correlation_matrix(df_v[num_cols[:14]], method="spearman")
        chart(figures.heatmap(r_mat, title="Correlaciones ambientales (Spearman)"))
        section_title("Variables más asociadas entre sí (p<0.05)")
        pares = []
        cols = r_mat.columns
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                if p_mat.iloc[i, j] < 0.05:
                    pares.append({"par_1": cols[i], "par_2": cols[j],
                                  "r": r_mat.iloc[i, j], "p": p_mat.iloc[i, j]})
        tabla(pd.DataFrame(pares).sort_values("p").head(25))
    else:
        resumen_box("Se requieren ≥3 variables numéricas con datos suficientes.", "warn")

footer()
