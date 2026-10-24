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
from omniplast.utils import fmt_num

page_header("Muestras y microplásticos",
            "Explorador de registros: polímero · morfología · color · banda de tamaño · método")

mp = dao.mp_df()
if mp.empty:
    st.stop()

f1, f2, f3, f4 = st.columns(4)
matriz = f1.selectbox("Matriz", sorted(mp["matriz"].unique()))
d = mp[mp["matriz"] == matriz]
sitio = f2.selectbox("Estación", ["Todas"] + sorted(d["sitio"].unique()))
polimero = f3.selectbox("Polímero", ["Todos"] + sorted(d["polimero"].unique()))
morf = f4.selectbox("Morfología", ["Todas"] + sorted(d["morfologia"].unique()))
if sitio != "Todas":
    d = d[d["sitio"] == sitio]
if polimero != "Todos":
    d = d[d["polimero"] == polimero]
if morf != "Todas":
    d = d[d["morfologia"] == morf]

per_sample = d.groupby(["sample_id", "codigo", "sitio", "zona", "campana", "temporada",
                        "profundidad_m"], observed=True)["abundancia"].sum().reset_index()
kpi_val = per_sample["abundancia"]
resumen_box(f"<b>{len(d):,}</b> registros filtrados · mediana <b>{fmt_num(kpi_val.median(), 0)}</b> "
            f"· rango {fmt_num(kpi_val.min(), 0)} – {fmt_num(kpi_val.max(), 0)} · "
            f"clase de carga: <b>{indices.mp_tier(kpi_val.median(), matriz)}</b>")

t_reg, t_comp, t_tam, t_qa = st.tabs(["Registros", "Composición", "Tamaños", "QA/QC"])

with t_reg:
    tabla(per_sample)
    tabla(d.head(300))
    exportar_dfs("registros_mp", {"por_muestra": per_sample, "registros": d})

with t_comp:
    c1, c2 = st.columns(2, gap="large")
    with c1:
        poly = d.groupby("polimero", observed=True)["abundancia"].sum().sort_values(ascending=False)
        chart(figures.pie(list(poly.index), list(poly.values), title="Polímeros"))
    with c2:
        morfs = d.groupby("morfologia", observed=True)["abundancia"].sum().sort_values(ascending=False)
        chart(figures.bars(morfs.reset_index(name="abundancia"), "morfologia", "abundancia",
                           title="Morfologías", ytitle="partículas"))
    section_title("Composición por estación")
    mix = d.groupby(["sitio", "polimero"], observed=True)["abundancia"].sum().reset_index()
    mix["pct"] = mix.groupby("sitio")["abundancia"].transform(lambda s: s / s.sum() * 100)
    chart(figures.stacked(mix, "sitio", "pct", "polimero", title="Firma polimérica por estación",
                          ytitle="% del total"))
    colors = d.groupby("color", observed=True)["abundancia"].sum().sort_values(ascending=False)
    chart(figures.bars(colors.reset_index(name="abundancia"), "color", "abundancia",
                       title="Colores (proxy de foto-degradación)", ytitle="partículas"))

with t_tam:
    bands = d.assign(banda=pd.cut(d["tamano_max_um"], [0, 100, 300, 500, 1000, 5000],
                                  labels=["1-100", "100-300", "300-500", "500-1000",
                                          "1000-5000"]))
    counts = bands.groupby("banda", observed=True)["abundancia"].sum().reset_index()
    counts["pct"] = counts["abundancia"] / counts["abundancia"].sum() * 100
    chart(figures.bars(counts, "banda", "abundancia", title="Distribución de tamaños (µm)",
                       ytitle="partículas"))
    frac = counts.loc[counts["banda"].astype(str).isin(["1-100", "100-300"]), "pct"].sum()
    resumen_box(f"Fracción ≤300 µm: <b>{frac:.1f}%</b> — rango biodisponible para "
                "zooplancton, bivalvos y alevines (riesgo de transferencia trófica).")
    section_title("Tamaño mediano por polímero")
    med = d.assign(um=(d["tamano_min_um"] + d["tamano_max_um"]) / 2).groupby("polimero", observed=True)["um"].median().sort_values()
    chart(figures.bars(med.reset_index(name="um"), "polimero", "um",
                       title="Tamaño mediano por polímero", ytitle="µm"))

with t_qa:
    samples = dao.samples_df()
    qc_col = "qa_qc" if "qa_qc" in samples.columns else None
    if qc_col:
        qc = samples.drop_duplicates("id" if "id" in samples.columns else "sample_id")[qc_col]
        qc = qc.apply(lambda x: x if isinstance(x, dict) else {})
        blancos = sum(1 for q in qc if q.get("blanco_proceso"))
        dups = sum(1 for q in qc if q.get("duplicado"))
        unknown_pct = (d["polimero"] == "UNKNOWN").mean() * 100
        metodos = d["metodo_id"].value_counts()
        resumen_box(f"Blancos de proceso: <b>{blancos}</b> · duplicados de campo: <b>{dups}</b> · "
                    f"polímeros sin identificar: <b>{unknown_pct:.1f}%</b>")
        tabla(metodos.reset_index().rename(columns={"index": "método", "metodo_id": "n"}))
        if unknown_pct > 20:
            resumen_box("Más del 20% de partículas sin identificación espectral: el PHI y el "
                        "análisis de fuentes pierden precisión (Koelmans et al. 2019).", "warn")
    exportar_dfs("qaqc_mp", {"metodos": metodos.reset_index()})

footer()
