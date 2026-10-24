import sys
from pathlib import Path

_here = Path(__file__).resolve()
for _p in _here.parents:
    if (_p / "app" / "ui" / "bootstrap.py").exists():
        for _s in (str(_p), str(_p / "app")):
            if _s not in sys.path:
                sys.path.insert(0, _s)
        break

import numpy as np
import pandas as pd
import streamlit as st
from ui.bootstrap import (chart, exportar_dfs, footer, page_header, resumen_box,
                          section_title, tabla)

from omniplast.core.eco import biomarkers
from omniplast.db import dao
from omniplast.reporting import figures

page_header("Biota y ecotoxicología",
            "Ingestión de MP, biomarcadores oxidativos, condición e histopatología")

bio = dao.bio_df()
mpb = dao.mp_df()
mpb = mpb[mpb["matriz"] == "biota"] if not mpb.empty else mpb

if bio.empty and mpb.empty:
    resumen_box("Sin datos de biota en el proyecto.", "warn")
    st.stop()

t_ing, t_bio, t_cond, t_hist = st.tabs(["Ingestión de MP", "Biomarcadores",
                                        "Índices de condición", "Histopatología"])

with t_ing:
    if mpb.empty:
        resumen_box("Sin registros de MP en biota.", "warn")
    else:
        r = biomarkers.ingestion_metrics(mpb)
        resumen_box(r["interpretacion"])
        tabla(r["resumen"])
        chart(figures.box(r["por_muestra"], x="sitio", y="carga_total",
                          title="MP por individuo y sitio", ytitle="partículas/ind"))
        if "subtipo" in mpb.columns:
            chart(figures.box(r["por_muestra"].merge(
                mpb.drop_duplicates("sample_id")[["sample_id", "subtipo"]],
                on="sample_id", how="left"), x="subtipo", y="carga_total",
                title="Ingestión por especie", ytitle="partículas/ind"))

with t_bio:
    markers = ["sod_u_mg_prot", "cat_u_mg_prot", "gst_u_mg_prot", "lpo_nmol_tba_mg"]
    sub = bio[bio["parametro"].isin(markers)]
    if sub.empty:
        resumen_box("Sin biomarcadores registrados (SOD/CAT/GST/LPO).", "warn")
    else:
        pivot = sub.pivot_table(index=["sample_id", "codigo", "zona", "subtipo"],
                                columns="parametro", values="valor",
                                aggfunc="mean").reset_index()
        tabla(pivot)
        chart(figures.box(pivot, x="zona", y="lpo_nmol_tba_mg",
                          title="Peroxidación lipídica (LPO) por zona",
                          ytitle="nmol TBARS/mg"))
        section_title("IBR/2 integrado por zona (Beliaeff & Burgeot 2002)")
        zonas = pivot["zona"].unique()
        if len(zonas) >= 2:
            ref_zona = zonas[0]
            ref = pivot[pivot["zona"] == ref_zona]
            rows = []
            for zona in zonas:
                if zona == ref_zona:
                    continue
                tr = pivot[pivot["zona"] == zona]
                ibr = biomarkers.ibr_beliaeff(
                    {m: float(tr[m].mean()) for m in markers if m in tr.columns},
                    {m: float(ref[m].mean()) for m in markers if m in ref.columns},
                    {m: float(ref[m].std(ddof=1)) if len(ref) > 1 else 1.0
                     for m in markers if m in ref.columns})
                if "ibr" in ibr:
                    rows.append({"zona": zona, "referencia": ref_zona, "ibr": ibr["ibr"]})
            if rows:
                tabla(pd.DataFrame(rows))
                resumen_box(f"Referencia usada: <b>{ref_zona}</b> (primera zona alfabética). "
                            "Cambie la referencia según su diseño (sitio menos impactado).")
        exportar_dfs("biomarcadores", {"pivot": pivot})

with t_cond:
    sub = bio[bio["parametro"].isin(["factor_condicion", "hsi_pct", "gsi_pct"])]
    if sub.empty:
        resumen_box("Sin CF/HSI/GSI (requiere peces).", "warn")
    else:
        pivot = sub.pivot_table(index=["sample_id", "codigo", "zona", "subtipo"],
                                columns="parametro", values="valor",
                                aggfunc="mean").reset_index()
        tabla(pivot)
        chart(figures.box(pivot, x="zona", y="factor_condicion",
                          title="Factor de condición por zona"))

with t_hist:
    resumen_box("Registre lesiones con: <code>lesion, prevalencia_pct, severidad (0-5), "
                "factor_importancia (1-4)</code>. Ejecute el agente AG-089 para el índice "
                "de Bernet et al. (1999) o cárguelas aquí:", "warn")
    n = st.number_input("Número de lesiones", 3, 15, 5)
    rows = []
    for i in range(n):
        c1, c2, c3, c4 = st.columns([3, 1, 1, 1])
        lesion = c1.text_input(f"Lesión {i+1}", key=f"les{i}")
        prev = c2.number_input("Prev %", 0, 100, 50, key=f"prev{i}")
        sev = c3.number_input("Sev", 0, 5, 2, key=f"sev{i}")
        imp = c4.number_input("Imp", 1, 4, 2, key=f"imp{i}")
        if lesion:
            rows.append({"lesion": lesion, "prevalencia_pct": prev, "severidad": sev,
                         "factor_importancia": imp})
    if st.button("Calcular índice histopatológico", type="primary") and len(rows) >= 3:
        r = biomarkers.histopathology_index(pd.DataFrame(rows))
        resumen_box(f"<b>{r['interpretacion']}</b>")
        tabla(r["tabla"])

footer()
