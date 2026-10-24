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
from ui.bootstrap import chart, footer, page_header, render_result, resumen_box, section_title, tabla

from omniplast.core.genetics import community, genotyping, primers, qpcr
from omniplast.db import dao
from omniplast.reporting import figures

page_header("Genética & molecular",
            "qPCR (ΔΔCt · eficiencia · fusión) · genotipado · cebadores · catálogo de genes")

q = dao.qpcr_df()
t_expr, t_eff, t_melt, t_genot, t_primers = st.tabs(
    ["Expresión ΔΔCt", "Curva estándar", "Curvas de fusión", "Genotipado", "Cebadores y panel"])

with t_expr:
    if q.empty:
        resumen_box("Sin datos qPCR en el proyecto.", "warn")
    else:
        targets = sorted(set(q["target"]) - {"16s"})
        c1, c2 = st.columns(2)
        gen = c1.selectbox("Gen objetivo", targets, index=0)
        control = c2.selectbox("Zona control (ΔΔCt)", sorted(q["zona"].unique()))
        r = qpcr.ddct_table(q, target=gen, ref="16s", control_group=control)
        if "error" in r:
            resumen_box(r["error"], "warn")
        else:
            g = qpcr.ddct_group_test(r)
            resumen_box(f"{r['interpretacion']} — {g.get('interpretacion', g.get('error', ''))}")
            tabla(r["tabla"].round(3))
            chart(figures.box(r["tabla"], x="zona", y="fc",
                              title=f"2^(-ΔΔCt) — {gen} (control: {control})", ytitle="FC"))
            audit = qpcr.miqe_audit(q)
            section_title("Auditoría MIQE")
            resumen_box(audit["interpretacion"] or "Sin alertas automáticas",
                        "warn" if audit["alertas"] else "ok")
            for a in audit["alertas"]:
                st.markdown(f"- ⚠️ {a}")

with t_eff:
    st.markdown("**Curva estándar:** ingrese 5+ puntos (dilución relativa y Ct). "
                "Ejemplo sintético precargado.")
    cols = st.columns(5)
    defaults_d = [1e6, 1e5, 1e4, 1e3, 1e2]
    defaults_c = [18.2, 21.6, 25.1, 28.4, 31.7]
    vals_d, vals_c = [], []
    for i, col in enumerate(cols):
        vals_d.append(col.number_input(f"Dilución {i+1}", 1.0, 1e8,
                                       float(defaults_d[i]), key=f"d{i}"))
        vals_c.append(col.number_input(f"Ct {i+1}", 5.0, 40.0, float(defaults_c[i]),
                                        step=0.1, key=f"c{i}"))
    if st.button("Calcular eficiencia", type="primary"):
        r = qpcr.standard_curve(vals_d, vals_c)
        resumen_box(r.get("interpretacion", r.get("error", "")),
                    "ok" if r.get("eficiencia", {}).get("valida_miqe") else "warn")
        chart(figures.scatter(pd.DataFrame({"log_copias": np.log10(vals_d), "ct": vals_c}),
                              x="log_copias", y="ct", trendline=True,
                              title="Curva estándar qPCR", xtitle="log₁₀ copias"))

with t_melt:
    temps = st.file_uploader("CSV con columnas: temperatura, fluorescencia", type="csv")
    if temps is not None:
        df = pd.read_csv(temps)
        if {"temperatura", "fluorescencia"} <= set(df.columns):
            r = qpcr.melt_curve_peaks(df["temperatura"], df["fluorescencia"])
            resumen_box(r.get("interpretacion", r.get("error", "")))
            chart(figures.line(pd.DataFrame({
                "temperatura_c": df["temperatura"],
                "dF_dT": -np.gradient(df["fluorescencia"], df["temperatura"])}),
                x="temperatura_c", y="dF_dT", title="-dF/dT", xtitle="°C"))
        else:
            resumen_box("Columnas requeridas: temperatura, fluorescencia.", "error")
    else:
        st.info("Estructura del CSV: `temperatura,fluorescencia` de 65 a 95 °C, "
                "fluorescencia decreciente con mesetas en los Tm.")

with t_genot:
    c1, c2 = st.columns(2)
    n_h1 = c1.number_input("n homocigotos canal X", 5, 100, 24)
    n_h2 = c1.number_input("n homocigotos canal Y", 5, 100, 22)
    n_het = c2.number_input("n heterocigotos", 5, 100, 28)
    if st.button("Llamar genotipos (clusters sintéticos de demostración)", type="primary"):
        rng = np.random.default_rng(7)
        x = np.concatenate([rng.normal(.9, .08, int(n_h1)), rng.normal(.1, .06, int(n_h2)),
                           rng.normal(.5, .07, int(n_het))])
        y = np.concatenate([rng.normal(.1, .06, int(n_h1)), rng.normal(.9, .08, int(n_h2)),
                           rng.normal(.55, .07, int(n_het))])
        r = genotyping.call_genotypes(x, y)
        resumen_box(r["interpretacion"])
        chart(figures.scatter(r["llamadas"], x="x", y="y", color="genotipo",
                              title="Discriminación alélica por clustering"))
        c = r["conteo"]
        if {"homocigoto_X", "heterocigoto", "homocigoto_Y"} <= set(c.index):
            hwe = genotyping.hardy_weinberg_test(
                [int(c["homocigoto_X"]), int(c["heterocigoto"]), int(c["homocigoto_Y"])])
            resumen_box(f"<b>Hardy-Weinberg:</b> {hwe['interpretacion']}")
        st.caption("Genotipado de demostración: cargue sus canales FAM/HEX reales por CSV "
                   "en producción (columnas x, y).")

with t_primers:
    section_title("Panel de genes de degradación (referencias primarias)")
    panel, disclaimer = primers.gene_panel_info()
    tabla(pd.DataFrame(panel))
    resumen_box(disclaimer, "warn")
    section_title("Cebadores universales clásicos (secuencias publicadas)")
    tabla(pd.DataFrame(primers.universal_primers_table()))
    seq = st.text_area("Diseño de cebadores: pegue una secuencia objetivo (≥120 nt)",
                       height=110)
    if len(seq.strip()) >= 120 and st.button("Proponer cebadores", type="primary"):
        r = primers.design_primers_from_sequence(seq)
        if r["candidatos"]:
            tabla(pd.DataFrame(r["candidatos"]).head(10))
        else:
            resumen_box("Sin parejas que cumplan Tm 57-63 °C / GC 45-60% / sin dímeros "
                        "3'. Relaje los criterios o alargue la secuencia.", "warn")

footer()
