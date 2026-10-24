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
from ui.bootstrap import chart, editor_seguro, footer, page_header, resumen_box, section_title, tabla

from omniplast.core import indices
from omniplast.core.eco import onehealth, toxicology
from omniplast.db import dao
from omniplast.reporting import figures

page_header("Riesgo ecológico & One Health",
            "Cocientes HQ (PEC/PNEC) · SSD y HC5 · exposición humana por consumo de pesca")

mp = dao.mp_df()
if mp.empty:
    st.stop()

t_hq, t_ssd, t_oh, t_social = st.tabs(["HQ por polímero", "SSD (HC5)",
                                        "Exposición humana", "Percepción social"])

with t_hq:
    mpa = mp[mp["matriz"] == "agua"]
    pec = (mpa.groupby("polimero", observed=True)["abundancia"].mean() / 1000.0).to_dict()
    r = toxicology.hazard_quotient_table(pec)
    resumen_box(f"<b>{r['interpretacion']}</b><br>{r['nota']}")
    tabla(r["tabla"])
    if r["tabla"]["hq"].notna().any():
        chart(figures.bars(r["tabla"].dropna(subset=["hq"]), "polimero", "hq",
                           title="HQ = PEC/PNEC (PNEC p/L referencial)", ytitle="HQ"))
    resumen_box("<b>Marco:</b> EU TGD-R.10: HQ<1 riesgo aceptable; 1-10 medio; >10 alto. "
                "Los PNEC por defecto son de literatura y editables — NO son valores "
                "regulatorios oficiales para microplásticos.", "warn")

with t_ssd:
    st.markdown("La **Distribution of Species Sensitivities** requiere ≥5 valores de "
                "toxicidad (NOEC/EC10/LC10) de especies distintas, en la MISMA unidad. "
                "Cárguela o edite la tabla de ejemplo (valores ilustrativos).")
    base = pd.DataFrame({
        "especie": ["Daphnia magna", "Artemia sp.", "Danio rerio", "Mytilus edulis",
                    "Skeletonema sp.", "Tigriopus japonicus", "Oryzias latipes"],
        "conc_mg_L": [1.2, 5.8, 8.1, 2.6, 0.9, 4.4, 12.5]})
    edited = editor_seguro(base, num_rows="dynamic")
    if st.button("Ajustar SSD lognormal y calcular HC5", type="primary"):
        data = dict(zip(edited["especie"], pd.to_numeric(edited["conc_mg_L"],
                                                         errors="coerce")))
        r = toxicology.ssd_fit(data)
        if "error" in r:
            resumen_box(r["error"], "warn")
        else:
            resumen_box(f"<b>HC5 = {r['hc5']:.3g}</b> — concentra por debajo de la cual "
                        "se protege al 95% de las especies. Compare contra sus PEC y "
                        "use PAF(PEC) para fracción afectada.")
            curve = pd.DataFrame({"concentracion": r["curve_x"], "PAF_pct": r["curve_y"]})
            chart(figures.line(curve, x="concentracion", y="PAF_pct",
                               title="Curva SSD (PAF acumulada)", xtitle="concentración",
                               ytitle="% especies afectadas"))
            pec_input = st.number_input("PEC a evaluar (misma unidad)", 0.0, 100.0, 5.0)
            st.metric("PAF en su PEC", f"{r['paf_func'](pec_input):.1f}%")

with t_oh:
    mpb = mp[mp["matriz"] == "biota"]
    if mpb.empty:
        resumen_box("Sin MP en biota para estimar exposición.", "warn")
    else:
        per = mpb.groupby(["sample_id", "subtipo"], observed=True)["abundancia"].sum()
        c1, c2 = st.columns(2)
        g_ostion = c1.number_input("Porción ostra (g)", 5, 60, 15)
        g_pez = c2.number_input("Porción pescado (g)", 50, 300, 120)
        conc = []
        for row in per.reset_index().itertuples():
            clave = "ostra" if "ostra" in str(row.subtipo).lower() else "pez"
            g = g_ostion if clave == "ostra" else g_pez
            conc.append(max(0.0, row.abundancia / max(g, 1e-9)))
        r = onehealth.exposure_seafood(np.array(conc),
                                       meal_g=(80, 150, 300),
                                       meals_per_week=(1, 3, 7))
        resumen_box(r.get("interpretacion", r.get("error", "")))
        if "sensibilidad" in r:
            section_title("Sensibilidad de la estimación (Spearman)")
            tabla(r["sensibilidad"])
        st.caption("Advertencia metodológica: la OMS (2019) concluye evidencia "
                   "insuficiente para límites de seguridad; esta estimación compara "
                   "ESCENARIOS entre sí (p. ej., antes/después de una intervención).")

with t_social:
    st.markdown("Cargue un CSV de encuesta Likert (columnas `p1..pn` en escala 1-5, "
                "opcional `estrato`). El sistema calcula alfa de Cronbach y comparaciones "
                "por grupo. También puede probar con datos de demostración.")
    up = st.file_uploader("Encuesta (CSV)", type="csv")
    if up is not None:
        df = pd.read_csv(up)
    else:
        rng = np.random.default_rng(11)
        df = pd.DataFrame({f"p{i}": rng.integers(1, 6, 24) for i in range(1, 6)})
        df["estrato"] = rng.choice(["pescador", "turista", "residente"], 24)
    likert = [c for c in df.columns if c.startswith("p")]
    if len(likert) >= 3:
        r = onehealth.social_survey_analysis(df, likert,
                                             group_col="estrato" if "estrato" in df.columns else None)
        resumen_box(r.get("interpretacion", r.get("error", "")))
        tabla(r.get("resumen", pd.DataFrame()))
        if "por_grupo" in r and len(r["por_grupo"]):
            section_title("Diferencias por estrato (Kruskal-Wallis)")
            tabla(r["por_grupo"])
    else:
        resumen_box("Se requieren ≥3 ítems Likert (p1..pn).", "warn")

footer()
