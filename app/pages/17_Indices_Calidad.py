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
from ui.bootstrap import chart, footer, kpi_row, page_header, resumen_box, section_title, tabla

from omniplast.config import QUALITY_COLORS
from omniplast.core import indices
from omniplast.db import dao
from omniplast.reporting import figures, report

page_header("Índices de calidad del ecosistema",
            "PLI · CF · PHI · WQI-NSF · TSI · **IICAM** (índice compuesto propuesto) con Monte Carlo")

mp = dao.mp_df()
env = dao.env_wide()
if mp.empty or env.empty:
    st.stop()

matriz = st.selectbox("Matriz para índices de MP", ["sedimento", "agua", "biota"])
mpm = mp[mp["matriz"] == matriz]
por_sitio = mpm.groupby("sitio")["abundancia"].sum()
fondo = por_sitio.min()

rows = []
agua = env[env["matriz"] == "agua"]
for sitio in sorted(por_sitio.index):
    total = por_sitio[sitio]
    cf = indices.contamination_factor(total, fondo)
    pli = indices.pli_site([cf])
    phi = indices.polymer_hazard_index(mpm[mpm["sitio"] == sitio])
    g = agua[agua["sitio"] == sitio]
    wqi = indices.wqi_nsf(g.select_dtypes("number").mean(numeric_only=True).to_dict()) if len(g) else {"wqi": np.nan, "clase": "—"}
    rows.append({"sitio": sitio, "zona": mpm[mpm["sitio"] == sitio]["zona"].iloc[0],
                 "mp_total": float(total), "cf": cf,
                 "cf_clase": indices.cf_classify(cf) if cf == cf else "—",
                 "pli": pli if pli == pli else None,
                 "phi": round(phi["phi"], 1) if phi["phi"] == phi["phi"] else None,
                 "phi_clase": phi["clase"],
                 "wqi": round(wqi["wqi"], 1) if wqi["wqi"] == wqi["wqi"] else None,
                 "wqi_clase": wqi["clase"]})
tab = pd.DataFrame(rows)

kpi_row([
    (f"{np.exp(np.mean(np.log(tab['pli'].clip(lower=1e-9)))):.2f}", "PLI del sistema",
     ">1 = deterioro (Tomlinson 1980)", "#0A2540"),
    (f"{tab['phi'].dropna().median():.0f}" if tab["phi"].notna().any() else "—",
     "PHI mediano", "peligro polimérico (Lithner 2011)", "#E9B44C"),
    (f"{tab['wqi'].dropna().mean():.0f}" if tab["wqi"].notna().any() else "—",
     "WQI-NSF medio", "calidad fisicoquímica", "#05BFDB"),
    (f"{fondo:,.0f}", "Fondo de referencia", "estación mínima del sistema", "#7C9EB2"),
])

t_tabla, t_iicam, t_reporte = st.tabs(["Índices por estación", "IICAM compuesto",
                                       "Informe ejecutivo (DOCX)"])

with t_tabla:
    tabla(tab)
    chart(figures.bars(tab, "sitio", "cf", title="Factor de contaminación por estación",
                       ytitle="CF"))
    chart(figures.bars(tab.dropna(subset=["phi"]).sort_values("phi"), "sitio", "phi",
                       title="Polymer Hazard Index por estación", ytitle="PHI"))
    resumen_box("<b>Cadena de citado:</b> CF (Hakanson 1980) → PLI (Tomlinson 1980) · "
                "PHI (Lithner 2011 + abundancias) · WQI (Brown 1970) · "
                "clases editables en <code>omniplast/config.py</code>.")

with t_iicam:
    st.markdown("**IICAM — Índice Integrado de Calidad Ambiental por Microplásticos** "
                "(propuesto en OMNIPLAST; metodología completa en docs/METODOLOGIA.md). "
                "Compone: carga de MP (percentil invertido) · PHI (log-normalizado) · WQI · "
                "subíndices opcionales de biota. Salida 0-100 (100 = excelente) con banda "
                "de incertidumbre Monte Carlo.")
    w1, w2, w3 = st.columns(3)
    peso_mp = w1.slider("Peso carga MP", 0.0, 1.0, 0.30, 0.05)
    peso_phi = w2.slider("Peso PHI", 0.0, 1.0, 0.20, 0.05)
    peso_wqi = w3.slider("Peso WQI", 0.0, 1.0, 0.20, 0.05)
    resultados = {}
    for sitio in sorted(por_sitio.index):
        row = tab[tab["sitio"] == sitio].iloc[0]
        r = indices.iicam({"mp_total": float(row["mp_total"])}, por_sitio.values,
                          phi_site=row["phi"], wqi_site=row["wqi"],
                          weights={"mp": peso_mp, "phi": peso_phi, "wqi": peso_wqi,
                                   "bio": 0.15, "tox": 0.15})
        resultados[sitio] = r
    rows_i = [{"sitio": s, "iicam": r["iicam"], "clase": r["clase"],
               "ic_inf": r["ic95_mc"][0], "ic_sup": r["ic95_mc"][1],
               "componentes": ", ".join(r["componentes_activos"])}
              for s, r in resultados.items()]
    tab_i = pd.DataFrame(rows_i).sort_values("iicam")
    chart(figures.bars(tab_i, "sitio", "iicam", title="IICAM por estación (0-100)",
                       ytitle="calidad"))
    for s, r in resultados.items():
        col = QUALITY_COLORS.get(r["clase"], "#00A896")
        chart(figures.gauge(r["iicam"], title=f"{s} — {r['clase']}", height=210))
    tabla(tab_i)
    resumen_box("El IICAM es una propuesta metodológica de OMNIPLASTAPLATAFORMA: valide "
                "ponderaciones con panel de expertos y análisis de sensibilidad antes de "
                "uso formal. La banda IC proviene de perturbación Monte Carlo (σ=5).", "warn")

with t_reporte:
    if st.button("Generar informe DOCX de calidad", type="primary"):
        secciones = [{"titulo": "Radiografía por estación",
                      "tablas": {"indices": tab},
                      "parrafos": ["Corte transversal de índices de contaminación por "
                                   "microplásticos (CF, PLI, PHI), calidad fisicoquímica "
                                   "(WQI-NSF) y clase de carga."]},
                     {"titulo": "IICAM",
                      "tablas": {"iicam": tab_i},
                      "parrafos": ["Índice compuesto propuesto con incertidumbre Monte "
                                   "Carlo. Clases: Excelente ≥80, Buena 60-79, Aceptable "
                                   "45-59, Deteriorada 30-44, Crítica <30."]}]
        path = report.generate_docx("Índices de calidad frente a microplásticos",
                                     secciones,
                                     meta={"proyecto": "CLV-DEMO (datos sintéticos)",
                                           "region": "Ciénaga de La Virgen, Cartagena",
                                           "advertencia": "DATOS SINTÉTICOS — demostración.",
                                           "referencias": [
                                               "Hakanson (1980) Hydrobiologia 91:832-848",
                                               "Tomlinson et al. (1980)",
                                               "Lithner et al. (2011) ES&T 45:4727-4735",
                                               "Brown et al. (1970) WQI-NSF",
                                               "Carlson (1977) Limnol. Oceanogr."]})
        with open(path, "rb") as f:
            st.download_button("Descargar informe DOCX", f.read(),
                               file_name=path.name,
                               mime="application/vnd.openxmlformats-officedocument."
                                    "wordprocessingml.document")
        st.success(f"Informe generado: {path.name}")

footer()
