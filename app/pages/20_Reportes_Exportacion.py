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
from ui.bootstrap import footer, page_header, resumen_box, section_title, tabla

from omniplast.db import dao
from omniplast.reporting import r_scripts, report

page_header("Reportes & exportación",
            "Informe DOCX académico · scripts R para RStudio (vegan/DESeq2/lavaan/glmmTMB) "
            "· bitácora de resultados")

t_docx, t_r, t_historial = st.tabs(["Informe integrado (DOCX)",
                                    "Scripts R (RStudio)", "Resultados guardados"])

with t_docx:
    st.markdown("El informe ejecutivo compone: resumen del sistema, tabla de índices por "
                "estación, composición de MP, IICAM y descargo metodológico. Use la página "
                "<i>Índices de calidad</i> para informes más específicos.")
    titulo = st.text_input("Título del informe", "Diagnóstico integral de microplásticos "
                           "— Ciénaga de La Virgen (demo)")
    if st.button("Generar informe DOCX", type="primary"):
        from omniplast.core import indices
        mp = dao.mp_df()
        env = dao.env_wide()
        por_sitio = mp[mp["matriz"] == "agua"].groupby("sitio")["abundancia"].median()
        rows = []
        agua = env[env["matriz"] == "agua"]
        for sitio, med in por_sitio.items():
            g = agua[agua["sitio"] == sitio]
            w = indices.wqi_nsf(g.select_dtypes("number").mean(numeric_only=True).to_dict()) if len(g) else {"wqi": float("nan"), "clase": "—"}
            rows.append({"sitio": sitio, "mp_mediana": float(med),
                         "clase_carga": indices.mp_tier(med, "agua"),
                         "wqi": w["wqi"], "clase_wqi": w["clase"]})
        tab = pd.DataFrame(rows)
        comp = mp.groupby("polimero", observed=True)["abundancia"].sum().sort_values(ascending=False).head(10)
        secciones = [
            {"titulo": "Carga por estación", "tablas": {"indices": tab}},
            {"titulo": "Composición polimérica del sistema",
             "tablas": {"polimeros": comp.reset_index(name="abundancia")}},
        ]
        path = report.generate_docx(titulo, secciones,
                                    meta={"proyecto": "CLV-DEMO (DATOS SINTÉTICOS)",
                                          "region": "Ciénaga de La Virgen, Cartagena",
                                          "advertencia": "DATOS SINTÉTICOS de demostración: "
                                          "no citar como evidencia ambiental.",
                                          "referencias": [
                                              "Hartmann et al. (2019) ES&T 53:1039-1047",
                                              "Koelmans et al. (2019) ES&T 53:9615-9623",
                                              "Lithner et al. (2011) ES&T 45:4727-4735",
                                              "Brown et al. (1970) WQI-NSF"]})
        with open(path, "rb") as f:
            st.download_button("Descargar DOCX", f.read(), file_name=path.name,
                               mime="application/vnd.openxmlformats-officedocument."
                                    "wordprocessingml.document")
        st.success(f"Informe generado en {path.parent}")

with t_r:
    st.markdown("OMNIPLAST genera los scripts canónicos de R para que el flujo doctoral "
                "pueda reproducirse en **RStudio** (vegan, DESeq2, lavaan, glmmTMB). "
                "Ejecútelos con `Rscript` o ábralos en RStudio.")
    if st.button("Generar los 4 scripts R", type="primary"):
        r = r_scripts.generate_all()
        for name, path in r["scripts"].items():
            st.markdown(f"**{name}** — `{path}`")
            with open(path, "rb", encoding=None) if False else open(path, "rb") as f:
                st.download_button(f"Descargar {name}", f.read(), file_name=name,
                                   mime="text/plain", key=f"dl_{name}")
        st.info(r["nota"])
    st.markdown("""**Paquetes R requeridos:**
```r
install.packages(c("vegan", "ggplot2", "readr", "glmmTMB", "DHARMa", "lavaan", "semPlot"))
BiocManager::install("DESeq2")
```""")

with t_historial:
    section_title("Plantillas de la base de datos")
    from omniplast.reporting import plantillas as _plt
    cA, cB = st.columns(2, gap="large")
    with cA:
        st.download_button("PLANTILLA vacía (estructura + validaciones)",
                           _plt.generar_plantilla_bytes(),
                           file_name="PLANTILLA_OMNIPLAST.xlsx",
                           mime="application/vnd.openxmlformats-officedocument."
                                "spreadsheetml.sheet")
    with cB:
        st.download_button("EJEMPLO lleno (dataset demo sintético)",
                           _plt.generar_ejemplo_bytes(),
                           file_name="EJEMPLO_OMNIPLAST.xlsx",
                           mime="application/vnd.openxmlformats-officedocument."
                                "spreadsheetml.sheet")
    section_title("Resultados guardados")
    res = dao.list_results(limit=30)
    if not res.empty:
        tabla(res)
    else:
        st.caption("Guarde resultados desde las páginas de análisis o agentes con el botón "
                   "'Guardar en informes'.")
    section_title("Exportar todo el proyecto (CSV)")
    if st.button("Exportar tablas maestras (Excel)"):
        from ui.bootstrap import exportar_dfs
        exportar_dfs("omniplast_proyecto", {
            "sitios": dao.sites_df(), "muestras": dao.samples_df(),
            "microplasticos": dao.mp_df(), "ambiente": dao.env_wide(),
            "biota": dao.bio_df(), "qpcr": dao.qpcr_df(),
            "microbioma": dao.microbiome_wide()})

footer()
