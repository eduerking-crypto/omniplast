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
from ui.bootstrap import (badge, footer, page_header, render_result, resumen_box,
                          section_title, tabla)
from ui.components import agente_card

from omniplast.agents import CATEGORIES, REGISTRY, run_agent
from omniplast.db import dao

page_header("Los 100 agentes especialistas",
            "Cada agente ejecuta un flujo analítico real sobre el proyecto activo con "
            "interpretación, tablas, figuras y referencias")

col1, col2 = st.columns([3, 2])
busqueda = col1.text_input("Buscar agente (nombre, técnica, variable…)").strip()
cat_sel = col2.multiselect("Categorías", list(CATEGORIES.keys()),
                            default=list(CATEGORIES.keys()))

resumen_box("Los agentes en <span class='omp-badge warn'>Doctoral</span> implementan "
            "técnicas de nivel doctoral (PERMANOVA, GLMM, SSD, kriging, FDI…); los "
            "<span class='omp-badge info'>Básico/Intermedio</span> cubren el flujo "
            "cotidiano. Todos logran su ejecución en <b>results_store</b>.", "ok" if True else "ok")

def agentes_filtrados():
    out = []
    for aid, a in REGISTRY.items():
        if a["categoria"] not in cat_sel:
            continue
        if busqueda:
            q = busqueda.lower()
            hay = (q in a["nombre"].lower() or q in a["desc"].lower()
                   or q in CATEGORIES[a["categoria"]][0].lower() or q in aid.lower())
            if not hay:
                continue
        out.append(a)
    return out

agentes = agentes_filtrados()
st.caption(f"{len(agentes)} agentes visibles")

if "ag_resultado" not in st.session_state:
    st.session_state["ag_resultado"] = {}

por_categoria = {}
for a in agentes:
    por_categoria.setdefault(a["categoria"], []).append(a)

tabs = st.tabs([f"{cat} · {CATEGORIES[cat][0]}" for cat in por_categoria])
for tab, (cat, lista) in zip(tabs, por_categoria.items()):
    with tab:
        for a in lista:
            if agente_card(a, key_idx=f"{cat}_{a['id']}"):
                with st.spinner(f"Ejecutando {a['nombre']}…"):
                    r = run_agent(a["id"])
                st.session_state["ag_resultado"][a["id"]] = r
            if a["id"] in st.session_state["ag_resultado"]:
                r = st.session_state["ag_resultado"][a["id"]]
                st.markdown(f"**Resultado de {a['id']} — {a['nombre']}** "
                            f"({r.get('duracion_s', 0)} s)")
                render_result(r)
                if r.get("estado") == "ok":
                    if st.button("Guardar en informes", key=f"save_{a['id']}",
                                 type="secondary"):
                        dao.save_result("Agentes", f"{a['id']} · {a['nombre']}",
                                        r["resumen"],
                                        {"refs": r.get("refs", []),
                                         "pasos": r.get("pasos", [])})
                        st.success("Guardado. Consúltelo en la página Reportes.")

section_title("Bitácora de ejecuciones")
bit = dao.agent_runs_df(limit=25)
if not bit.empty and "agent_id" in bit.columns:
    tabla(bit[["agent_id", "agent_name", "status", "duracion_s", "creado"]])
else:
    st.caption("Aún no hay ejecuciones registradas en esta sesión de base.")

footer()
