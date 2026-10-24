import sys
from pathlib import Path

_here = Path(__file__).resolve()
for _p in _here.parents:
    if (_p / "app" / "ui" / "bootstrap.py").exists():
        for _s in (str(_p), str(_p / "app")):
            if _s not in sys.path:
                sys.path.insert(0, _s)
        break

import streamlit as st
from ui.bootstrap import footer, page_header, resumen_box, section_title

from omniplast.assistant import (KB, ask, ollama_available, recommend_test,
                                 sampling_design, search_kb)

page_header("Asistente PLASTI",
            "Asistente metodológico doctoral: responde con citas primarias y sugiere la "
            "ruta en OMNIPLAST (usa Ollama local si está disponible)")

estado = "Ollama ACTIVO (respuestas generativas con contexto KB)" if ollama_available() \
    else "Ollama NO disponible — modo base de conocimiento local (reglas + KB)"
resumen_box(f"<b>Motor:</b> {estado}", "ok" if ollama_available() else "warn")

t1, t2, t3 = st.tabs(["Chat con PLASTI", "Asesor estadístico", "Base de conocimiento"])

with t1:
    if "plasti_msgs" not in st.session_state:
        st.session_state["plasti_msgs"] = []
    for msg in st.session_state["plasti_msgs"]:
        with st.chat_message(msg["role"]):
            st.markdown(msg["content"])
    pregunta = st.chat_input("Pregúntale sobre protocolos, índices, estadística, genes…")
    if pregunta:
        st.session_state["plarti"] = None
        st.session_state["plasti_msgs"].append({"role": "user", "content": pregunta})
        with st.chat_message("user"):
            st.markdown(pregunta)
        with st.chat_message("assistant"):
            with st.spinner("PLASTI pensando…" if ollama_available() else "Consultando KB…"):
                r = ask(pregunta)
            st.markdown(r["respuesta"])
            st.caption(f"Fuente: {r.get('fuente', '?')} · modelo: {r.get('modelo', '?')}")
        st.session_state["plasti_msgs"].append(
            {"role": "assistant",
             "content": r["respuesta"] + f"\n\n_(fuente: {r.get('fuente')})_"})

with t2:
    st.markdown("El asesor decide la prueba según su diseño y supuestos "
                "(Zar 2010; Anderson 2001; Cameron & Trivedi 2013).")
    c1, c2, c3, c4 = st.columns(4)
    n_grupos = c1.selectbox("N° de grupos a comparar", [1, 2, 3, "4+"])
    pareado = c2.checkbox("Diseño pareado (mismo sitio/individuo)")
    conteos = c3.checkbox("Variable = conteo (abundancias)")
    multiv = c4.checkbox("Respuesta multivariada (comunidad)")
    c5, c6, c7 = st.columns(3)
    normalidad = None
    if not (conteos or multiv):
        normalidad = c5.selectbox("Normalidad (Shapiro p>0.05)", ["desconocida", "sí", "no"])
        normalidad = {"desconocida": None, "sí": True, "no": False}[normalidad]
    hetero = None
    if normalidad:
        hetero = c6.selectbox("Varianzas (Levene p<0.05)", ["homogéneas", "heterogéneas"])
        hetero = {"homogéneas": False, "heterogéneas": True}[hetero]
    repetido = c7.checkbox("Medidas repetidas (>2 por unidad)")
    if st.button("Recomendar prueba", type="primary"):
        rec = recommend_test(objetivo="comparar", n_grupos=int(n_grupos) if str(n_grupos).isdigit() else 4,
                             pareado=pareado, normalidad=normalidad,
                             heterocedasticidad=hetero, conteos=conteos,
                             multivariado=multiv, repetido=repetido)
        resumen_box(f"<b>Prueba recomendada:</b> {rec['prueba']}<br>"
                    f"<b>Justificación:</b> {rec['justificacion']}<br>"
                    f"<b>Alternativa:</b> {rec['alternativa']}<br>"
                    f"<b>Supuestos a verificar:</b> {rec['supuestos']}<br>"
                    f"<b>Referencia:</b> {rec['ref']}")
    section_title("Diseño muestral a priori")
    d = sampling_design()
    st.json(d, expanded=False)

with t3:
    busq = st.text_input("Buscar en la base de conocimiento…").strip().lower()
    if busq:
        hits = search_kb(busq, top=6)
    else:
        hits = [{"clave": k, **v} for k, v in KB.items()][:6]
    for h in hits:
        with st.expander(f"**{h['tema']}** — {h['resumen']}"):
            st.markdown(h["detalle"])
            st.caption("Fuentes: " + "; ".join(h["fuentes"]))
    st.caption(f"{len(KB)} entradas curadas · regla de honestidad: sin valores normativos "
               "inventados; todo lo normativo se marca 'verificar vigencia'.")

footer()
