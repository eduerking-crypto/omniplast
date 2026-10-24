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

from omniplast.core.connectors import connectors
from omniplast.db import dao
from omniplast.reporting import figures

page_header("Literatura & síntesis de evidencia",
            "Búsqueda en Crossref/PubMed, red de co-citas por palabras clave y biblioteca local")

t_busq, t_biblio, t_red = st.tabs(["Búsqueda viva", "Biblioteca local",
                                    "Mapa de evidencia (co-palabras)"])

with t_busq:
    c1, c2 = st.columns([3, 2])
    q = c1.text_input("Consulta", "microplastics mangrove OR estuary Caribbean")
    desde = c2.number_input("Desde el año", 2015, 2026, 2019)
    if st.button("Buscar en Crossref y PubMed", type="primary"):
        cr = connectors.search_crossref(q, rows=40, from_year=int(desde))
        pm = connectors.search_pubmed(q, rows=20)
        if cr.get("ok"):
            st.session_state["lit_cr"] = cr["df"]
            st.success(f"Crossref: {cr['n']} resultados (CC0, citas incluidas)")
        else:
            resumen_box(f"Crossref: {cr.get('error')}", "error")
        if pm.get("ok"):
            st.session_state["lit_pm"] = pm["df"]
            st.info(f"PubMed: {pm['n']} resultados")
        else:
            resumen_box(f"PubMed: {pm.get('error')}", "error")
    if "lit_cr" in st.session_state:
        df = st.session_state["lit_cr"]
        tabla(df[["titulo", "anio", "revista", "citas", "doi"]].head(40))
        if st.button("Guardar estos resultados en la biblioteca"):
            refs = [{"titulo": r.titulo, "autores": r.autores, "anio": r.anio,
                     "revista": r.revista, "doi": r.doi, "fuente": "Crossref",
                     "resumen": "", "etiquetas": q} for r in df.itertuples()]
            dao.save_literature(refs)
            st.success(f"{len(refs)} referencias guardadas.")

with t_biblio:
    lit = dao.literature_df()
    if lit.empty:
        resumen_box("Biblioteca vacía: guarde resultados desde la pestaña anterior o la "
                    "página Bases de datos públicas.", "warn")
    else:
        tabla(lit)
        by_year = lit.groupby("anio").size().reset_index(name="n")
        chart(figures.bars(by_year, "anio", "n", title="Literatura por año",
                           ytitle="artículos"))

with t_red:
    st.markdown("Red conceptual de los títulos recuperados: se extraen los términos más "
                "frecuentes (excluyendo stopwords) y se cuentan co-ocurrencias en títulos.")
    if "lit_cr" in st.session_state:
        df = st.session_state["lit_cr"]
        import re
        from collections import Counter
        stopwords = {"microplastics", "microplastic", "plastic", "study", "effect",
                     "effects", "using", "based", "new", "review", "the", "and", "for",
                     "with", "from", "marine", "water", "environmental", "environment",
                     "pollution", "a", "in", "of", "on", "to", "by", "an"}
        docs = []
        for t in df["titulo"].fillna(""):
            words = set(w for w in re.findall(r"[a-zíéóá]{4,}", t.lower())
                        if w not in stopwords)
            docs.append(words)
        freq = Counter(w for d in docs for w in d)
        top = [w for w, _ in freq.most_common(18)]
        if len(top) >= 6:
            import plotly.graph_objects as go
            edges = Counter()
            for d in docs:
                for i, w1 in enumerate(top):
                    for w2 in top[i + 1:]:
                        if w1 in d and w2 in d:
                            edges[(w1, w2)] += 1
            fig = go.Figure()
            angles = [2 * 3.14159 * i / len(top) for i in range(len(top))]
            xs = [round(__import__("math").cos(a), 3) for a in angles]
            ys = [round(__import__("math").sin(a), 3) for a in angles]
            for (w1, w2), n in edges.items():
                i1, i2 = top.index(w1), top.index(w2)
                fig.add_trace(go.Scatter(x=[xs[i1], xs[i2]], y=[ys[i1], ys[i2]],
                                         mode="lines", showlegend=False,
                                         line=dict(width=min(1 + n, 5), color="#00A896"),
                                         hoverinfo="text",
                                         hovertext=f"{w1} × {w2}: {n}"))
            fig.add_trace(go.Scatter(x=xs, y=ys, mode="markers+text", text=top,
                                      textposition="top center",
                                      marker=dict(size=14, color="#0A2540"),
                                      showlegend=False))
            fig.update_layout(title="Red de co-palabras en títulos (evidencia recuperada)",
                              height=560, template="plotly_white")
            st.plotly_chart(fig, use_container_width=True)
            tabla(pd.DataFrame(freq.most_common(15), columns=["termino", "frecuencia"]))
    else:
        resumen_box("Ejecute primero una búsqueda en la pestaña 'Búsqueda viva'.", "warn")

footer()
