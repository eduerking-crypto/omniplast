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
from ui.bootstrap import chart, footer, page_header, resumen_box, section_title, tabla

from omniplast.core.genetics import community
from omniplast.core.statistics import inferential, multivariate
from omniplast.db import dao
from omniplast.reporting import figures

page_header("Microbioma & plastísfera",
            "Diversidad alfa/beta · abundancia diferencial · LEfSe · core · redes")

mb = dao.microbiome_wide()
if mb.empty:
    resumen_box("Sin tabla de microbioma (16S) en el proyecto.", "warn")
    st.stop()

taxa = [c for c in mb.columns if c not in
        ("sample_id", "codigo", "sitio", "zona", "campana", "temporada", "matriz", "anio")]

t_alfa, t_beta, t_da, t_redes, t_core = st.tabs(
    ["Diversidad alfa", "Beta (PCoA)", "Diferencial", "Redes de co-ocurrencia", "Core"])

with t_alfa:
    alpha = community.alpha_diversity(mb)
    alpha["zona"] = mb.drop_duplicates("sample_id").set_index("sample_id")["zona"].reindex(
        alpha["muestra"]).values
    tabla(alpha)
    metrica = st.selectbox("Métrica", ["shannon", "chao1", "simpson", "riqueza", "pielou"])
    chart(figures.box(alpha, x="zona", y=metrica, title=f"{metrica} por zona",
                      ytitle=metrica))
    kw = inferential.kruskal_wallis(alpha, metrica, "zona")
    resumen_box(kw.get("interpretacion", kw.get("error", "")))

with t_beta:
    idx = mb.set_index("sample_id")
    dist = multivariate.distance_matrix(idx[taxa], "bray")
    p = multivariate.pcoa(dist)
    coords = p["coordenadas"].copy()
    coords["zona"] = idx.loc[coords.index, "zona"].values
    e1 = p["varianza_explicada"].iloc[0]
    e2 = p["varianza_explicada"].iloc[1] if len(p["varianza_explicada"]) > 1 else 0
    chart(figures.scatter(coords, x="MDS1", y="MDS2", color="zona",
                          title=f"PCoA Bray-Curtis ({e1:.0f}%+{e2:.0f}%)"))
    per = multivariate.permanova(dist, idx["zona"], permutations=999)
    ani = multivariate.anosim(dist, idx["zona"], permutations=999)
    resumen_box(f"<b>PERMANOVA:</b> {per.get('interpretacion', '')}<br>"
                f"<b>ANOSIM:</b> {ani.get('interpretacion', '')}")
    tabla(multivariate.betadisper_groups(dist, idx["zona"]).round(3))

with t_da:
    modo = st.radio("Método", ["Abundancia diferencial (DESeq2-lite)", "LEfSe-lite"],
                    horizontal=True)
    if modo.startswith("Abundancia"):
        r = community.differential_abundance(mb, mb["zona"])
        resumen_box(r.get("interpretacion", r.get("error", "")))
        if "tabla" in r:
            chart(figures.bars(r["tabla"].head(12), "taxon", "log2fc",
                              title="log2FC por taxón", ytitle="log2FC"))
            tabla(r["tabla"].head(20).round(4))
    else:
        r = community.lefse_effect(mb, mb["zona"])
        resumen_box(r.get("interpretacion", r.get("error", "")))
        if "tabla" in r:
            tabla(r["tabla"].head(20).round(4))
            chart(figures.bars(r["tabla"].head(10), "taxon", "efecto",
                              title="Efecto estandarizado por taxón", ytitle="d"))

with t_redes:
    r = community.cooccurrence_network(mb)
    if "error" in r:
        resumen_box(r.get("error", ""), "warn")
    else:
        resumen_box(r["interpretacion"])
        try:
            import networkx as nx
            import plotly.graph_objects as go
            G = nx.Graph()
            for _, e in r["edges"].iterrows():
                G.add_edge(e["source"], e["target"])
            pos = nx.spring_layout(G, seed=42, k=1.8 / np.sqrt(max(G.number_of_nodes(), 1)))
            ex, ey = [], []
            for a, b in G.edges():
                ex += [pos[a][0], pos[b][0], None]
                ey += [pos[a][1], pos[b][1], None]
            deg = dict(G.degree())
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=ex, y=ey, mode="lines",
                                     line=dict(color="#D8E2EC"), showlegend=False))
            fig.add_trace(go.Scatter(
                x=[pos[n][0] for n in G.nodes()], y=[pos[n][1] for n in G.nodes()],
                mode="markers+text", text=list(G.nodes()), textposition="top center",
                textfont=dict(size=9),
                marker=dict(size=[8 + deg[n] * 2.2 for n in G.nodes()], color="#00A896"),
                showlegend=False))
            fig.update_layout(title=f"Red ({G.number_of_edges()} aristas), modularidad="
                                    f"{r.get('modularidad', 'n/d')}",
                               height=520, template="plotly_white")
            st.plotly_chart(fig, use_container_width=True)
        except Exception:
            tabla(r["edges"].round(3))
        if "comunidades" in r and r["comunidades"]:
            section_title("Comunidades detectadas (módulos putativos de consortios)")
            for i, com in enumerate(r["comunidades"][:4]):
                st.markdown(f"**Módulo {i+1}:** {', '.join(com)}")

with t_core:
    prev = st.slider("Prevalencia mínima (fracción de muestras)", 0.3, 1.0, 0.7)
    r = community.core_microbiome(mb, prevalence=prev)
    resumen_box(r["interpretacion"])
    tabla(r["prevalencias"].to_frame("prevalencia").reset_index(names="taxon").round(3))

footer()
