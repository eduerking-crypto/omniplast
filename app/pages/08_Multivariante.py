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

from omniplast.core.statistics import multivariate
from omniplast.db import dao
from omniplast.reporting import figures

page_header("Estadística III — Multivariante comunitaria",
            "PCA · RDA · PCoA · NMDS · PERMANOVA + betadisper · ANOSIM · SIMPER · clusters")

env = dao.env_wide()
mb = dao.microbiome_wide()
mp = dao.mp_df()
if env.empty:
    st.stop()

t_pca, t_rdapcoa, t_nmds, t_perm, t_cluster = st.tabs(
    ["PCA ambiental", "RDA / PCoA", "NMDS", "PERMANOVA · ANOSIM · SIMPER", "Clustering"])

num_cols = [c for c in env.select_dtypes("number").columns
            if c not in ("sample_id", "lat", "lon", "anio", "profundidad_m")
            and env[c].notna().sum() >= 15]

with t_pca:
    sel = st.multiselect("Variables para el PCA", num_cols, default=num_cols[:8],
                         key="pca_vars")
    if len(sel) >= 3:
        sub = env.dropna(subset=sel)
        if len(sub) < 5:
            resumen_box(f"Solo {len(sub)} filas completas para {len(sel)} variables: "
                        "reduzca variables o impute (la versión actual no imputa).", "warn")
        else:
            r = multivariate.pca(sub[sel])
            scores = r["scores"].copy()
            scores["zona"] = sub["zona"].values
            chart(figures.scatter(scores, x="PC1", y="PC2", color="zona",
                                  title=f"PCA — {r['explicada_total']:.0f}% var (PC1+PC2)"))
            tabla(r["loadings"].round(3))
            resumen_box(f"Kaiser (λ>1): <b>{r['kaiser_componentes']}</b> componentes retenibles.")
            if st.checkbox("Aplicar rotación Varimax"):
                rot = multivariate.varimax(r["loadings"])
                tabla(rot.round(3))

with t_rdapcoa:
    if not mb.empty:
        taxa = [c for c in mb.columns if c not in
                ("sample_id", "codigo", "sitio", "zona", "campana", "temporada",
                 "matriz", "anio")]
        idx = mb.set_index("sample_id")
        c1, c2 = st.columns(2)
        target = c1.selectbox("Matriz a ordenar", ["microbioma 16S", "ambiente"])
        if target == "microbioma 16S":
            r = multivariate.pcoa(multivariate.distance_matrix(idx[taxa], "bray"))
            coords = r["coordenadas"].copy()
            coords["zona"] = idx.loc[coords.index, "zona"].values
            chart(figures.scatter(coords, x="MDS1", y="MDS2", color="zona",
                                  title="PCoA Bray-Curtis del microbioma"))
            tabla(r["varianza_explicada"].round(1))
        else:
            sub = env.dropna(subset=num_cols[:6]).set_index("sample_id")
            r = multivariate.pcoa(multivariate.distance_matrix(sub[num_cols[:6]],
                                                              "euclidean"))
            coords = r["coordenadas"].copy()
            coords["zona"] = sub.loc[coords.index, "zona"].values
            chart(figures.scatter(coords, x="MDS1", y="MDS2", color="zona",
                                  title="PCoA ambiental (euclídea estandarizada)"))
        if not mp.empty:
            mpm = mp[mp["matriz"] == "sedimento"].groupby("sample_id")["abundancia"].sum()
            sed_env = env[env["matriz"] == "sedimento"].merge(
                mpm.rename("mp_total"), on="sample_id")
            sed_env = sed_env.merge(
                dao.qpcr_df().pivot_table(index="sample_id", columns="target",
                                          values="valor", aggfunc="max")
                .add_prefix("gen_"), on="sample_id", how="left")
            env_vars = [c for c in ("toc_pct", "arena_pct", "limo_arcilla_pct")
                        if c in sed_env.columns]
            sp_vars = [c for c in sed_env.columns if c.startswith("gen_")][:4] or ["mp_total"]
            cc = st.multiselect("Restricciones ambientales (RDA)", env_vars,
                                default=env_vars[:2])
            if len(cc) >= 2:
                clean = sed_env.dropna(subset=cc + sp_vars)
                if len(clean) >= 8:
                    rda = multivariate.rda(clean[cc], clean[sp_vars], permutations=499)
                    if "error" not in rda:
                        resumen_box(rda["interpretacion"])
                        st_sc = rda["site_scores"].copy()
                        st_sc["zona"] = clean["zona"].values
                        chart(figures.scatter(st_sc, x="RDA1", y="RDA2", color="zona",
                                              title="RDA triplot (sitios)"))
                        tabla(rda["species_scores"].round(3))

with t_nmds:
    if mb.empty:
        resumen_box("El NMDS requiere la tabla de microbioma.", "warn")
    else:
        taxa = [c for c in mb.columns if c not in
                ("sample_id", "codigo", "sitio", "zona", "campana", "temporada",
                 "matriz", "anio")]
        idx = mb.set_index("sample_id")
        runs = st.slider("Reinicios aleatorios", 4, 30, 12)
        r = multivariate.nmds(multivariate.distance_matrix(idx[taxa], "bray"),
                              n_init=runs)
        coords = r["coordenadas"].copy()
        coords["zona"] = idx.loc[coords.index, "zona"].values
        chart(figures.scatter(coords, x="NMDS1", y="NMDS2", color="zona",
                              title=f"NMDS — stress={r['stress']:.3f} ({r['calidad_ajuste']})"))
        resumen_box(f"Reglas de Kruskal: <0.05 excelente · <0.10 buena · <0.20 útil con "
                    "cautela · >0.30 cuestionable. Ajuste actual: "
                    f"<b>{r['calidad_ajuste']}</b>.")

with t_perm:
    if mb.empty:
        resumen_box("Se requiere microbioma para PERMANOVA comunitario.", "warn")
    else:
        taxa = [c for c in mb.columns if c not in
                ("sample_id", "codigo", "sitio", "zona", "campana", "temporada",
                 "matriz", "anio")]
        idx = mb.set_index("sample_id")
        dist = multivariate.distance_matrix(idx[taxa], "bray")
        grupos = idx["zona"]
        per = multivariate.permanova(dist, grupos, permutations=999)
        resumen_box(f"<b>PERMANOVA:</b> {per.get('interpretacion', per.get('error', ''))}")
        ani = multivariate.anosim(dist, grupos, permutations=999)
        resumen_box(f"<b>ANOSIM:</b> {ani.get('interpretacion', '')}")
        tabla(multivariate.betadisper_groups(dist, grupos).round(3))
        sim = multivariate.simper(mb[taxa], mb["zona"])
        first = list(sim["por_par"].keys())[0]
        section_title(f"SIMPER — {first}")
        tabla(sim["por_par"][first].head(12))

with t_cluster:
    sub = env.dropna(subset=num_cols[:8])
    site_mean = sub.groupby("sitio")[num_cols[:8]].mean().dropna()
    if len(site_mean) >= 4:
        r = multivariate.hierarchical_clustering(site_mean)
        resumen_box(f"Cophenetic = {r['cophenetic']:.2f} (>0.7 = dendrograma fiel a las "
                    "distancias originales).")
        from scipy.cluster.hierarchy import dendrogram
        import plotly.graph_objects as go
        fig = go.Figure()
        den = dendrogram(r["linkage"], labels=list(r["labels"]), no_plot=True)
        icoord = den["icoord"]
        dcoord = den["dcoord"]
        for i, d in zip(icoord, dcoord):
            fig.add_trace(go.Scatter(x=i, y=d, mode="lines",
                                     line=dict(color="#0A2540"), showlegend=False))
        labels = den["ivl"]
        for x0, y0, lab in zip([sum(ic[1:3]) / 2 for ic in icoord],
                               [max(dc) for dc in dcoord], labels):
            fig.add_annotation(x=x0, y=y0 + 1, text=lab, showarrow=False)
        fig.update_layout(title="Clustering jerárquico de estaciones (Ward)",
                          height=430, template="plotly_white")
        st.plotly_chart(fig, use_container_width=True)

footer()
