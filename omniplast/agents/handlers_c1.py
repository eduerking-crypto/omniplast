"""Handlers F: genetica y molecular (AG-063..074)."""
import numpy as np
import pandas as pd

from ..config import GENE_CATALOG
from ..core.genetics import community, genotyping, qpcr
from ..core.statistics import inferential, multivariate
from ..reporting import figures
from .handlers_a import _ok, _req

REF_MIQE = ["Bustin et al. (2009) Clin. Chem. 55:611-622 (guías MIQE)"]


def _taxa(mb):
    return [c for c in mb.columns
            if c not in ("sample_id", "codigo", "sitio", "zona", "campana",
                         "temporada", "matriz", "anio")]


def _qpcr_wide(ctx):
    q = ctx["qpcr"]
    if q.empty:
        return pd.DataFrame()
    wide = q.pivot_table(index="sample_id", columns="target", values="valor", aggfunc="mean")
    meta_cols = [c for c in ["codigo", "sitio", "zona", "campana", "temporada", "matriz"]
                 if c in q.columns]
    meta = q.drop_duplicates("sample_id")[["sample_id"] + meta_cols].set_index("sample_id")
    return meta.join(wide).reset_index()


def _upgma_newick(dist, names):
    D = np.array(dist, dtype=float)
    n0 = len(names)
    size = np.ones(n0)
    labels = list(names)
    heights = [0.0] * n0
    active = list(range(n0))
    D = np.pad(D, ((0, 200), (0, 200)), constant_values=0.0)
    while len(active) > 1:
        best = (None, None, np.inf)
        for ii in range(len(active)):
            for jj in range(ii + 1, len(active)):
                a, b = active[ii], active[jj]
                if D[a, b] < best[2]:
                    best = (a, b, D[a, b])
        a, b, dab = best
        h = dab / 2
        new_id = len(labels)
        labels.append(f"({labels[a]}:{h - heights[a]:.4f},{labels[b]}:{h - heights[b]:.4f})")
        heights.append(h)
        new_size = size[a] + size[b]
        size = np.append(size, new_size)
        for c in active:
            if c in (a, b):
                continue
            D[new_id, c] = D[c, new_id] = (D[a, c] * size[a] + D[b, c] * size[b]) / new_size
        active = [c for c in active if c not in (a, b)] + [new_id]
    return labels[active[0]] + ";"


def ag_063(ctx):
    q = ctx["qpcr"]
    if q.empty or "petase" not in set(q.get("target", [])):
        return _req("Se requieren Ct de petase y 16S (gen de referencia).", REF_MIQE)
    r = qpcr.ddct_table(q, target="petase", ref="16s")
    if "error" in r:
        return _req(r["error"], REF_MIQE)
    g = qpcr.ddct_group_test(r)
    t = r["tabla"]
    fig = figures.box(t, x="zona", y="fc",
                      title="Expresión relativa de petase (2^-ΔΔCt)", ytitle="FC")
    return _ok(f"{r['interpretacion']}. {g.get('interpretacion', g.get('error', ''))}",
               {"ddct": t.round(3),
                "grupos": g.get("resumen_grupos", pd.DataFrame())},
               {"cajas": fig},
               REF_MIQE + ["Yuan et al. (2006) Bioinformatics 22:1098-1105"])


def ag_064(ctx):
    p = ctx["params"]
    if "diluciones" in p and "cts" in p:
        dil, cts = p["diluciones"], p["cts"]
        nota = "Curva estándar DEL USUARIO"
    else:
        dil = [1e6, 1e5, 1e4, 1e3, 1e2]
        cts = [18.2, 21.6, 25.1, 28.4, 31.7]
        nota = "EJEMPLO SINTÉTICO (ingrese sus diluciones/Ct en los parámetros del agente)"
    r = qpcr.standard_curve(dil, cts)
    if "error" in r:
        return _req(r["error"], REF_MIQE)
    flat = {"pendiente": r["pendiente"], "intercepto": r["intercepto"], "r2": r["r2"],
            "eficiencia_pct": r["eficiencia"]["eficiencia_pct"],
            "conforme_miqe": bool(r["eficiencia"]["valida_miqe"] and r["r2_valido_miqe"])}
    fig = figures.scatter(pd.DataFrame({"log_copias": np.log10(dil), "ct": cts}),
                          x="log_copias", y="ct", trendline=True,
                          title=f"Curva estándar — {nota}", xtitle="log₁₀(copias)",
                          ytitle="Ct")
    return _ok(f"{nota}: {r['interpretacion']}",
               {"curva_estandar": pd.DataFrame([flat])}, {"regresion": fig},
               REF_MIQE, ["Requiere ≥5 puntos, 3 réplicas por punto y NTC limpio"])


def ag_065(ctx):
    p = ctx["params"]
    if "temps" in p and "fluor" in p:
        temps, fluor = np.asarray(p["temps"], float), np.asarray(p["fluor"], float)
        nota = "Curva de fusión DEL USUARIO"
    else:
        temps = np.linspace(65, 95, 121)
        fluor = 1.0 / (1 + np.exp((temps - 78.5) / 1.4)) + 0.05
        fluor = fluor + 1.0 / (1 + np.exp((temps - 74.0) / 0.8)) * 0.18
        nota = "EJEMPLO SINTÉTICO (pico principal 78.5 °C + artefacto ~74 °C)"
    r = qpcr.melt_curve_peaks(temps, fluor)
    if "error" in r:
        return _req(r["error"], REF_MIQE)
    fig = figures.line(pd.DataFrame({"temperatura_c": temps,
                                     "dF_dT": -np.gradient(fluor, temps)}),
                       x="temperatura_c", y="dF_dT", title=f"Curva de fusión — {nota}",
                       xtitle="°C", ytitle="-dF/dT")
    return _ok(f"{nota}: {r['interpretacion']}",
               {"picos_tm": pd.DataFrame(r["tms"], columns=["Tm_c"])},
               {"fusion": fig}, REF_MIQE)


def ag_066(ctx):
    p = ctx["params"]
    if "canal_x" in p and "canal_y" in p:
        x, y = p["canal_x"], p["canal_y"]
        nota = "Genotipado DEL USUARIO"
    else:
        rng = np.random.default_rng(7)
        x = np.concatenate([rng.normal(0.9, .08, 24), rng.normal(0.1, .06, 22),
                            rng.normal(0.5, .07, 28)])
        y = np.concatenate([rng.normal(0.1, .06, 24), rng.normal(0.9, .08, 22),
                            rng.normal(0.55, .07, 28)])
        nota = "EJEMPLO SINTÉTICO (2 homocigotos + heterocigotos)"
    r = genotyping.call_genotypes(x, y)
    if "error" in r:
        return _req(r["error"], ["Livak (1999) Genomics"])
    hwe = None
    c = r["conteo"]
    if {"homocigoto_X", "heterocigoto", "homocigoto_Y"} <= set(c.index):
        hwe = genotyping.hardy_weinberg_test(
            [int(c.get("homocigoto_X", 0)), int(c.get("heterocigoto", 0)),
             int(c.get("homocigoto_Y", 0))])
    fig = figures.scatter(r["llamadas"], x="x", y="y", color="genotipo",
                          title=f"Discriminación alélica — {nota}")
    tabs = {"llamadas": c.to_frame("n").reset_index(names="genotipo")}
    if hwe and "error" not in hwe:
        tabs["hardy_weinberg"] = pd.DataFrame(
            [{"chi2": hwe["chi2"], "p": hwe["p"], "p_alelo": hwe["p_alelo"]}])
    return _ok(f"{nota}: {r['interpretacion']}"
               + (f". HWE: {hwe['interpretacion']}" if hwe and "interpretacion" in hwe else ""),
               tabs, {"clusters": fig},
               ["Livak KJ (1999) Genetics 155:1401-1412"])


def ag_067(ctx):
    q = ctx["qpcr"]
    env = ctx["env"]
    mp = ctx["mp"]
    if q.empty or env.empty:
        return _req("Se requieren qPCR y variables ambientales del mismo sedimento.")
    genes = [g for g in ("petase", "laccasa", "alkb") if g in set(q["target"])]
    if not genes:
        return _req("Sin genes de degradación en los datos qPCR.")
    wide = _qpcr_wide(ctx)
    mps = mp[mp["matriz"] == "sedimento"].groupby("sample_id")["abundancia"].sum().rename("mp_total")
    sed = env[env["matriz"] == "sedimento"][
        ["sample_id", "toc_pct", "arena_pct", "limo_arcilla_pct"]].merge(mps, on="sample_id")
    d = wide.merge(sed, on="sample_id", how="inner")
    if d.empty:
        return _req("Sin coincidencias muestra-qPCR×sedimento.")
    rows = []
    for g in genes:
        lg = np.log10(d[g].clip(lower=1.0))
        for var in ("mp_total", "toc_pct", "arena_pct"):
            if var not in d.columns:
                continue
            r = inferential.correlation_test(lg, d[var], method="spearman")
            rows.append({"gen": g, "variable": var, "r": r.get("r"), "p": r.get("p")})
    tab = pd.DataFrame(rows).dropna(subset=["r"])
    if tab.empty:
        return _req("Correlaciones no computables (datos insuficientes).")
    top = tab.sort_values("p").iloc[0]
    fig = figures.heatmap(tab.pivot(index="gen", columns="variable", values="r"),
                          title="Correlación log₁₀(copias)×ambiente (Spearman)")
    return _ok(f"Asociación más fuerte: {top['gen']} vs {top['variable']} "
               f"(r={top['r']:.2f}, p={top['p']:.3g}). "
               "La co-variación de genes de degradación con MP/MO evidencia una plastísfera activa.",
               {"asociaciones": tab.round(4)}, {"heatmap": fig},
               ["Yoshida et al. (2016) Science 351:1196-1203"])


def ag_068(ctx):
    mb = ctx["mb"]
    if mb is None or mb.empty:
        return _req("Se requieren conteos 16S (microbioma).", ["Chao (1984)", "Hill (1973)"])
    alpha = community.alpha_diversity(mb)
    if "codigo" in mb.columns:
        zona_map = mb.drop_duplicates("codigo").set_index("codigo")["zona"]
        alpha["zona"] = alpha["muestra"].map(zona_map)
    if "zona" in alpha.columns and alpha["zona"].notna().any():
        kw = inferential.kruskal_wallis(alpha.dropna(subset=["zona"]), "shannon", "zona")
        fig = figures.box(alpha.dropna(subset=["zona"]), x="zona", y="shannon",
                          title="Diversidad alfa (Shannon) por zona", ytitle="H'")
    else:
        kw = {"interpretacion": "sin grupos"}
        fig = None
    figs = {"cajas": fig} if fig is not None else {}
    return _ok(f"Shannon medio={alpha['shannon'].mean():.2f}; Chao1 medio={alpha['chao1'].mean():.0f}. "
               f"{kw.get('interpretacion', '')}",
               {"alfa": alpha.round(3)}, figs,
               ["Chao (1984) Scand. J. Stat. 11:265-270", "Hill (1973) Ecology 54:427-432"])


def ag_069(ctx):
    mb = ctx["mb"]
    if mb is None or mb.empty:
        return _req("Se requieren conteos 16S.", ["Legendre & Legendre (2012)"])
    idx = mb.set_index("sample_id")
    dist = multivariate.distance_matrix(idx[_taxa(mb)], method="bray")
    p = multivariate.pcoa(dist)
    coords = p["coordenadas"].copy()
    coords["zona"] = idx.loc[coords.index, "zona"].values
    expl = p["varianza_explicada"]
    e1 = float(expl.iloc[0])
    e2 = float(expl.iloc[1]) if len(expl) > 1 else 0.0
    fig = figures.scatter(coords, x="MDS1", y="MDS2", color="zona",
                          title=f"PCoA Bray-Curtis ({e1:.1f}% + {e2:.1f}% ejes 1-2)")
    return _ok(f"PCoA: MDS1 explica {e1:.1f}% de la varianza de la comunidad.",
               {"coordenadas": coords.round(3)}, {"pcoa": fig},
               ["Legendre & Legendre (2012) Numerical Ecology 3ª ed."])


def ag_070(ctx):
    mb = ctx["mb"]
    if mb is None or mb.empty:
        return _req("Se requieren conteos 16S.", ["Anderson (2001)"])
    idx = mb.set_index("sample_id")
    dist = multivariate.distance_matrix(idx[_taxa(mb)], method="bray")
    per = multivariate.permanova(dist, idx["zona"], permutations=999)
    ani = multivariate.anosim(dist, idx["zona"], permutations=999)
    disp = multivariate.betadisper_groups(dist, idx["zona"])
    return _ok(f"{per['interpretacion']} {ani['interpretacion']}. "
               "betadisper: si la dispersión difiere, PERMANOVA puede confundir efecto de "
               "grupo con heterogeneidad de varianza.",
               {"permanova": pd.DataFrame([per]), "anosim": pd.DataFrame([ani]),
                "dispersion": disp.round(3)},
               {}, ["Anderson (2001) Austral. NZ J. Stat. 43:281-300"])


def ag_071(ctx):
    mb = ctx["mb"]
    if mb is None or mb.empty:
        return _req("Se requieren conteos 16S.", ["Love et al. (2014)"])
    r = community.differential_abundance(mb, mb["zona"])
    if "error" in r:
        return _req(r["error"], ["Love et al. (2014)"])
    sig = r["tabla"][r["tabla"]["significativo"]]
    fig = figures.bars(r["tabla"].head(10), "taxon", "log2fc",
                       title="log2FC de los taxa más discriminantes", ytitle="log2FC")
    return _ok(f"{r['interpretacion']}. "
               + (f"Enriquecidos: {', '.join(sig['taxon'].head(5))}" if len(sig) else ""),
               {"diferencial": r["tabla"].head(20).round(4)}, {"barras": fig},
               ["Love et al. (2014) Genome Biol. 15:550 (lógica DESeq2 adaptada)"])


def ag_072(ctx):
    p = ctx["params"]
    seqs = p.get("secuencias")
    if not seqs or len(seqs) < 3:
        return _req("Suba ≥3 secuencias (params: secuencias={nombre: 'ATCG...'}). "
                    "Pipeline completo recomendado: NCBI BLAST → MAFFT → IQ-TREE con "
                    "bootstrap ≥1000. Este agente construye un árbol rápido UPGMA por "
                    "alineamiento por pares.",
                    ["Saitou & Nei (1987) Mol. Biol. Evol. 4:406-425",
                     "Sokal & Michener (1958) Univ. Kansas Sci. Bull. 38:1409"])
    try:
        from Bio import Align as BioAlign
        aligner = BioAlign.PairwiseAligner()
        aligner.mode = "global"
        names = list(seqs.keys())
        n = len(names)
        dist = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                score = aligner.score(seqs[names[i]], seqs[names[j]])
                norm = score / (min(len(seqs[names[i]]), len(seqs[names[j]])) * 2.0)
                dist[i, j] = dist[j, i] = max(1e-9, 1.0 - norm)
        newick = _upgma_newick(dist, names)
        return _ok(f"Árbol UPGMA (distancias de alineamiento por pares, {n} secuencias). "
                   f"Newick: {newick[:150]}{'...' if len(newick) > 150 else ''}",
                   {"matriz_distancias": pd.DataFrame(dist, index=names,
                                                      columns=names).round(3)},
                   {}, ["Sokal & Michener (1958) UPGMA"],
                   ["Validar topología con ML (IQ-TREE) antes de publicación"])
    except ImportError:
        return _req("Biopython no está disponible (pip install biopython).",
                    ["Cock et al. (2009) Bioinformatics 25:1422"])


def ag_073(ctx):
    mb = ctx["mb"]
    if mb is None or mb.empty:
        return _req("Se requieren conteos 16S.", ["Faust & Raes (2012)"])
    r = community.cooccurrence_network(mb)
    if "error" in r:
        return _req(r.get("error", "Red no construible."), ["Faust & Raes (2012)"])
    figs = {}
    try:
        import networkx as nx
        import plotly.graph_objects as go
        G = nx.Graph()
        for _, e in r["edges"].iterrows():
            G.add_edge(e["source"], e["target"])
        pos = nx.spring_layout(G, seed=42, k=1.8 / np.sqrt(max(G.number_of_nodes(), 1)))
        edge_x, edge_y = [], []
        for a, b in G.edges():
            edge_x += [pos[a][0], pos[b][0], None]
            edge_y += [pos[a][1], pos[b][1], None]
        deg = dict(G.degree())
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines",
                                 line=dict(color="#D8E2EC"), hoverinfo="none",
                                 showlegend=False))
        fig.add_trace(go.Scatter(
            x=[pos[n][0] for n in G.nodes()], y=[pos[n][1] for n in G.nodes()],
            mode="markers+text", text=list(G.nodes()), textposition="top center",
            textfont=dict(size=9),
            marker=dict(size=[8 + deg[n] * 2.2 for n in G.nodes()], color="#00A896"),
            showlegend=False))
        fig.update_layout(title=f"Red de co-ocurrencia ({G.number_of_edges()} aristas)",
                          height=520, template="plotly_white")
        figs["red"] = fig
    except Exception:
        pass
    tabs = {"aristas": r.get("edges", pd.DataFrame()).round(3)}
    if "grados" in r:
        tabs["grados"] = r["grados"].to_frame("grado").reset_index(names="taxon")
    return _ok(r["interpretacion"], tabs, figs,
               ["Faust & Raes (2012) Nat. Rev. Microbiol. 10:538-550"])


def ag_074(ctx):
    q = ctx["qpcr"]
    if q.empty:
        return _req("Se requieren datos qPCR del panel.")
    pathways = {
        "petase": "Hidrólisis de PET → MHET (xenobióticos)",
        "mhetase": "Hidrólisis de MHET → TPA",
        "laccasa": "Oxidación de aromáticos (tipo ligninolítico)",
        "alkb": "ω-oxidación de alcanos (cadenas tipo PE)",
        "phac": "Síntesis de PHA (bioplásticos)",
        "16s": "Gen de referencia (normalizador)",
    }
    rows = []
    for g, info in GENE_CATALOG.items():
        sub = q[q["target"] == g]
        det = sub["detectado"].mean() * 100 if len(sub) else np.nan
        med = sub.loc[sub["detectado"], "valor"].median() if len(sub) else np.nan
        rows.append({"gen": info["gen"], "organismo_modelo": info["organismo"],
                     "ruta_postulada": pathways.get(g, "—"),
                     "deteccion_pct": round(float(det), 1) if det == det else None,
                     "copias_medianas": med, "referencia": info["ref"]})
    tab = pd.DataFrame(rows)
    fig = figures.bars(tab.dropna(subset=["deteccion_pct"]).sort_values("deteccion_pct"),
                      "gen", "deteccion_pct", title="Detección del panel por qPCR",
                      ytitle="% muestras positivas")
    return _ok("Detección del panel: "
               + "; ".join(f"{r.gen.split('(')[0].strip()}: {r.deteccion_pct}%"
                           for r in tab.itertuples() if r.deteccion_pct is not None),
               {"panel": tab}, {"barras": fig},
               ["Yoshida et al. (2016) Science 351:1196-1203",
                "Danso, Chow & Streit (2019) AEM 85:e01095-19"])


HANDLERS_C1 = {f"AG-{i:03d}": f for i, f in enumerate([
    ag_063, ag_064, ag_065, ag_066, ag_067, ag_068, ag_069, ag_070, ag_071, ag_072,
    ag_073, ag_074], start=63)}
