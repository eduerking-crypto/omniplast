"""Ecologia de comunidades microbianas (16S): diversidad alfa/beta, DA, core, LEfSe."""
import numpy as np
import pandas as pd
from scipy import stats

from ...utils import benjamini_hochberg
from . import qpcr  # noqa: F401


def alpha_diversity(counts_wide, taxa_cols=None):
    """counts_wide: muestras x taxa (lecturas)."""
    data = counts_wide.copy()
    cols = taxa_cols if taxa_cols is not None else [c for c in data.columns
                                                    if c not in ("sample_id", "codigo", "sitio",
                                                                 "zona", "campana", "temporada",
                                                                 "matriz", "anio")]
    X = data[cols].apply(pd.to_numeric, errors="coerce").fillna(0).values
    rows = []
    for i, (_, r) in enumerate(data.iterrows()):
        counts = X[i]
        present = counts[counts > 0]
        n_reads = counts.sum()
        richness = len(present)
        if len(present) > 0 and n_reads > 0:
            p = present / present.sum()
            shannon = float(-(p * np.log(p)).sum())
            simpson = float(1 - (p ** 2).sum())
            pielou = float(shannon / np.log(len(present))) if len(present) > 1 else np.nan
        else:
            shannon = simpson = pielou = np.nan
        s_obs = richness
        f1 = float((present == 1).sum())
        f2 = float((present == 2).sum())
        chao1 = s_obs + (f1 ** 2 / (2 * f2) if f2 > 0 else f1 * (f1 - 1) / 2)
        rows.append({"muestra": data.iloc[i].get("codigo", i), "lecturas": float(n_reads),
                    "riqueza": richness, "chao1": float(chao1),
                    "shannon": shannon, "simpson": simpson, "pielou": pielou})
    return pd.DataFrame(rows)


def rarefaction_curve(counts, max_depth=None, steps=25, seed=42):
    """Curva de rarefaccion individual por muestreo sin reemplazo (Hurlbert 1971)."""
    rng = np.random.default_rng(seed)
    counts = np.asarray(counts, float)
    counts = counts[counts > 0]
    if len(counts) == 0:
        return pd.DataFrame()
    total = int(counts.sum())
    pool = np.repeat(np.arange(len(counts)), counts.astype(int))
    max_depth = max_depth or total
    depths = np.unique(np.linspace(1, max_depth, steps).astype(int))
    rows = []
    for d in depths:
        draws = [len(np.unique(rng.choice(pool, size=min(d, total), replace=False)))
                for _ in range(20)]
        rows.append({"profundidad": int(d), "riqueza_media": float(np.mean(draws)),
                     "sd": float(np.std(draws))})
    return pd.DataFrame(rows)


def differential_abundance(counts_wide, groups, taxa_cols=None, fdr=0.05):
    """DESeq2-lite: size factors median-of-ratios + prueba por taxon + BH FDR."""
    data = counts_wide.copy()
    cols = taxa_cols if taxa_cols is not None else [c for c in data.columns
                                                    if c not in ("sample_id", "codigo", "sitio",
                                                                 "zona", "campana", "temporada",
                                                                 "matriz", "anio")]
    X = data[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    X = X.clip(lower=0)
    ref = np.exp(np.log(X.clip(lower=1e-9)).median(axis=0))
    geo = np.exp(np.log(X.clip(lower=1e-9)).mean(axis=1))
    size_factors = (X.clip(lower=1e-9) / ref).median(axis=1).replace(0, 1)
    norm = X.div(size_factors, axis=0)
    lab = pd.Series(groups).astype(str).values
    levels = np.unique(lab)
    if len(levels) < 2:
        return {"test": "Abundancia diferencial", "error": "≥2 grupos requeridos"}
    rows = []
    for taxon in cols:
        vals = norm[taxon].values
        a = vals[lab == levels[0]]
        b = vals[lab == levels[1]]
        a_clean = a[a > 0]
        b_clean = b[b > 0]
        if len(a_clean) < 3 or len(b_clean) < 3:
            rows.append({"taxon": taxon, "log2fc": np.nan, "p": np.nan,
                         "p_adj": np.nan, "media_1": float(np.mean(a)) if len(a) else np.nan,
                         "media_2": float(np.mean(b)) if len(b) else np.nan})
            continue
        try:
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        except ValueError:
            p = np.nan
        m1 = float(np.mean(np.log2(a_clean)))
        m2 = float(np.mean(np.log2(b_clean)))
        rows.append({"taxon": taxon, "log2fc": m2 - m1, "p": float(p),
                     "media_1": float(np.mean(a)), "media_2": float(np.mean(b))})
    tab = pd.DataFrame(rows)
    tab["p_adj"] = benjamini_hochberg(tab["p"].fillna(1).values)
    tab["significativo"] = (tab["p_adj"] < fdr) & (tab["log2fc"].abs() > 1)
    tab = tab.sort_values("p_adj")
    n_sig = int(tab["significativo"].sum())
    return {"test": "Abundancia diferencial (DESeq2-lite: size factors + MWU + BH)",
            "tabla": tab, "n_significativos": n_sig,
            "interpretacion": f"{n_sig} taxa con abundancia diferencial (FDR<{fdr}, |log2FC|>1)"}


def lefse_effect(counts_wide, groups, taxa_cols=None):
    """Analogo LEfSe: KW por taxon + efecto como diferencia estandarizada."""
    data = counts_wide.copy()
    cols = taxa_cols if taxa_cols is not None else [c for c in data.columns
                                                    if c not in ("sample_id", "codigo", "sitio",
                                                                 "zona", "campana", "temporada",
                                                                 "matriz", "anio")]
    X = data[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    rel = X.div(X.sum(axis=1).replace(0, 1), axis=0)
    lab = pd.Series(groups).astype(str).values
    levels = np.unique(lab)
    rows = []
    for taxon in cols:
        vals = rel[taxon].values
        groups_v = [vals[lab == lv] for lv in levels]
        groups_v = [g for g in groups_v if len(g) >= 3]
        if len(groups_v) < 2:
            continue
        try:
            h, p = stats.kruskal(*groups_v)
        except ValueError:
            continue
        pooled_sd = np.sqrt(np.mean([np.var(g, ddof=1) for g in groups_v if len(g) > 1])) or 1e-9
        effect = (np.mean(groups_v[0]) - np.mean(groups_v[1])) / pooled_sd if len(groups_v) == 2 else np.nan
        rows.append({"taxon": taxon, "efecto": float(effect), "kw_p": float(p),
                     "enriquecido_en": levels[0] if effect and effect > 0 else
                     (levels[1] if len(levels) > 1 and effect and effect < 0 else "—")})
    tab = pd.DataFrame(rows)
    tab["kw_p_adj"] = benjamini_hochberg(tab["kw_p"].fillna(1).values)
    tab = tab.sort_values("kw_p_adj")
    return {"test": "LEfSe-lite (KW + tamaño de efecto)",
            "tabla": tab,
            "interpretacion": f"{int((tab['kw_p_adj'] < .05).sum())} taxa discriminantes entre grupos"}


def core_microbiome(counts_wide, taxa_cols=None, prevalence=0.7, rel_threshold=0.001):
    data = counts_wide.copy()
    cols = taxa_cols if taxa_cols is not None else [c for c in data.columns
                                                    if c not in ("sample_id", "codigo", "sitio",
                                                                 "zona", "campana", "temporada",
                                                                 "matriz", "anio")]
    X = data[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    rel = X.div(X.sum(axis=1).replace(0, 1), axis=0)
    prev = (rel > rel_threshold).mean()
    core = prev[prev >= prevalence].sort_values(ascending=False)
    return {"test": f"Core microbioma (prevalencia ≥{prevalence:.0%}, rel≥{rel_threshold})",
            "core": core, "n_core": int(len(core)),
            "prevalencias": prev.sort_values(ascending=False),
            "interpretacion": f"{len(core)} taxa forman el núcleo: {', '.join(core.index[:8])}"
                              + ("..." if len(core) > 8 else "")}


def cooccurrence_network(counts_wide, taxa_cols=None, method="spearman", p_cut=0.05,
                         fdr=True, min_prevalence=0.4, top_n=25):
    """Red de co-ocurrencia con matriz de correlaciones esparsa (RMT-lite)."""
    data = counts_wide.copy()
    cols = taxa_cols if taxa_cols is not None else [c for c in data.columns
                                                    if c not in ("sample_id", "codigo", "sitio",
                                                                 "zona", "campana", "temporada",
                                                                 "matriz", "anio")]
    X = data[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    rel = X.div(X.sum(axis=1).replace(0, 1), axis=0)
    prev = (rel > 0).mean()
    keep = prev[prev >= min_prevalence].index.tolist()[:top_n]
    if len(keep) < 4:
        return {"test": "Red de co-ocurrencia", "error": "≥4 taxones prevalentes"}
    sub = rel[keep]
    n = len(sub)
    edges = []
    for i in range(len(keep)):
        for j in range(i + 1, len(keep)):
            a, b = sub[keep[i]].values, sub[keep[j]].values
            if method == "spearman":
                r, p = stats.spearmanr(a, b)
            else:
                r, p = stats.pearsonr(a, b)
            edges.append({"source": keep[i], "target": keep[j], "r": float(r), "p": float(p)})
    tab = pd.DataFrame(edges)
    if fdr:
        tab["p_adj"] = benjamini_hochberg(tab["p"].fillna(1).values)
        tab = tab[tab["p_adj"] < p_cut]
    else:
        tab = tab[tab["p"] < p_cut]
    if len(tab) < 2:
        return {"test": "Red de co-ocurrencia", "error": "Pocas asociaciones significativas",
                "tabla": tab}
    try:
        import networkx as nx
        G = nx.Graph()
        for _, r in tab.iterrows():
            G.add_edge(r["source"], r["target"], weight=abs(r["r"]), sign="+" if r["r"] > 0 else "-")
        deg = pd.Series(dict(G.degree())).sort_values(ascending=False)
        try:
            communities = nx.community.greedy_modularity_communities(G)
            mod = nx.community.modularity(G, communities)
        except Exception:
            mod, communities = np.nan, []
        top_hub = deg.index[0] if len(deg) else "—"
        return {"test": f"Red de co-ocurrencia ({method}, FDR={fdr})",
                "edges": tab, "n_nodos": int(len(keep)), "n_aristas": int(len(tab)),
                "grados": deg, "modularidad": float(mod) if mod == mod else np.nan,
                "comunidades": [list(c) for c in list(communities)[:8]],
                "hub_principal": str(top_hub),
                "interpretacion": f"{len(tab)} aristas; hub={top_hub}; modularidad={mod if mod != mod else round(mod, 3)}"}
    except ImportError:
        return {"test": "Red de co-ocurrencia", "tabla": tab,
                "interpretacion": "networkx no disponible; se devuelven las aristas"}
