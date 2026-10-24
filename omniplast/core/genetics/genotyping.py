"""Genotipado: llamadas alélicas, HWE, Fst, asociación gen-ambiente-fenotipo."""
import numpy as np
import pandas as pd
from scipy import stats

from ...utils import benjamini_hochberg


def call_genotypes(x_channel, y_channel, n_clusters=3, seed=42):
    """Discriminacion alelica por agrupamiento 2D de canales de fluorescencia."""
    from sklearn.cluster import KMeans
    data = pd.DataFrame({"x": pd.to_numeric(x_channel, errors="coerce"),
                         "y": pd.to_numeric(y_channel, errors="coerce")}).dropna()
    if len(data) < n_clusters * 3:
        return {"test": "Genotipado", "error": "Puntos insuficientes para clustering"}
    km = KMeans(n_clusters=n_clusters, n_init=25, random_state=seed).fit(data.values)
    labels = km.labels_
    centers = km.cluster_centers_
    no_template = centers[:, 0].argmin() if n_clusters >= 3 else None
    def _call(i):
        if i == no_template:
            return "NTC"
        return "homocigoto_X" if centers[i, 0] > centers[i, 1] else "homocigoto_Y"
    calls = []
    for lab in labels:
        calls.append(_call(lab))
    het_cluster = None
    if n_clusters >= 3:
        angles = np.degrees(np.arctan2(centers[:, 1], centers[:, 0]))
        het_cluster = int(np.argmin(np.abs(angles - 45)))
    for i, lab in enumerate(labels):
        if het_cluster is not None and lab == het_cluster and calls[i] != "NTC":
            calls[i] = "heterocigoto"
    counts = pd.Series(calls).value_counts()
    call_rate = float((pd.Series(calls) != "sin llamada").mean())
    return {"test": "Discriminación alélica (clustering 2D)",
            "llamadas": pd.DataFrame({"x": data["x"].values, "y": data["y"].values,
                                      "genotipo": calls}),
            "conteo": counts, "centros": centers.tolist(),
            "call_rate": call_rate,
            "interpretacion": f"Distribución: {counts.to_dict()}; call rate={call_rate:.1%}"}


def hardy_weinberg_test(genotype_counts):
    """genotype_counts: [n_hom1, n_het, n_hom2]. Chi² con frecuencia esperada p²:2pq:q²."""
    n_hom1, n_het, n_hom2 = genotype_counts
    n = n_hom1 + n_het + n_hom2
    if n < 5:
        return {"test": "Hardy-Weinberg", "error": "n≥5"}
    p = (2 * n_hom1 + n_het) / (2 * n)
    q = 1 - p
    exp = np.array([p ** 2 * n, 2 * p * q * n, q ** 2 * n])
    obs = np.array([n_hom1, n_het, n_hom2])
    mask = exp > 0
    chi2 = float(((obs[mask] - exp[mask]) ** 2 / exp[mask]).sum())
    dof = int(mask.sum() - 1)
    pval = float(stats.chi2.sf(chi2, dof)) if dof > 0 else np.nan
    return {"test": "Equilibrio Hardy-Weinberg", "chi2": chi2, "p": pval,
            "p_alelo": float(p), "q_alelo": float(q),
            "esperados": exp.tolist(), "observados": obs.tolist(),
            "interpretacion": ("Desviación del HWE (p<0.05): posible selección, "
                               "estructura poblacional o error de genotipado"
                               if pval < .05 else "Población en equilibrio Hardy-Weinberg")}


def fst_pairwise(pop1_genotypes, pop2_genotypes):
    """Fst simplificado de Weir & Cockerham (1984) para alelos bialélicos codificados 0/1/2."""
    g1 = np.asarray(pd.Series(pop1_genotypes).astype(float).dropna(), float) / 2
    g2 = np.asarray(pd.Series(pop2_genotypes).astype(float).dropna(), float) / 2
    if len(g1) < 3 or len(g2) < 3:
        return {"test": "Fst", "error": "n≥3 por población"}
    p1, p2 = g1.mean(), g2.mean()
    p_bar = (p1 * len(g1) + p2 * len(g2)) / (len(g1) + len(g2))
    n_bar = (len(g1) + len(g2)) / 2
    s2 = ((len(g1) * (p1 - p_bar) ** 2 + len(g2) * (p2 - p_bar) ** 2) /
          ((2 - 1) * n_bar))
    h1 = p1 * (1 - p1) * 2 * len(g1) / (len(g1) - 1) if len(g1) > 1 else 0
    h2 = p2 * (1 - p2) * 2 * len(g2) / (len(g2) - 1) if len(g2) > 1 else 0
    h_bar = (h1 * len(g1) + h2 * len(g2)) / (len(g1) + len(g2))
    num = s2 - h_bar / (2 * n_bar) / 1
    den = s2 + h_bar * (1 - 1 / (2 * n_bar)) if s2 + h_bar > 0 else np.nan
    fst = num / den if den and den > 0 else 0.0
    return {"test": "Fst (Weir & Cockerham 1984, simplificado)", "fst": float(fst),
            "interpretacion": f"Fst={fst:.4f} ({'diferenciación considerable' if fst > .15 else 'diferenciación moderada' if fst > .05 else 'estructura débil o nula'})"}


def snp_association(genotypes, phenotype, model="additive", logistic=False):
    """genotypes: 0/1/2 (alelos); phenotype: continuo o binario."""
    g = pd.to_numeric(pd.Series(genotypes), errors="coerce")
    y = pd.to_numeric(pd.Series(phenotype), errors="coerce")
    ok = g.notna() & y.notna()
    g, y = g[ok], y[ok]
    if model == "dominant":
        g = (g > 0).astype(float)
    elif model == "recessive":
        g = (g > 1).astype(float)
    if len(g) < 10 or g.nunique() < 2:
        return {"test": "Asociación SNP", "error": "n≥10 con variabilidad genética"}
    X = np.column_stack([np.ones(len(g)), g.values])
    if logistic or y.nunique() <= 2:
        try:
            import statsmodels.api as sm
            res = sm.Logit((y > 0).astype(float).values, X).fit(disp=0)
            beta = float(res.params[1])
            se = float(res.bse[1])
            pval = float(res.pvalues[1])
            or_ = float(np.exp(beta))
            ci = (float(np.exp(beta - 1.96 * se)), float(np.exp(beta + 1.96 * se)))
            return {"test": "Asociación SNP (regresión logística)", "odds_ratio": or_,
                    "ic95": ci, "p": pval, "modelo": model,
                    "interpretacion": f"OR={or_:.2f} (IC95% {ci[0]:.2f}-{ci[1]:.2f}), p={pval:.4g}"}
        except Exception:
            pass
    res = stats.linregress(g.values, y.values)
    return {"test": "Asociación SNP (regresión lineal)", "beta": float(res.slope),
            "r": float(res.rvalue), "p": float(res.pvalue), "modelo": model,
            "interpretacion": f"β={res.slope:.3g} por alelo, r={res.rvalue:.3f}, p={res.pvalue:.4g}"}


def panel_association(panel_df, phenotype, correction="BH"):
    """panel_df: muestras x SNPs (0/1/2); corrección multiple por panel."""
    results = []
    for snp in panel_df.columns:
        if panel_df[snp].nunique() < 2:
            continue
        res = snp_association(panel_df[snp], pd.Series(phenotype).values, model="additive")
        if "error" in res:
            continue
        res["snp"] = snp
        results.append(res)
    tab = pd.DataFrame(results)
    if tab.empty:
        return {"test": "Asociación de panel", "error": "Sin SNPs informativos"}
    tab["p_adj"] = benjamini_hochberg(tab["p"].fillna(1).values)
    tab["-log10_p"] = -np.log10(tab["p"].clip(lower=1e-12))
    tab = tab.sort_values("p")
    return {"test": "Asociación de panel genético (GWAS-lite)",
            "tabla": tab[["snp", "beta", "r", "p", "p_adj", "-log10_p"]],
            "interpretacion": f"{int((tab['p_adj'] < .05).sum())} marcadores significativos tras corrección {correction}"}
