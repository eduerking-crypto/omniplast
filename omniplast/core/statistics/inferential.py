"""Pruebas inferenciales parametricas y no parametricas con tamanos de efecto y post-hoc."""
import numpy as np
import pandas as pd
from scipy import stats

from ...utils import holm


def cohen_d(a, b, paired=False):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if paired and len(a) == len(b):
        d = a - b
        return float(np.mean(d) / (np.std(d, ddof=1) or 1))
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return np.nan
    sp = np.sqrt(((n1 - 1) * np.var(a, ddof=1) + (n2 - 1) * np.var(b, ddof=1)) / (n1 + n2 - 2))
    return float((np.mean(a) - np.mean(b)) / sp) if sp > 0 else np.nan


def cliffs_delta(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) == 0 or len(b) == 0:
        return np.nan
    gt = sum((a_i > b).sum() for a_i in a)
    lt = sum((a_i < b).sum() for a_i in a)
    return float((gt - lt) / (len(a) * len(b)))


def _interpret_d(d):
    d = abs(d) if d == d else 0
    if d < 0.2:
        return "trivial"
    if d < 0.5:
        return "pequeño"
    if d < 0.8:
        return "mediano"
    return "grande"


def ttest_two(df, value, group, g1=None, g2=None, paired=False, welch=True):
    gcol = df[group].astype(str)
    levels = sorted(gcol.unique())
    g1 = g1 or levels[0]
    g2 = g2 or levels[-1]
    sub = df[df[group].astype(str).isin([str(g1), str(g2)])].copy()
    sub["val"] = pd.to_numeric(sub[value], errors="coerce")
    sub = sub.dropna(subset=["val"])
    a = sub[sub[group].astype(str) == str(g1)]["val"].values
    b = sub[sub[group].astype(str) == str(g2)]["val"].values
    if len(a) < 2 or len(b) < 2:
        return {"test": "t de Student/Welch", "error": "Al menos 2 observaciones por grupo"}
    if paired and len(a) == len(b):
        t, p = stats.ttest_rel(a, b)
        d = cohen_d(a, b, paired=True)
    else:
        t, p = (stats.ttest_ind(a, b, equal_var=False) if welch
                else stats.ttest_ind(a, b, equal_var=True))
        d = cohen_d(a, b)
    return {"test": "t de Welch" if welch else "t de Student", "grupo_1": g1, "grupo_2": g2,
            "estadistico": float(t), "p": float(p), "df": None,
            "media_1": float(np.mean(a)), "media_2": float(np.mean(b)),
            "d_cohen": d, "magnitud": _interpret_d(d),
            "interpretacion": f"Media {g1}={np.mean(a):.3g} vs {g2}={np.mean(b):.3g}; "
                              f"d de Cohen={d:.2f} ({_interpret_d(d)})"}


def anova_one(df, value, group):
    sub = df[[group, value]].copy()
    sub[value] = pd.to_numeric(sub[value], errors="coerce")
    sub = sub.dropna()
    groups = [g[value].values for _, g in sub.groupby(group, observed=True)]
    groups = [g for g in groups if len(g) >= 2]
    if len(groups) < 2:
        return {"test": "ANOVA", "error": "Se requieren ≥2 grupos con n≥2"}
    f, p = stats.f_oneway(*groups)
    grand = sub[value].mean()
    ss_b = sum(len(g) * (g.mean() - grand) ** 2 for g in groups)
    ss_t = ((sub[value] - grand) ** 2).sum()
    eta2 = ss_b / ss_t if ss_t > 0 else np.nan
    k, n = len(groups), len(sub)
    omega2 = (ss_b - (k - 1) * ((ss_t - ss_b) / (n - k))) / (ss_t + (k - 1) * ((ss_t - ss_b) / (n - k))) if n > k else np.nan
    return {"test": "ANOVA de un factor", "estadistico": float(f), "p": float(p),
            "k_grupos": k, "n": n, "eta2": float(eta2), "omega2": float(omega2),
            "interpretacion": f"F({k-1},{n-k})={f:.3f}, p={p:.4g}; η²={eta2:.3f} "
                              f"({'alta' if eta2 > .14 else 'moderada' if eta2 > .06 else 'baja'} varianza explicada)"}


def kruskal_wallis(df, value, group):
    sub = df[[group, value]].copy()
    sub[value] = pd.to_numeric(sub[value], errors="coerce")
    sub = sub.dropna()
    groups = [g[value].values for _, g in sub.groupby(group, observed=True)]
    groups = [g for g in groups if len(g) >= 2]
    if len(groups) < 2:
        return {"test": "Kruskal-Wallis", "error": "Se requieren ≥2 grupos con n≥2"}
    h, p = stats.kruskal(*groups)
    k, n = len(groups), len(sub)
    epsilon2 = (h - k + 1) / (n - k) if n > k else np.nan
    return {"test": "Kruskal-Wallis", "estadistico": float(h), "p": float(p),
            "k_grupos": k, "n": n, "epsilon2": float(epsilon2),
            "interpretacion": f"H={h:.3f}, p={p:.4g}; ε²={epsilon2:.3f}"}


def dunn_posthoc(df, value, group):
    sub = df[[group, value]].copy()
    sub[value] = pd.to_numeric(sub[value], errors="coerce")
    sub = sub.dropna()
    sub[group] = sub[group].astype(str)
    if sub[group].nunique() < 3:
        return pd.DataFrame()
    ranks = stats.rankdata(sub[value].values)
    sub["rank"] = ranks
    n = len(sub)
    _, ties_counts = np.unique(sub[value].values, return_counts=True)
    tie_term = np.sum(ties_counts ** 3 - ties_counts) if len(ties_counts) else 0
    s2 = n * (n + 1) / 12 - tie_term / (12 * (n - 1)) if n > 1 else np.nan
    summary = sub.groupby(group)["rank"].agg(["count", "mean"])
    levels = list(summary.index)
    rows = []
    for i in range(len(levels)):
        for j in range(i + 1, len(levels)):
            l1, l2 = levels[i], levels[j]
            n1, n2 = summary.loc[l1, "count"], summary.loc[l2, "count"]
            se = np.sqrt(s2 * (1 / n1 + 1 / n2)) if s2 > 0 else np.nan
            z = (summary.loc[l1, "mean"] - summary.loc[l2, "mean"]) / se
            rows.append({"comparacion": f"{l1} vs {l2}", "z": float(z),
                         "p_sin_ajustar": float(2 * stats.norm.sf(abs(z)))})
    rows = pd.DataFrame(rows)
    rows["p_holm"] = holm(rows["p_sin_ajustar"].values)
    return rows


def tukey_hsd(df, value, group):
    sub = df[[group, value]].copy()
    sub[value] = pd.to_numeric(sub[value], errors="coerce")
    sub = sub.dropna()
    if sub[group].nunique() < 2:
        return pd.DataFrame()
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    res = pairwise_tukeyhsd(sub[value].values, sub[group].astype(str).values, alpha=0.05)
    return pd.DataFrame(res.summary().data[1:], columns=res.summary().data[0])


def games_howell(df, value, group):
    sub = df[[group, value]].copy()
    sub[value] = pd.to_numeric(sub[value], errors="coerce")
    sub = sub.dropna()
    sub[group] = sub[group].astype(str)
    stats_g = sub.groupby(group)[value].agg(["count", "mean", "var"])
    levels = list(stats_g.index)
    k = len(levels)
    rows = []
    for i in range(k):
        for j in range(i + 1, k):
            n1, n2 = stats_g.loc[levels[i], "count"], stats_g.loc[levels[j], "count"]
            v1, v2 = stats_g.loc[levels[i], "var"], stats_g.loc[levels[j], "var"]
            if n1 < 2 or n2 < 2:
                continue
            se = np.sqrt(v1 / n1 + v2 / n2)
            dfree = (v1 / n1 + v2 / n2) ** 2 / ((v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1))
            t = (stats_g.loc[levels[i], "mean"] - stats_g.loc[levels[j], "mean"]) / se
            try:
                p = stats.studentized_range.sf(abs(t) / np.sqrt(2), k, dfree)
            except Exception:
                p = 2 * stats.t.sf(abs(t), dfree)
            rows.append({"comparacion": f"{levels[i]} vs {levels[j]}",
                         "diferencia": float(stats_g.loc[levels[i], "mean"] - stats_g.loc[levels[j], "mean"]),
                         "q": float(abs(t) / np.sqrt(2)), "df": float(dfree), "p": float(p)})
    return pd.DataFrame(rows)


def pairwise_u(df, value, group):
    sub = df[[group, value]].copy()
    sub[value] = pd.to_numeric(sub[value], errors="coerce")
    sub = sub.dropna()
    sub[group] = sub[group].astype(str)
    levels = sorted(sub[group].unique())
    rows = []
    for i in range(len(levels)):
        for j in range(i + 1, len(levels)):
            a = sub.loc[sub[group] == levels[i], value].values
            b = sub.loc[sub[group] == levels[j], value].values
            if len(a) < 2 or len(b) < 2:
                continue
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            rows.append({"comparacion": f"{levels[i]} vs {levels[j]}", "U": float(u),
                         "cliffs_delta": cliffs_delta(a, b), "p_sin_ajustar": float(p)})
    rows = pd.DataFrame(rows)
    if not rows.empty:
        rows["p_holm"] = holm(rows["p_sin_ajustar"].values)
    return rows


def wilcoxon_paired(df, value, group, g1, g2):
    sub = df[df[group].astype(str).isin([str(g1), str(g2)])]
    a = sub[sub[group].astype(str) == str(g1)][value].astype(float).values
    b = sub[sub[group].astype(str) == str(g2)][value].astype(float).values
    n = min(len(a), len(b))
    if n < 3:
        return {"test": "Wilcoxon", "error": "Se requieren ≥3 pares"}
    w, p = stats.wilcoxon(a[:n], b[:n])
    return {"test": "Wilcoxon (muestras relacionadas)", "estadistico": float(w), "p": float(p), "n_pares": n}


def friedman_test(df, value, group, block):
    sub = df[[value, group, block]].copy()
    sub[value] = pd.to_numeric(sub[value], errors="coerce")
    wide = sub.pivot_table(index=block, columns=group, values=value, aggfunc="mean")
    wide = wide.dropna()
    if wide.shape[0] < 3 or wide.shape[1] < 3:
        return {"test": "Friedman", "error": "Se requieren ≥3 bloques y ≥3 condiciones"}
    chi2, p = stats.friedmanchisquare(*[wide[c].values for c in wide.columns])
    kendall_w = chi2 / (wide.shape[0] * (wide.shape[1] - 1)) if wide.shape[1] > 1 else np.nan
    return {"test": "Friedman", "estadistico": float(chi2), "p": float(p),
            "n_bloques": wide.shape[0], "k_condiciones": wide.shape[1],
            "W_kendall": float(kendall_w)}


def chi_square_test(table):
    table = np.asarray(table, float)
    chi2, p, dof, exp = stats.chi2_contingency(table)
    n = table.sum()
    cramers_v = float(np.sqrt(chi2 / (n * (min(table.shape) - 1)))) if min(table.shape) > 1 else np.nan
    return {"test": "Chi-cuadrado de independencia", "estadistico": float(chi2), "p": float(p),
            "df": int(dof), "v_cramer": cramers_v,
            "interpretacion": f"χ²({dof})={chi2:.3f}, p={p:.4g}; V de Cramer={cramers_v:.3f}"}


def correlation_test(x, y, method="pearson"):
    x = pd.to_numeric(pd.Series(x), errors="coerce")
    y = pd.to_numeric(pd.Series(y), errors="coerce")
    ok = x.notna() & y.notna()
    x, y = x[ok].values, y[ok].values
    if len(x) < 3:
        return {"test": "correlación", "error": "n < 3"}
    if method == "pearson":
        r, p = stats.pearsonr(x, y)
        fisher_z = np.arctanh(np.clip(r, -0.9999, 0.9999))
        se = 1 / np.sqrt(len(x) - 3)
        ci = (np.tanh(fisher_z - 1.96 * se), np.tanh(fisher_z + 1.96 * se))
    elif method == "spearman":
        r, p = stats.spearmanr(x, y)
        ci = (np.nan, np.nan)
    else:
        r, p = stats.kendalltau(x, y)
        ci = (np.nan, np.nan)
    strength = ("muy fuerte" if abs(r) > .8 else "fuerte" if abs(r) > .6 else
                "moderada" if abs(r) > .4 else "débil" if abs(r) > .2 else "trivial")
    return {"test": f"Correlación {method}", "r": float(r), "p": float(p), "n": int(len(x)),
            "ic95": (float(ci[0]), float(ci[1])),
            "interpretacion": f"r={r:.3f} ({strength}, {'positiva' if r > 0 else 'negativa'}), p={p:.4g}"}


def correlation_matrix(df, columns=None, method="pearson"):
    cols = columns or df.select_dtypes(include=[np.number]).columns.tolist()
    data = df[cols].apply(pd.to_numeric, errors="coerce")
    r = data.corr(method=method)
    p = pd.DataFrame(np.ones((len(cols), len(cols))), index=cols, columns=cols)
    for i in cols:
        for j in cols:
            if i == j:
                continue
            res = correlation_test(data[i], data[j], method="pearson" if method == "pearson" else method)
            p.loc[i, j] = res.get("p", np.nan)
    return r, p


def partial_correlation(df, x, y, covariates):
    data = df[[x, y] + list(covariates)].apply(pd.to_numeric, errors="coerce").dropna()
    if len(data) < len(covariates) + 3:
        return {"test": "correlación parcial", "error": "n insuficiente"}
    Xc = sm_add_constant(data[covariates])
    rx = sm_residuals(Xc, data[x])
    ry = sm_residuals(Xc, data[y])
    r, p = stats.pearsonr(rx, ry)
    return {"test": "Correlación parcial", "x": x, "y": y, "controlada_por": covariates,
            "r": float(r), "p": float(p), "n": int(len(data)),
            "interpretacion": f"r parcial={r:.3f}, p={p:.4g} (controlando {', '.join(covariates)})"}


def sm_add_constant(frame):
    import statsmodels.api as sm
    return sm.add_constant(frame)


def sm_residuals(Xc, target):
    import statsmodels.api as sm
    model = sm.OLS(target.values, Xc.values).fit()
    return model.resid
