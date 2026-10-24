"""Analisis multivariante comunitario: ordenaciones, distancias, PERMANOVA, SIMPER, RDA.

Convenciones tipo vegan (Oksanen et al.): transformaciones antes de distancias,
PERMANOVA por permutacion de etiquetas (Anderson 2001), RDA con pseudo-F por
permutacion de las restricciones (ter Braak & Smilauer 2012).
"""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.spatial import procrustes as _procrustes


def distance_matrix(df, method="bray", transform="sqrt"):
    """df: observaciones x variables (abundancias o quimica)."""
    data = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    if method in ("bray", "jaccard") and transform == "sqrt":
        data = np.sqrt(data)
    if method == "bray":
        d = pdist(data.values, metric="braycurtis")
    elif method == "jaccard":
        d = pdist((data > 0).astype(float).values, metric="jaccard")
    elif method == "euclidean":
        d = pdist(data.values, metric="euclidean")
    elif method == "aitchison":
        data = data.clip(lower=1e-9)
        clr = np.log(data).sub(np.log(data).mean(axis=1), axis=0)
        d = pdist(clr.values, metric="euclidean")
    elif method == "hellinger":
        hell = np.sqrt(data.div(data.sum(axis=1).replace(0, 1), axis=0))
        d = pdist(hell.values, metric="euclidean")
    else:
        d = pdist(data.values, metric="euclidean")
    return pd.DataFrame(squareform(d), index=df.index, columns=df.index)


def pca(df, n_components=None):
    data = df.select_dtypes(include=[np.number]).apply(pd.to_numeric, errors="coerce")
    data = data.dropna(axis=1, how="all").dropna(axis=0)
    cols = data.columns.tolist()
    z = (data - data.mean()) / (data.std(ddof=0).replace(0, 1))
    cov = np.cov(z.values.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    eigvals = np.clip(eigvals, 0, None)
    total = eigvals.sum()
    expl = eigvals / total * 100
    n = n_components or min(5, len(cols))
    n = min(n, len(cols))
    scores = z.values @ eigvecs[:, :n] * np.sqrt(eigvals[:n])
    loadings = pd.DataFrame(eigvecs[:, :n], index=cols,
                            columns=[f"PC{i+1}" for i in range(n)])
    kaiser = int((eigvals > 1).sum())
    return {"scores": pd.DataFrame(scores, index=data.index,
                                   columns=[f"PC{i+1}" for i in range(n)]),
            "loadings": loadings,
            "explicada": pd.Series(expl[:n], index=[f"PC{i+1}" for i in range(n)]),
            "explicada_total": float(expl[:2].sum()) if n >= 2 else float(expl[0]),
            "kaiser_componentes": kaiser,
            "variables": cols}


def varimax(loadings, gamma=1.0, max_iter=100, tol=1e-6):
    L = loadings.values.copy()
    p, k = L.shape
    R = np.eye(k)
    d_old = 0
    for _ in range(max_iter):
        LR = L @ R
        u, s, vt = np.linalg.svd(L.T @ (LR ** 3 - (gamma / p) * LR @ np.diag(np.diag(LR.T @ LR))))
        R = u @ vt
        d = s.sum()
        if d_old != 0 and abs(d - d_old) / d_old < tol:
            break
        d_old = d
    return pd.DataFrame(L @ R, index=loadings.index, columns=loadings.columns)


def factor_analysis(df, n_factors=2):
    from sklearn.decomposition import FactorAnalysis
    data = df.select_dtypes(include=[np.number]).apply(pd.to_numeric, errors="coerce").dropna()
    z = (data - data.mean()) / (data.std(ddof=0).replace(0, 1))
    fa = FactorAnalysis(n_components=n_factors, random_state=42).fit(z.values)
    load = pd.DataFrame(fa.components_.T, index=data.columns,
                        columns=[f"AF{i+1}" for i in range(n_factors)])
    communal = (load ** 2).sum(axis=1)
    load["comunalidad"] = communal
    return {"loadings": load, "n_factors": n_factors, "n_obs": len(data)}


def correspondence_analysis(table):
    """table: filas x columnas de conteos."""
    t = np.asarray(table, float) + 1e-12
    n = t.sum()
    P = t / n
    row_m = P.sum(axis=1)
    col_m = P.sum(axis=0)
    E = np.outer(row_m, col_m)
    S = (P - E) / np.sqrt(E)
    u, s, vt = np.linalg.svd(S, full_matrices=False)
    inertia = s ** 2
    row_scores = (u * s) / np.sqrt(row_m)[:, None]
    col_scores = (vt.T * s) / np.sqrt(col_m)[:, None]
    return {"inercia_total": float(inertia.sum()),
            "inercia": pd.Series(inertia[:min(5, len(s))],
                                 index=[f"Dim{i+1}" for i in range(min(5, len(s)))]),
            "score_filas": pd.DataFrame(row_scores[:, :min(5, row_scores.shape[1])],
                                       index=getattr(table, "index", range(t.shape[0]))),
            "score_columnas": pd.DataFrame(col_scores[:, :min(5, col_scores.shape[1])],
                                           columns=getattr(table, "columns", range(t.shape[1])))}


def hierarchical_clustering(df, method="ward", metric="euclidean"):
    data = df.select_dtypes(include=[np.number]).apply(pd.to_numeric, errors="coerce")
    data = data.dropna(axis=1, how="all")
    data = data.fillna(data.median(numeric_only=True))
    z = linkage(data.values, method=method, metric=metric if method != "ward" else "euclidean")
    cop = float(cophenet(z, pdist(data.values))[0])
    return {"linkage": z, "cophenetic": cop, "labels": list(data.index)}


def kmeans_clustering(df, k_min=2, k_max=8, features=None):
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    data = df[features] if features else df.select_dtypes(include=[np.number])
    data = data.apply(pd.to_numeric, errors="coerce").dropna()
    z = (data - data.mean()) / (data.std(ddof=0).replace(0, 1))
    rows = []
    models = {}
    for k in range(k_min, min(k_max, len(z) - 1) + 1):
        km = KMeans(n_clusters=k, n_init=25, random_state=42).fit(z.values)
        sil = silhouette_score(z.values, km.labels_) if len(np.unique(km.labels_)) > 1 else np.nan
        rows.append({"k": k, "inercia": float(km.inertia_), "silueta": float(sil)})
        models[k] = km
    return {"tabla": pd.DataFrame(rows), "models": models, "data": z, "columns": list(data.columns)}


def pcoa(dist):
    """dist: matriz de distancias (DataFrame o array cuadrado)."""
    D = np.asarray(dist, float)
    n = D.shape[0]
    A = -0.5 * D ** 2
    G = A - A.mean(axis=0, keepdims=True) - A.mean(axis=1, keepdims=True) + A.mean()
    eigvals, eigvecs = np.linalg.eigh(G)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    pos = eigvals > 1e-10
    coords = eigvecs[:, pos] * np.sqrt(eigvals[pos])
    perc = eigvals[pos] / eigvals[pos].sum() * 100
    labels = list(dist.index) if hasattr(dist, "index") else list(range(n))
    return {"coordenadas": pd.DataFrame(coords[:, :min(5, coords.shape[1])],
                                        index=labels,
                                        columns=[f"MDS{i+1}" for i in range(min(5, coords.shape[1]))]),
            "varianza_explicada": pd.Series(perc[:min(5, coords.shape[1])],
                                            index=[f"MDS{i+1}" for i in range(min(5, coords.shape[1]))])}


def nmds(dist, ndim=2, n_init=12, max_iter=400, tol=1e-7, seed=42):
    """SMACOF clasico (stress-1 de Kruskal)."""
    rng = np.random.default_rng(seed)
    D = np.asarray(dist, float)
    n = D.shape[0]
    starts = []
    p0 = pcoa(pd.DataFrame(D))
    starts.append(np.asarray(p0["coordenadas"].values[:, :ndim], float))
    best_stress, best_x = np.inf, None
    for i_init in range(n_init):
        if i_init < len(starts):
            x = starts[i_init].copy()
        else:
            x = rng.normal(size=(n, ndim))
            x = (x - x.mean(axis=0)) * np.sqrt((D ** 2).sum() / 2 / n) / np.sqrt((x ** 2).sum() / n + 1e-12)
        old_stress = np.inf
        for _it in range(max_iter):
            dist_hat = squareform(pdist(x))
            with np.errstate(divide="ignore", invalid="ignore"):
                b = -D / np.where(dist_hat == 0, 1e-12, dist_hat)
            np.fill_diagonal(b, 0)
            np.fill_diagonal(b, -b.sum(axis=1))
            x_new = b @ x / n
            x_new = x_new - x_new.mean(axis=0)
            dist_hat = squareform(pdist(x_new))
            stress = np.sqrt(((dist_hat - D) ** 2).sum() / (D ** 2).sum())
            x = x_new
            if abs(old_stress - stress) < tol:
                break
            old_stress = stress
        if stress < best_stress:
            best_stress, best_x = stress, x
    base = pcoa(pd.DataFrame(D))
    ref = base["coordenadas"].values[:, :ndim]
    if best_x is not None:
        try:
            _, m1, m2 = _procrustes(ref, best_x)
            best_x = best_x * m2 + m1
        except Exception:
            pass
    labels = list(dist.index) if hasattr(dist, "index") else list(range(n))
    quality = ("excelente" if best_stress < .05 else "buena" if best_stress < .10 else
               "aceptable" if best_stress < .20 else "pobre")
    return {"coordenadas": pd.DataFrame(best_x, index=labels,
                                         columns=[f"NMDS{i+1}" for i in range(ndim)]),
            "stress": float(best_stress), "calidad_ajuste": quality,
            "regla": "Kruskal: <0.05 excelente, <0.10 buena, <0.20 util con cautela"}


def _ss_permanova(D, labels):
    n = D.shape[0]
    total = D[np.triu_indices(n, 1)].sum() / n * 2 / 2
    ss_t = (D[np.triu_indices(n, 1)] ** 2).sum() / n
    ss_w = 0.0
    for g in np.unique(labels):
        idx = np.where(labels == g)[0]
        if len(idx) < 2:
            continue
        sub = D[np.ix_(idx, idx)]
        ss_w += (sub[np.triu_indices(len(idx), 1)] ** 2).sum() / len(idx)
    ss_a = ss_t - ss_w
    a = len(np.unique(labels))
    f = (ss_a / (a - 1)) / (ss_w / (n - a)) if ss_w > 0 else np.nan
    return f, ss_a / ss_t if ss_t > 0 else np.nan


def permanova(dist, groups, permutations=999, seed=42):
    """PERMANOVA de un factor (Anderson 2001). groups: serie alineada a dist."""
    rng = np.random.default_rng(seed)
    D = np.asarray(dist, float)
    labels = np.asarray(pd.Series(groups).astype(str).values)
    n = len(labels)
    a = len(np.unique(labels))
    if n < 4 or a < 2:
        return {"test": "PERMANOVA", "error": "n≥4 y ≥2 grupos"}
    f_obs, r2 = _ss_permanova(D, labels)
    count = 0
    for _ in range(permutations):
        perm = rng.permutation(labels)
        f_perm, _ = _ss_permanova(D, perm)
        if f_perm >= f_obs:
            count += 1
    p = (count + 1) / (permutations + 1)
    return {"test": "PERMANOVA (Anderson 2001)", "pseudo_F": float(f_obs), "p": float(p),
            "R2": float(r2), "permutaciones": permutations, "k_grupos": a, "n": n,
            "interpretacion": f"Pseudo-F={f_obs:.3f}, R²={r2:.3f}, p={p:.4f} "
                               f"(la estructura multivariada {'difiere' if p < .05 else 'no difiere'} entre grupos)"}


def anosim(dist, groups, permutations=999, seed=42):
    rng = np.random.default_rng(seed)
    D = np.asarray(dist, float)
    labels = np.asarray(pd.Series(groups).astype(str).values)
    n = len(labels)
    iu = np.triu_indices(n, 1)
    ranks = stats.rankdata(D[iu])
    labels_pair = labels[iu[0]]
    within = ranks[labels_pair == labels[iu[1]]]
    between = ranks[labels_pair != labels[iu[1]]]
    n_total = n * (n - 1) / 2
    r_obs = (between.mean() - within.mean()) / (n_total / 2)
    count = 0
    for _ in range(permutations):
        perm = rng.permutation(labels)
        lp = perm[iu[0]]
        w = ranks[lp == perm[iu[1]]]
        b = ranks[lp != perm[iu[1]]]
        if (b.mean() - w.mean()) / (n_total / 2) >= r_obs:
            count += 1
    p = (count + 1) / (permutations + 1)
    return {"test": "ANOSIM", "R": float(r_obs), "p": float(p), "permutaciones": permutations,
            "interpretacion": f"R={r_obs:.3f} ({'grupos bien separados' if r_obs > .5 else 'solapamiento parcial' if r_obs > .25 else 'grupos solapados'}), p={p:.4f}"}


def simper(df_counts, groups, transform="sqrt"):
    """Contribucion de cada variable a la disimilitud entre grupos (Clarke 1993)."""
    data = df_counts.apply(pd.to_numeric, errors="coerce").fillna(0)
    if transform == "sqrt":
        data = np.sqrt(data)
    lab = pd.Series(groups).astype(str).values
    levels = sorted(np.unique(lab))
    rows = []
    overall = {}
    for i in range(len(levels)):
        for j in range(i + 1, len(levels)):
            ia = np.where(lab == levels[i])[0]
            ib = np.where(lab == levels[j])[0]
            if len(ia) == 0 or len(ib) == 0:
                continue
            contrib = {}
            for col in data.columns:
                vals_a, vals_b = data[col].values[ia], data[col].values[ib]
                pairs = [(a, b) for a in vals_a for b in vals_b]
                num = np.mean([abs(a - b) for a, b in pairs])
                den = np.mean([(a + b) for a, b in pairs])
                contrib[col] = (num / den) if den > 0 else 0.0
            total_c = sum(contrib.values())
            for col, c in contrib.items():
                rows.append({"grupo_par": f"{levels[i]} vs {levels[j]}", "variable": col,
                             "contribucion_pct": c / total_c * 100 if total_c > 0 else 0.0})
            overall[f"{levels[i]} vs {levels[j]}"] = total_c / len(data.columns) * 100
    df = pd.DataFrame(rows)
    out = {}
    for par, sub in df.groupby("grupo_par"):
        sub = sub.sort_values("contribucion_pct", ascending=False)
        sub["acumulado_pct"] = sub["contribucion_pct"].cumsum()
        out[par] = sub
    return {"por_par": out, "disimilitud_media_pct": overall}


def rda(df_env, df_species, permutations=999, seed=42, scale_env=True):
    """RDA: especies ~ ambiente. Devuelve scores y pseudo-F por permutacion."""
    rng = np.random.default_rng(seed)
    Y = df_species.apply(pd.to_numeric, errors="coerce").fillna(0).values
    X = df_env.apply(pd.to_numeric, errors="coerce").fillna(0).values
    n, p = X.shape
    if n < 5 or p < 1:
        return {"test": "RDA", "error": "n≥5 observaciones"}
    if scale_env:
        sd = X.std(axis=0)
        X = (X - X.mean(axis=0)) / np.where(sd == 0, 1.0, sd)
    Xc = np.column_stack([np.ones(n), X])
    beta, *_ = np.linalg.lstsq(Xc, Y, rcond=None)
    Y_fit = Xc @ beta
    Y_res = Y - Y_fit
    def _pca(M):
        C = M.T @ M / len(M)
        vals, vecs = np.linalg.eigh(C)
        order = np.argsort(vals)[::-1]
        return vals[order], vecs[:, order]
    vals_c, vecs_c = _pca(Y_fit)
    vals_r, _ = _pca(Y_res)
    k = int((vals_c > 1e-10).sum())
    k = max(1, k)
    vals_c_pos = np.clip(vals_c[:k], 0, None)
    total_var = vals_c_pos.sum() + np.clip(vals_r, 0, None).sum()
    site_scores = Y_fit @ vecs_c[:, :k]
    sp_scores = vecs_c[:, :k] * np.sqrt(vals_c[:k])
    f_obs = (vals_c_pos.sum() / max(p, 1)) / (np.clip(vals_r, 0, None).sum() / max(n - p - 1, 1)) if vals_r.sum() != 0 else np.nan
    count = 0
    for _ in range(permutations):
        idx = rng.permutation(n)
        Xp = Xc[idx]
        beta_p, *_ = np.linalg.lstsq(Xp, Y, rcond=None)
        Yfp = Xp @ beta_p
        Yrp = Y - Yfp
        vp, _ = _pca(Yfp)
        rp, _ = _pca(Yrp)
        fp = (np.clip(vp, 0, None).sum() / max(p, 1)) / (np.clip(rp, 0, None).sum() / max(n - p - 1, 1))
        if fp >= f_obs:
            count += 1
    pval = (count + 1) / (permutations + 1)
    return {"test": "RDA (redundancy analysis)",
            "varianza_total": float(total_var),
            "varianza_restringida": float(vals_c_pos.sum()),
            "pct_restringida": float(vals_c_pos.sum() / total_var * 100) if total_var else np.nan,
            "pseudo_F": float(f_obs) if f_obs == f_obs else np.nan,
            "p": float(pval), "permutaciones": permutations,
            "eig_canonicas": pd.Series(vals_c[:min(5, k)] / total_var if total_var else vals_c[:k],
                                       index=[f"RDA{i+1}" for i in range(min(5, k))]),
            "site_scores": pd.DataFrame(site_scores[:, :min(5, site_scores.shape[1])],
                                         index=df_species.index,
                                         columns=[f"RDA{i+1}" for i in range(min(5, site_scores.shape[1]))]),
            "species_scores": pd.DataFrame(sp_scores[:, :min(5, sp_scores.shape[1])],
                                           index=df_species.columns,
                                           columns=[f"RDA{i+1}" for i in range(min(5, sp_scores.shape[1]))]),
            "interpretacion": f"{vals_c_pos.sum() / total_var * 100:.1f}% de la varianza de las "
                              f"especies explicada por el ambiente; pseudo-F={f_obs:.2f}, p={pval:.4f}"}


def betadisper_groups(dist, groups):
    """Distancia al centroide por grupo (para evaluar dispersion detras de PERMANOVA)."""
    from sklearn.metrics import pairwise_distances
    D = np.asarray(dist, float)
    lab = pd.Series(groups).astype(str).values
    rows = []
    for g in np.unique(lab):
        idx = np.where(lab == g)[0]
        if len(idx) < 2:
            continue
        sub = D[np.ix_(idx, idx)]
        centroid = sub.mean(axis=1)
        rows.append({"grupo": g, "n": len(idx),
                     "dist_media_al_centroide": float(sub.mean(axis=1).mean()),
                     "dist_sd": float(np.sqrt(((centroid - centroid.mean()) ** 2).mean()))})
    return pd.DataFrame(rows)
