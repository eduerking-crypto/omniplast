"""Estadistica descriptiva: resumen, normalidad, homocedasticidad, transformaciones."""
import numpy as np
import pandas as pd
from scipy import stats


def describe(series):
    x = pd.to_numeric(series, errors="coerce").dropna().values
    if len(x) == 0:
        return {}
    q1, q3 = np.percentile(x, [25, 75])
    out = {
        "n": int(len(x)), "media": float(np.mean(x)), "mediana": float(np.median(x)),
        "sd": float(np.std(x, ddof=1)) if len(x) > 1 else 0.0,
        "min": float(np.min(x)), "max": float(np.max(x)),
        "q1": float(q1), "q3": float(q3), "iqr": float(q3 - q1),
        "cv_pct": float(np.std(x, ddof=1) / np.mean(x) * 100) if np.mean(x) != 0 and len(x) > 1 else np.nan,
        "sesgo": float(stats.skew(x)) if len(x) > 2 else np.nan,
        "curtosis": float(stats.kurtosis(x)) if len(x) > 3 else np.nan,
    }
    out["se"] = out["sd"] / np.sqrt(out["n"]) if out["n"] > 1 else np.nan
    return out


def describe_by_group(df, value, group):
    rows = []
    for g, sub in df.groupby(group, observed=True):
        d = describe(sub[value])
        d[group] = g
        rows.append(d)
    cols = ["n", "media", "sd", "se", "mediana", "q1", "q3", "iqr", "min", "max", "cv_pct", "sesgo", "curtosis"]
    base = {c: np.nan for c in cols}
    base[group] = ""
    rows = [{**base, **r} for r in rows]
    return pd.DataFrame(rows).set_index(group).reset_index()


def normality_tests(series):
    x = pd.to_numeric(series, errors="coerce").dropna().values
    rows = []
    if len(x) >= 3:
        w, p = stats.shapiro(x)
        rows.append({"prueba": "Shapiro-Wilk", "estadístico": float(w), "p": float(p),
                     "recomendado_n": "3 ≤ n ≤ 2000"})
    if len(x) >= 8:
        z = (x - np.mean(x)) / (np.std(x, ddof=1) if np.std(x, ddof=1) > 0 else 1)
        ks, p = stats.kstest(z, "norm")
        rows.append({"prueba": "Kolmogorov-Smirnov (Lilliefors aprox.)", "estadístico": float(ks),
                     "p": float(p), "recomendado_n": "n > 50, con precaución"})
    if len(x) >= 5:
        ad = stats.anderson(x, dist="norm")
        rows.append({"prueba": "Anderson-Darling", "estadístico": float(ad.statistic),
                     "p": float(ad.statistic < ad.critical_values[2]),
                     "recomendado_n": f"A²={ad.statistic:.3f} vs valor crítico 5%={ad.critical_values[2]:.3f}"})
    for r in rows:
        r["interpretación"] = ("Distribución no normal (p<0.05)" if r["p"] < 0.05
                               else "Sin evidencia de no normalidad")
    return pd.DataFrame(rows)


def homogeneity_tests(df, value, group):
    groups = [pd.to_numeric(g, errors="coerce").dropna().values
              for _, g in df.groupby(group, observed=True)[value]]
    groups = [g for g in groups if len(g) >= 2]
    rows = []
    if len(groups) >= 2:
        w, p = stats.levene(*groups, center="median")
        rows.append({"prueba": "Levene (mediana, robusta)", "estadístico": float(w), "p": float(p)})
        try:
            b, pb = stats.bartlett(*groups)
            rows.append({"prueba": "Bartlett (sensible a no normalidad)", "estadístico": float(b), "p": float(pb)})
        except Exception:
            pass
    for r in rows:
        r["interpretación"] = ("Varianzas heterogéneas (p<0.05)" if r["p"] < 0.05
                               else "Sin evidencia de heterocedasticidad")
    return pd.DataFrame(rows)


def suggest_transformation(series):
    x = pd.to_numeric(series, errors="coerce").dropna().values
    if len(x) < 5:
        return "n insuficiente"
    results = {"original": stats.shapiro(x)[1] if 3 <= len(x) <= 5000 else np.nan}
    pos = x[x > 0]
    if len(pos) == len(x):
        results["log"] = stats.shapiro(np.log(pos))[1] if 3 <= len(pos) <= 5000 else np.nan
        results["sqrt"] = stats.shapiro(np.sqrt(pos))[1] if 3 <= len(pos) <= 5000 else np.nan
        results["1/x"] = stats.shapiro(1 / pos)[1] if 3 <= len(pos) <= 5000 else np.nan
    from scipy.stats import boxcox
    if len(pos) == len(x):
        try:
            _, lam = boxcox(x)
            results[f"boxcox (λ={lam:.2f})"] = stats.shapiro(boxcox(x)[0])[1]
        except Exception:
            pass
    best = max(results, key=lambda k: (results[k] if results[k] == results[k] else -1))
    return {"p_valores": {k: round(v, 4) for v, k in [(v, k) for k, v in results.items()]},
            "mejor_transformacion": best}


def iqr_outliers(df, value, group=None):
    df = df.copy()
    vals = pd.to_numeric(df[value], errors="coerce")

    def flag(s):
        q1, q3 = np.nanpercentile(s, [25, 75])
        iqr = q3 - q1
        return (s < q1 - 1.5 * iqr) | (s > q3 + 1.5 * iqr)

    if group and group in df.columns:
        df["outlier"] = vals.groupby(df[group], observed=True).transform(flag)
    else:
        df["outlier"] = flag(vals)
    return df
