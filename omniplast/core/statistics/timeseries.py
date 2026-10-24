"""Series temporales: Mann-Kendall, pendiente de Sen, STL, ACF, Holt."""
import numpy as np
import pandas as pd
from scipy import stats


def _mk_core(x):
    """S, varianza con correccion por empates y Z de Mann-Kendall."""
    x = np.asarray(x, float)
    n = len(x)
    s = 0
    for i in range(n - 1):
        s += np.sign(x[i + 1:] - x[i]).sum()
    _, counts = np.unique(x, return_counts=True)
    tie = np.sum(counts * (counts - 1) * (2 * counts + 5))
    var = (n * (n - 1) * (2 * n + 5) - tie) / 18
    if s > 0:
        z = (s - 1) / np.sqrt(var) if var > 0 else 0
    elif s < 0:
        z = (s + 1) / np.sqrt(var) if var > 0 else 0
    else:
        z = 0
    return s, var, z


def mann_kendall(series):
    """Mann-Kendall para tendencia (MODIFIED: Wu et al. 2008? no; Mann 1945, Kendall 1975)."""
    x = pd.to_numeric(pd.Series(series), errors="coerce").dropna().values
    if len(x) < 6:
        return {"test": "Mann-Kendall", "error": "Recomendado n ≥ 6"}
    s, var, z = _mk_core(x)
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    trend = "creciente" if z > 0 else "decreciente" if z < 0 else "sin tendencia"
    if abs(z) > 2.576:
        sig = "altamente significativa (α=0.01)"
    elif abs(z) > 1.96:
        sig = "significativa (α=0.05)"
    elif abs(z) > 1.645:
        sig = "marginal (α=0.10)"
    else:
        sig = "no significativa"
    return {"test": "Mann-Kendall", "S": float(s), "Z": float(z), "p": float(p),
            "tendencia": trend, "significancia": sig,
            "interpretacion": f"Z={z:.3f}, p={p:.4f} → tendencia {trend}, {sig}"}


def sens_slope(series, time=None):
    """Pendiente de Sen (1968): mediana de pendientes por pares."""
    s = pd.to_numeric(pd.Series(series), errors="coerce").dropna()
    x = s.values
    n = len(x)
    if n < 3:
        return {"test": "Pendiente de Sen", "error": "n < 3"}
    slopes = []
    for i in range(n - 1):
        slopes.extend((x[i + 1:] - x[i]) / (np.arange(i + 1, n) - i))
    med = float(np.median(slopes))
    nts = np.array(slopes)
    se = 1.253 * nts.std(ddof=1) / np.sqrt(n) if n > 2 else np.nan
    ci = (float(np.percentile(nts, 2.5)), float(np.percentile(nts, 97.5))) if n > 5 else (np.nan, np.nan)
    return {"test": "Pendiente de Sen", "pendiente": med, "ic95": ci,
            "interpretacion": f"Cambio mediano por paso temporal = {med:.4g} "
                              f"({(f'IC95%: {ci[0]:.3g} a {ci[1]:.3g}' if ci[0] == ci[0] else '')})"}


def seasonal_mann_kendall(df, date_col, value_col, season_col=None, freq="M"):
    """Mann-Kendall estacional (Hirsch et al. 1982)."""
    d = df.copy()
    d[date_col] = pd.to_datetime(d[date_col], errors="coerce")
    d[value_col] = pd.to_numeric(d[value_col], errors="coerce")
    d = d.dropna(subset=[date_col, value_col]).sort_values(date_col)
    if season_col is None:
        season_col = "estacion"
        d[season_col] = d[date_col].dt.month
    results = []
    s_tot, var_tot = 0.0, 0.0
    for season, sub in d.groupby(season_col, observed=True):
        if len(sub) < 3:
            continue
        s_i, var_i, z_i = _mk_core(sub[value_col].values)
        s_tot += s_i
        var_tot += var_i
        results.append({"estacion": season, "n": len(sub), "S": float(s_i), "Z": float(z_i),
                        "p": float(2 * (1 - stats.norm.cdf(abs(z_i))))})
    if not results or var_tot <= 0:
        return {"test": "MK estacional", "error": "Datos insuficientes por estación"}
    tab = pd.DataFrame(results)
    z_tot = s_tot / np.sqrt(var_tot) if s_tot != 0 else 0
    p_tot = 2 * (1 - stats.norm.cdf(abs(z_tot)))
    return {"test": "Mann-Kendall estacional (Hirsch 1982)", "por_estacion": tab,
            "k_estaciones": len(tab), "S_total": float(s_tot), "Z_total": float(z_tot),
            "p_total": float(p_tot),
            "interpretacion": f"Combinado estacional: Z={z_tot:.3f}, p={p_tot:.4f} → "
                              f"tendencia {'creciente' if z_tot > 0 else 'decreciente' if z_tot < 0 else 'nula'}"}


def stl_decompose(df, date_col, value_col, period=12):
    d = df.copy()
    d[date_col] = pd.to_datetime(d[date_col], errors="coerce")
    d[value_col] = pd.to_numeric(d[value_col], errors="coerce")
    d = d.dropna(subset=[date_col, value_col]).sort_values(date_col)
    s = pd.Series(d[value_col].values, index=pd.DatetimeIndex(d[date_col]))
    s = s.groupby(level=0).mean()
    freq = pd.infer_freq(s.index) or "MS"
    s = s.asfreq(freq)
    s = s.interpolate()
    if len(s.dropna()) < 2 * period:
        return {"test": "STL", "error": f"Se requieren ≥ {2 * period} observaciones"}
    from statsmodels.tsa.seasonal import STL
    res = STL(s, period=period, robust=True).fit()
    out = pd.DataFrame({"fecha": s.index, "observado": res.observed,
                        "tendencia": res.trend, "estacional": res.seasonal, "residuo": res.resid}).reset_index(drop=True)
    return {"test": "STL (descomposición robusta)", "componentes": out,
            "fuerza_estacional": float(max(0, 1 - np.var(res.resid) /
                                          max(np.var(res.seasonal + res.resid), 1e-12)))}


def acf_pacf(series, nlags=20):
    from statsmodels.tsa.stattools import acf, pacf
    x = pd.to_numeric(pd.Series(series), errors="coerce").dropna().values
    if len(x) < 4 * nlags:
        nlags = max(2, len(x) // 4)
    a = acf(x, nlags=nlags, fft=True)
    try:
        p = pacf(x, nlags=min(nlags, len(x) // 2 - 1))
    except Exception:
        p = np.array([np.nan])
    return {"acf": pd.DataFrame({"lag": range(len(a)), "acf": a}),
            "pacf": pd.DataFrame({"lag": range(len(p)), "pacf": p}),
            "ic_bartlett": 1.96 / np.sqrt(len(x))}


def holt_forecast(series, periods=6):
    from statsmodels.tsa.holtwinters import ExponentialSmoothing
    x = pd.to_numeric(pd.Series(series), errors="coerce").dropna().values
    if len(x) < 8:
        return {"test": "Holt", "error": "n < 8"}
    try:
        fit = ExponentialSmoothing(x, trend="add", damped_trend=True,
                                   initialization_method="estimated").fit()
        pred = fit.forecast(periods)
        return {"test": "Holt amortiguado", "pronostico": [float(v) for v in pred],
                "fitted": [float(v) for v in fit.fittedvalues],
                "aic": float(fit.aic)}
    except Exception as e:
        return {"test": "Holt", "error": str(e)}
