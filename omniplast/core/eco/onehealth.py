"""One Health: exposicion humana a MP via consumo de pesca y agua (modelos probabilísticos).

No existe umbral regulatorio internacional para exposicion humana a MP (WHO 2019:
evidencia insuficiente); los resultados se reportan como estimaciones con bandas de
incertidumbre y se comparan unicamente contra valores de contexto publicados.
"""
import numpy as np
import pandas as pd
from scipy import stats


def exposure_seafood(conc_p_per_g, meal_g=(80, 150, 300), meals_per_week=(1, 3, 7),
                    body_weight_kg=(55, 65, 80), n=50000, seed=42):
    """Monte Carlo de ingesta de MP por consumo de productos pesqueros.

    conc_p_per_g: serie de concentraciones medidas (particulas/g tejido) o (mu, sigma)
    para asumir lognormal.
    """
    rng = np.random.default_rng(seed)
    if isinstance(conc_p_per_g, tuple):
        mu, sigma = conc_p_per_g
        conc = rng.lognormal(mu, sigma, n)
    else:
        data = pd.to_numeric(pd.Series(conc_p_per_g), errors="coerce").dropna().values
        if len(data) < 5:
            return {"test": "Exposición", "error": "≥5 mediciones o parámetros lognormal (mu, sigma)"}
        ln = np.log(np.clip(data, 1e-6, None))
        mu, sigma = ln.mean(), ln.std()
        conc = rng.lognormal(mu, sigma, n)
    meal = rng.triangular(meal_g[0], meal_g[1], meal_g[2], n)
    freq = rng.triangular(meals_per_week[0], meals_per_week[1], meals_per_week[2], n) / 7
    bw = rng.triangular(body_weight_kg[0], body_weight_kg[1], body_weight_kg[2], n)
    p_day = conc * meal * freq
    p_kg_day = p_day / bw
    def _pct(x):
        return {p: float(np.percentile(x, p)) for p in (50, 75, 90, 95, 99)}
    inputs = {"concentracion": conc, "porcion_g": meal, "frecuencia_dia": freq, "peso_corporal": bw}
    sens = {}
    for k, v in inputs.items():
        try:
            rho, p = stats.spearmanr(v, p_kg_day)
            sens[k] = {"rho": float(rho), "p": float(p)}
        except Exception:
            pass
    sens_df = (pd.DataFrame(sens).T.sort_values("rho", key=abs, ascending=False)
               if sens else pd.DataFrame())
    mediana = float(np.median(p_kg_day))
    return {"test": "Exposición humana por consumo de pesca (Monte Carlo, One Health)",
            "p_kg_bw_dia": _pct(p_kg_day), "p_dia": _pct(p_day),
            "mediana_p_kg_bw_dia": mediana, "sensibilidad": sens_df, "n_iteraciones": n,
            "interpretacion": f"Mediana={mediana:.3g} partículas/kg peso corporal/día "
                              f"(P95={np.percentile(p_kg_day, 95):.3g}); contexto: no hay "
                              "valor toxicológico de referencia regulatorio para MP (WHO 2019); "
                              "usar para comparación entre escenarios, no como límite de seguridad"}


def exposure_drinking_water(conc_p_per_l, liters_per_day=(1.0, 2.0, 3.0), n=50000, seed=42):
    rng = np.random.default_rng(seed)
    if isinstance(conc_p_per_l, tuple):
        mu, sigma = conc_p_per_l
        conc = rng.lognormal(mu, sigma, n)
    else:
        data = pd.to_numeric(pd.Series(conc_p_per_l), errors="coerce").dropna().values
        if len(data) < 5:
            return {"test": "Exposición", "error": "≥5 mediciones o parámetros (mu, sigma)"}
        ln = np.log(np.clip(data, 1e-6, None))
        conc = rng.lognormal(ln.mean(), ln.std(), n)
    liters = rng.triangular(liters_per_day[0], liters_per_day[1], liters_per_day[2], n)
    p_day = conc * liters
    return {"test": "Exposición por agua de consumo",
            "p_dia": {p: float(np.percentile(p_day, p)) for p in (50, 75, 90, 95, 99)},
            "mediana_p_dia": float(np.median(p_day)), "n_iteraciones": n,
            "interpretacion": f"Mediana={np.median(p_day):.3g} partículas/día; el Reglamento (UE) "
                              "2024/1448 fija metodología de medida para agua de consumo: comparar "
                              "solo metodológicamente"}


def social_survey_analysis(df, likert_cols, group_col=None):
    """Análisis de encuestas de percepción social (escala Likert)."""
    data = df[likert_cols].apply(pd.to_numeric, errors="coerce")
    if data.empty:
        return {"test": "Percepción social", "error": "Sin respuestas"}
    cronbach = _cronbach_alpha(data)
    summary = data.describe().T[["count", "mean", "std"]].round(2)
    out = {"test": "Análisis de percepción social", "resumen": summary,
           "cronbach_alpha": cronbach,
           "interpretacion": f"Alfa de Cronbach={cronbach:.2f} "
                             f"({'consistencia aceptable' if cronbach > .7 else 'revisar consistencia del instrumento'})"}
    if group_col and group_col in df.columns:
        rows = []
        for col in likert_cols:
            groups = [g[col].dropna().values for _, g in df.groupby(group_col, observed=True)
                      if g[col].notna().sum() >= 3]
            if len(groups) >= 2:
                h, p = stats.kruskal(*groups)
                rows.append({"pregunta": col, "H": round(float(h), 2), "p": round(float(p), 4)})
        out["por_grupo"] = pd.DataFrame(rows)
    return out


def _cronbach_alpha(items):
    items = items.dropna()
    k = items.shape[1]
    if k < 2 or len(items) < 3:
        return np.nan
    var_total = items.sum(axis=1).var(ddof=1)
    var_items = items.var(ddof=1).sum()
    return float(k / (k - 1) * (1 - var_items / var_total)) if var_total > 0 else np.nan
