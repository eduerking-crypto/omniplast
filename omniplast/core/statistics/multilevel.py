"""Modelos multinivel, GEE, regresion cuantilica, mediacion y bootstrap."""
import numpy as np
import pandas as pd


def ols_regression(df, target, predictors, robust=True):
    import statsmodels.api as sm
    data = df[[target] + list(predictors)].apply(pd.to_numeric, errors="coerce").dropna()
    if len(data) < len(predictors) + 3:
        return {"modelo": "OLS", "error": "n insuficiente"}
    X = sm.add_constant(data[predictors])
    model = sm.OLS(data[target], X).fit(cov_type="HC3" if robust else None)
    coef = pd.DataFrame({"coef": model.params, "se": model.bse, "t": model.tvalues,
                         "p": model.pvalues, "ic_inf": model.conf_int()[0],
                         "ic_sup": model.conf_int()[1]})
    vif = _vif(data, predictors)
    from scipy import stats as _st
    jb_p = float(_st.jarque_bera(model.resid)[1])
    return {"modelo": "Regresión lineal múltiple (HC3)" if robust else "OLS",
            "coeficientes": coef, "r2": float(model.rsquared),
            "r2_adj": float(model.rsquared_adj), "aic": float(model.aic),
            "f_pvalue": float(model.f_pvalue), "n": int(model.nobs), "vif": vif,
            "jarque_bera_p": jb_p,
            "interpretacion": f"R²aj={model.rsquared_adj:.3f}; F-test p={model.f_pvalue:.4g}"}


def _vif(data, predictors):
    from statsmodels.stats.outliers_influence import variance_inflation_factor
    X = data[predictors].values
    if X.shape[1] < 2 or len(X) < 10:
        return pd.Series(dtype=float)
    vifs = {pred: float(variance_inflation_factor(X, i)) for i, pred in enumerate(predictors)}
    return pd.Series(vifs).sort_values(ascending=False)


def glm_counts(df, target, predictors, family="negativebinomial", offset=None):
    """GLM para conteos de MP: Poisson o Binomial Negativa (efecto de sobredispersion)."""
    import statsmodels.api as sm
    from statsmodels.discrete.discrete_model import NegativeBinomial
    data = df[[target] + list(predictors)].apply(pd.to_numeric, errors="coerce").dropna()
    if offset is not None:
        data = data.assign(_off=pd.to_numeric(df.loc[data.index, offset], errors="coerce").fillna(1))
    if len(data) < len(predictors) + 4:
        return {"modelo": "GLM", "error": "n insuficiente"}
    y = data[target].clip(lower=0)
    X = sm.add_constant(data[predictors])
    off = np.log(data["_off"].values) if offset is not None else None
    try:
        if family == "poisson":
            model = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit()
            alpha = np.nan
        else:
            nb = NegativeBinomial(y, X, loglike_method="nb2", offset=off)
            res_nb = nb.fit(disp=False)
            alpha = float(res_nb.params.get("alpha", np.nan))
            model = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=max(alpha, 1e-8)),
                           offset=off).fit()
        if not np.isfinite(model.params.values).all():
            raise ValueError("coeficientes no finitos")
    except Exception:
        model = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit()
        alpha = np.nan
        family = "poisson (fallback: NB no convergió)"
    coef = pd.DataFrame({"coef": model.params, "se": model.bse, "z": model.tvalues,
                         "p": model.pvalues, "ic_inf": model.conf_int()[0],
                         "ic_sup": model.conf_int()[1]})
    ratio = np.exp(model.params.drop("const", errors="ignore"))
    return {"modelo": f"GLM {family} (log-link)", "coeficientes": coef,
            "razones_de_cambio": ratio, "alpha_sobredispersion": alpha,
            "aic": float(model.aic), "devianza": float(model.deviance),
            "n": int(model.nobs),
            "interpretacion": ("Sobredispersión α>0 confirma necesidad de Binomial Negativa "
                               f"(α={alpha:.3f})" if family != "poisson" and alpha == alpha
                               else f"AIC={model.aic:.1f}")}


def mixedlm(df, target, predictors, group):
    """Modelo mixto con intercepto aleatorio por grupo (sitio/campana)."""
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    data = df[[target] + list(predictors) + [group]].apply(pd.to_numeric, errors="coerce")
    data[group] = df[group].astype(str)
    data = data.dropna()
    if data[group].nunique() < 2 or len(data) < 8:
        return {"modelo": "Mixto", "error": "Se requieren ≥2 grupos y n≥8"}
    formula = f"{target} ~ " + " + ".join(predictors)
    model = smf.mixedlm(formula, data, groups=data[group]).fit(reml=False)
    fe = pd.DataFrame({"coef": model.fe_params, "se": model.bse_fe, "p": model.pvalues})
    vc = float(model.cov_re.iloc[0, 0])
    resid = float(model.scale)
    icc = vc / (vc + resid) if (vc + resid) > 0 else np.nan
    return {"modelo": "Modelo lineal mixto (intercepto aleatorio)", "efectos_fijos": fe,
            "var_intercepto_grupo": vc, "var_residual": resid, "icc": icc,
            "interpretacion": f"ICC={icc:.3f} → {icc * 100:.1f}% de la varianza atribuible a {group}"}


def gee_population(df, target, predictors, group, family="poisson"):
    """GEE (ecuaciones de estimacion generalizadas) con estructura intercambiable."""
    import statsmodels.api as sm
    data = df[[target] + list(predictors) + [group]].copy()
    data[target] = pd.to_numeric(data[target], errors="coerce")
    for p in predictors:
        data[p] = pd.to_numeric(data[p], errors="coerce")
    data[group] = df[group].astype(str)
    data = data.dropna()
    if data[group].nunique() < 2:
        return {"modelo": "GEE", "error": "Se requieren ≥2 grupos (clusters)"}
    X = sm.add_constant(data[predictors])
    fam = sm.families.Poisson() if family == "poisson" else sm.families.Gaussian()
    model = sm.GEE(data[target], X, groups=data[group], family=fam,
                   cov=sm.cov_struct.Exchangeable()).fit()
    coef = pd.DataFrame({"coef": model.params, "se_robusta": model.bse, "p": model.pvalues})
    return {"modelo": f"GEE ({family}, intercambiable)", "coeficientes": coef,
            "n_clusters": int(data[group].nunique()),
            "interpretacion": "SE robustas por cluster; inferencia poblacional (no por-cluster)"}


def quantile_regression(df, target, predictor, qs=(0.25, 0.5, 0.75)):
    """Regresion cuantilica: efecto del predictor en extremos de la distribucion."""
    import statsmodels.formula.api as smf
    data = df[[target, predictor]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(data) < 10:
        return {"modelo": "Cuantílica", "error": "n≥10"}
    rows = []
    for q in qs:
        res = smf.quantreg(f"{target} ~ {predictor}", data).fit(q=q)
        slope = float(res.params[predictor])
        rows.append({"cuantil": q, "intercepto": float(res.params["Intercept"]),
                     "pendiente": slope, "p": float(res.pvalues[predictor])})
    return {"modelo": "Regresión cuantílica", "tabla": pd.DataFrame(rows),
            "interpretacion": "Compare las pendientes entre cuantiles: diferencias sugieren "
                              "efecto heterogéneo del predictor en la respuesta"}


def mediation_analysis(df, x, m, y, n_boot=2000, seed=42):
    """Mediacion de Baron-Kenny con bootstrap del efecto indirecto."""
    import statsmodels.api as sm
    rng = np.random.default_rng(seed)
    data = df[[x, m, y]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(data) < 12:
        return {"modelo": "Mediación", "error": "n≥12"}
    def _paths(d):
        a = sm.OLS(d[m], sm.add_constant(d[[x]])).fit().params[x]
        b = sm.OLS(d[y], sm.add_constant(d[[x, m]])).fit().params[m]
        c = sm.OLS(d[y], sm.add_constant(d[[x]])).fit().params[x]
        cpr = sm.OLS(d[y], sm.add_constant(d[[x, m]])).fit().params[x]
        return a, b, c, cpr
    a, b, c, cpr = _paths(data)
    indirect = []
    idx = np.arange(len(data))
    for _ in range(n_boot):
        sample = data.iloc[rng.choice(idx, len(idx), replace=True)]
        try:
            ai, bi, _, _ = _paths(sample)
            indirect.append(ai * bi)
        except Exception:
            continue
    lo, hi = np.percentile(indirect, [2.5, 97.5])
    ab = a * b
    pct = ab / c * 100 if c != 0 else np.nan
    return {"modelo": "Análisis de mediación (bootstrap 2000)", "a": float(a), "b": float(b),
            "c_total": float(c), "c_directo": float(cpr), "ab_indirecto": float(ab),
            "ic95_indirecto": (float(lo), float(hi)), "pct_mediado": float(pct),
            "interpretacion": f"Efecto indirecto a×b={ab:.4f} (IC95% bootstrap {lo:.4f} a {hi:.4f}); "
                              f"{'mediación significativa' if lo * hi > 0 else 'sin mediación significativa'}"}


def bootstrap_stat(series, stat="mean", n_boot=2000, seed=42, stat_fn=None):
    rng = np.random.default_rng(seed)
    x = pd.to_numeric(pd.Series(series), errors="coerce").dropna().values
    if len(x) < 3:
        return {"test": "Bootstrap", "error": "n≥3"}
    fn = stat_fn or {"mean": np.mean, "median": np.median, "sd": lambda v: np.std(v, ddof=1),
                     "cv": lambda v: np.std(v, ddof=1) / np.mean(v) if np.mean(v) else np.nan}.get(stat, np.mean)
    boots = np.array([fn(x[rng.integers(0, len(x), len(x))]) for _ in range(n_boot)])
    point = fn(x)
    return {"test": f"Bootstrap ({stat}, B={n_boot})", "estimado": float(point),
            "ic95": (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))),
            "se_bootstrap": float(boots.std()),
            "interpretacion": f"{stat}={point:.4g}, IC95%=[{np.percentile(boots, 2.5):.4g}, "
                              f"{np.percentile(boots, 97.5):.4g}]"}
