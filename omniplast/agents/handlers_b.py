"""Handlers de agentes D (estadistica) y E (machine learning): AG-033..062."""
import numpy as np
import pandas as pd

from ..core.statistics import (bayes, descriptive, inferential, ml, multilevel,
                               multivariate, timeseries)
from ..reporting import figures
from .handlers_a import _ok, _req


def _merged(ctx, matriz="agua"):
    mp = ctx["mp"]
    env = ctx["env"]
    if mp.empty or env.empty:
        return pd.DataFrame()
    mpm = mp[mp["matriz"] == matriz].groupby("sample_id")["abundancia"].sum().rename("mp_total")
    d = env[env["matriz"] == matriz].merge(mpm, on="sample_id", how="left")
    return d.dropna(subset=["mp_total"])


def _features(d, exclude=("mp_total", "sample_id", "lat", "lon", "anio", "profundidad_m")):
    num = d.select_dtypes("number").columns.tolist()
    return [c for c in num if c not in exclude and d[c].notna().sum() >= 10]


def ag_033(ctx):
    d = _merged(ctx)
    if d.empty:
        return _req("Se requieren datos de MP fusionados con ambiente.")
    nt = descriptive.normality_tests(d["mp_total"])
    rec = ("Distribución no normal y conteos → preferir pruebas no paramétricas o GLM "
           "Poisson/NB (nunca ANOVA directo sobre conteos)." if nt["p"].min() < .05
           else "Aproximación paramétrica admisible con cautela.")
    fig = figures.box(d, x="zona", y="mp_total", title="Distribución de MP", ytitle="part/m³")
    return _ok(f"Shapiro-Wilk p={nt.iloc[0]['p']:.4g}; sesgo={descriptive.describe(d['mp_total'])['sesgo']:.2f}. "
               f"Recomendación: {rec}", {"normalidad": nt.round(4)}, {"cajas": fig},
               ["Zar (2010) Biostatistical Analysis"])


def ag_034(ctx):
    d = _merged(ctx)
    if d.empty:
        return _req("Se requieren datos de MP.")
    st = descriptive.suggest_transformation(d["mp_total"])
    return _ok(f"Mejor transformación para MP: {st['mejor_transformacion']}. "
               "Los conteos ambientales suelen estabilizarse con log(x+1) o sqrt(x).",
               {"p_valores_shapiro": pd.DataFrame([st["p_valores"]])},
               [], ["Box & Cox (1964) J. R. Stat. Soc. B 26:211-252"])


def ag_035(ctx):
    d = _merged(ctx)
    if d.empty:
        return _req("Se requieren MP y al menos 2 zonas.")
    r = inferential.ttest_two(d, "mp_total", "zona", "impacto", "referencia")
    if "error" in r:
        return _req("Se requieren ≥2 observaciones en impacto y referencia.")
    fig = figures.box(d[d["zona"].isin(["impacto", "referencia"])], x="zona", y="mp_total",
                      title="impacto vs referencia", ytitle="part/m³")
    return _ok(r["interpretacion"], {"prueba": pd.DataFrame([r])}, {"cajas": fig},
               ["Welch (1947) Biometrika 34:28-35"])


def ag_036(ctx):
    d = _merged(ctx)
    if d.empty:
        return _req("Se requieren datos de MP.")
    r = inferential.anova_one(d, "mp_total", "zona")
    if "error" in r:
        return _req("Se requieren ≥2 grupos con n≥2.")
    tukey = inferential.tukey_hsd(d, "mp_total", "zona")
    return _ok(r["interpretacion"], {"anova": pd.DataFrame([r]), "tukey": tukey},
               [], ["Fisher (1925)", "Tukey (1949)"])


def ag_037(ctx):
    d = _merged(ctx)
    if d.empty:
        return _req("Se requieren datos de MP.")
    r = inferential.kruskal_wallis(d, "mp_total", "zona")
    dunn = inferential.dunn_posthoc(d, "mp_total", "zona")
    return _ok(r["interpretacion"], {"kruskal": pd.DataFrame([r]), "dunn": dunn.round(4)},
               [], ["Kruskal & Wallis (1952)", "Dunn (1964)"])


def ag_038(ctx):
    d = _merged(ctx)
    if d.empty:
        return _req("Se requieren datos de MP.")
    lev = descriptive.homogeneity_tests(d, "mp_total", "zona")
    gh = inferential.games_howell(d, "mp_total", "zona")
    return _ok("Games-Howell es válido bajo heterocedasticidad. "
               + ("Varianzas heterogéneas confirmadas → GH es la elección correcta."
                  if len(lev) and lev.iloc[0]["p"] < .05 else "Varianzas homogéneas; "
                  "ANOVA+Tukey también sería admisible."),
               {"levene": lev.round(4), "games_howell": gh.round(4)},
               [], ["Games & Howell (1976) J. Educ. Stat. 1:113-125"])


def ag_039(ctx):
    mp = ctx["mp"]
    per = mp[mp["matriz"] == "agua"].groupby(["sitio", "temporada"], observed=True)["abundancia"].sum().reset_index()
    piv = per.pivot_table(index="sitio", columns="temporada", values="abundancia").dropna()
    if piv.shape[1] < 2 or len(piv) < 3:
        return _req("Se requieren ≥3 sitios con ambas temporadas.", ["Zar (2010)"])
    cols = sorted(piv.columns)
    melted = piv.reset_index().melt(id_vars="sitio", value_vars=cols,
                                    var_name="temporada", value_name="abundancia")
    w = inferential.wilcoxon_paired(melted, "abundancia", "temporada", cols[0], cols[1])
    fig = figures.box(per, x="temporada", y="abundancia", title="Pareado por sitio",
                      ytitle="part/m³")
    return _ok(f"{w.get('interpretacion', w.get('error', ''))} sobre {len(piv)} sitios pareados.",
               {"pareado": piv.round(0)}, {"cajas": fig}, ["Zar (2010)"])


def ag_040(ctx):
    d = _merged(ctx)
    feats = _features(d)
    if d.empty or not feats:
        return _req("Se requieren variables numéricas ambientales.")
    r_mat, p_mat = inferential.correlation_matrix(d[["mp_total"] + feats[:12]], method="spearman")
    fig = figures.heatmap(r_mat, title="Matriz de correlaciones (Spearman)")
    top = (r_mat["mp_total"].drop("mp_total").dropna()).sort_values()
    return _ok("Variables más asociadas a MP: " +
               ", ".join(f"{k} (r={v:.2f})" for k, v in top.head(3).items()),
               {"r": r_mat.round(3), "p": p_mat.round(4)}, {"heatmap": fig}, ["Zar (2010)"])


def ag_041(ctx):
    d = _merged(ctx)
    if d.empty:
        return _req("Se requieren datos fusionados.")
    covs = [c for c in ("salinidad_psu", "turbidez_ntu", "tss_mg_l") if c in d.columns]
    if len(covs) < 1:
        return _req("Se requieren covariables (salinidad, turbidez, SST).")
    r = inferential.partial_correlation(d, "dbo_mg_l", "mp_total", covs)
    return _ok(f"{r.get('interpretacion', r.get('error', ''))}. "
               "La correlación parcial aísla el efecto de la covariable de confusión.",
               {"parcial": pd.DataFrame([r])}, [], ["Zar (2010)"])


def ag_042(ctx):
    d = _merged(ctx)
    feats = _features(d)
    if d.empty or not feats:
        return _req("Se requieren predictores ambientales.")
    r = multilevel.ols_regression(d, "mp_total", feats[:5])
    if "error" in r:
        return _req(r["error"])
    return _ok(f"{r['interpretacion']}; AIC={r['aic']:.1f}. "
               "VIF>5 indica colinealidad a corregir (retirar variable).",
               {"coeficientes": r["coeficientes"].round(4), "vif": r["vif"].round(2)},
               {}, ["MacKinnon & White (1985) J. Econometrics 21:53-70"])


def ag_043(ctx):
    d = _merged(ctx)
    if d.empty or "dbo_mg_l" not in d:
        return _req("Se requiere un predictor continuo (DBO).", ["Koenker & Bassett (1978)"])
    r = multilevel.quantile_regression(d.dropna(subset=["dbo_mg_l"]), "mp_total", "dbo_mg_l")
    fig = figures.bars(r["tabla"], "cuantil", "pendiente", title="Pendiente por cuantil",
                       ytitle="β (part/m³ por mg DBO/L)")
    return _ok(r["interpretacion"], {"cuantiles": r["tabla"].round(4)}, {"barras": fig},
               ["Koenker & Bassett (1978) Econometrica 46:33-50"])


def ag_044(ctx):
    d = _merged(ctx)
    if d.empty:
        return _req("Se requieren conteos y predictores.")
    r = multilevel.glm_counts(d, "mp_total", ["dbo_mg_l", "tss_mg_l"], family="poisson")
    if "error" in r:
        return _req(r["error"])
    return _ok(f"GLM Poisson: {r['interpretacion']}. "
               "Si la devianza/df >> 1 hay sobredispersión → usar Binomial Negativa (AG-045).",
               {"coeficientes": r["coeficientes"].round(4),
                "razones_de_cambio": r["razones_de_cambio"].round(3)},
               {}, ["McCullagh & Nelder (1989) GLM 2ª ed."])


def ag_045(ctx):
    d = _merged(ctx)
    if d.empty:
        return _req("Se requieren conteos y predictores.")
    r = multilevel.glm_counts(d, "mp_total", ["dbo_mg_l", "tss_mg_l"], family="negativebinomial")
    if "error" in r:
        return _req(r["error"])
    return _ok(f"GLM NB: {r['interpretacion']}. α>0 cuantifica la sobredispersión "
               "(agregación espacial de MP: esperada en puntos de vertimiento).",
               {"coeficientes": r["coeficientes"].round(4),
                "razones_de_cambio": r["razones_de_cambio"].round(3)},
               {}, ["Cameron & Trivedi (2013) Regression Analysis of Count Data 2ª ed."])


def ag_046(ctx):
    d = _merged(ctx)
    if d.empty or d["sitio"].nunique() < 2:
        return _req("Se requieren ≥2 sitios.", ["Zuur et al. (2009)"])
    r = multilevel.mixedlm(d, "mp_total", ["dbo_mg_l"], "sitio")
    if "error" in r:
        return _req(r["error"])
    return _ok(f"{r['interpretacion']}. El intercepto aleatorio por sitio absorbe la "
               "pseudorreplicación espacial.", {"efectos_fijos": r["efectos_fijos"].round(4)},
               {}, ["Zuur et al. (2009) Mixed Effects Models in R"])


def ag_047(ctx):
    d = _merged(ctx)
    if d.empty or d["sitio"].nunique() < 2:
        return _req("Se requieren ≥2 clusters (sitios).", ["Liang & Zeger (1986)"])
    r = multilevel.gee_population(d, "mp_total", ["dbo_mg_l"], "sitio", family="gaussian")
    if "error" in r:
        return _req(r["error"])
    return _ok(f"{r['interpretacion']}",
               {"coeficientes": r["coeficientes"].round(4)}, {},
               ["Liang & Zeger (1986) Biometrika 73:13-22"])


def ag_048(ctx):
    d = _merged(ctx)
    feats = _features(d)
    if d.empty or len(feats) < 3:
        return _req("Se requieren ≥3 variables ambientales.", ["Pearson (1901)"])
    r = multivariate.pca(d[feats].dropna())
    scores = r["scores"].copy()
    scores["zona"] = d.loc[scores.index, "zona"].values
    fig = figures.scatter(scores, x="PC1", y="PC2", color="zona",
                          title=f"PCA ambiental ({r['explicada_total']:.0f}% PC1+PC2)")
    return _ok(f"PC1 explica {r['explicada'].iloc[0]:.1f}%; Kaiser sugiere {r['kaiser_componentes']} componentes. "
               f"Variables más influyentes en PC1: "
               + ", ".join(r["loadings"]["PC1"].abs().sort_values(ascending=False).head(3).index),
               {"cargas": r["loadings"].round(3), "explicada": r["explicada"].round(2)},
               {"biplot": fig}, ["Pearson (1901)"])


def ag_049(ctx):
    d = _merged(ctx)
    feats = _features(d)
    if d.empty or not feats:
        return _req("Se requieren variables ambientales.")
    site_mean = d.groupby("sitio")[feats].mean().dropna()
    if len(site_mean) < 3:
        return _req("Se requieren ≥3 sitios.")
    r = multivariate.hierarchical_clustering(site_mean)
    return _ok(f"Cophenetic correlation={r['cophenetic']:.2f} (>0.7 = dendrograma fiel). "
               "Los pares de sitios que se unen primero comparten firma fisicoquímica.",
               {"medias_sitio": site_mean.round(2)}, {},
               ["Sokal & Rohlf (1962)"]) if "cophenetic" in r else _req("Clustering falló.")


def ag_050(ctx):
    mb = ctx["mb"]
    d = _merged(ctx)
    if mb is not None and not mb.empty and mb.shape[1] > 8:
        taxa = [c for c in mb.columns if c not in ("sample_id", "codigo", "sitio", "zona", "campana", "temporada", "matriz", "anio")]
        dist = multivariate.distance_matrix(mb.set_index("sample_id")[taxa], method="bray")
        groups = mb.set_index("sample_id")["zona"]
        per = multivariate.permanova(dist, groups, permutations=999)
        ani = multivariate.anosim(dist, groups, permutations=999)
        disp = multivariate.betadisper_groups(dist, groups)
        return _ok(f"Microbioma: {per['interpretacion']} {ani['interpretacion']}",
                   {"permanova": pd.DataFrame([per]), "anosim": pd.DataFrame([ani]),
                    "dispersion": disp.round(3)},
                   {}, ["Anderson (2001) Austral. NZ J. Stat. 43:281-300",
                        "Clarke & Gorley (2006) PRIMER v6"])
    if d.empty:
        return _req("Se requiere microbioma 16S o ambiente multivariado.",
                    ["Anderson (2001)"])
    feats = _features(d)
    dist = multivariate.distance_matrix(d.set_index("sample_id")[feats].fillna(0), method="euclidean")
    per = multivariate.permanova(dist, d.set_index("sample_id")["zona"], permutations=999)
    return _ok(f"PERMANOVA ambiental: {per['interpretacion']}",
               {"permanova": pd.DataFrame([per])}, {}, ["Anderson (2001)"])


def ag_051(ctx):
    mb = ctx["mb"]
    if mb is None or mb.empty:
        return _req("Se requieren abundancias de taxa.", ["Clarke (1993)"])
    taxa = [c for c in mb.columns if c not in ("sample_id", "codigo", "sitio", "zona", "campana", "temporada", "matriz", "anio")]
    sim = multivariate.simper(mb[taxa], mb["zona"])
    first = list(sim["por_par"].keys())[0]
    top = sim["por_par"][first].head(8)[["variable", "contribucion_pct", "acumulado_pct"]]
    return _ok(f"SIMPER {first}: disimilitud media {sim['disimilitud_media_pct'].get(first, np.nan):.1f}%. "
               f"Principales discriminadores: {', '.join(top['variable'].head(3))}",
               {"simper": top.round(2)}, {}, ["Clarke (1993) Mar. Ecol. Prog. Ser. 92:205-219"])


def ag_052(ctx):
    mb = ctx["mb"]
    if mb is None or mb.empty:
        return _req("Se requieren abundancias de taxa.", ["Kruskal (1964)"])
    taxa = [c for c in mb.columns if c not in ("sample_id", "codigo", "sitio", "zona", "campana", "temporada", "matriz", "anio")]
    dist = multivariate.distance_matrix(mb.set_index("sample_id")[taxa], method="bray")
    nmds = multivariate.nmds(dist, n_init=8)
    coords = nmds["coordenadas"].copy()
    coords["zona"] = mb.set_index("sample_id").loc[coords.index, "zona"].values
    fig = figures.scatter(coords, x="NMDS1", y="NMDS2", color="zona",
                          title=f"NMDS (stress={nmds['stress']:.3f}, {nmds['calidad_ajuste']})")
    return _ok(f"Stress={nmds['stress']:.3f} → ajuste {nmds['calidad_ajuste']}. {nmds['regla']}",
               {"coordenadas": coords.round(3)}, {"nmds": fig},
               ["Kruskal (1964) Psychometrika 29:115-129"])


def ag_053(ctx):
    d = _merged(ctx)
    feats = _features(d)
    if d.empty or len(feats) < 3:
        return _req("Se requieren ≥3 predictores.", ["Breiman (2001)"])
    d = d.assign(impacto=(d["zona"] == "impacto").astype(int))
    r = ml.random_forest_classifier(d, feats, "impacto")
    if "error" in r:
        return _req(r["error"])
    return _ok(r["interpretacion"], {"metricas": pd.DataFrame([r["metrics"]]),
                                     "confusion": r["matriz_confusion"],
                                     "importancias": r["importancias"].round(4)},
               {}, ["Breiman (2001) Mach. Learn. 45:5-32"],
               ["Reportar siempre CV, no solo test: evitaroptimismo por división única"])


def ag_054(ctx):
    d = _merged(ctx)
    feats = _features(d)
    if d.empty or len(feats) < 2:
        return _req("Se requieren predictores ambientales.", ["Breiman (2001)"])
    r = ml.random_forest_regressor(d, feats, "mp_total")
    if "error" in r:
        return _req(r["error"])
    return _ok(r["interpretacion"], {"metricas": pd.DataFrame([r["metrics"]]),
                                     "importancias": r["importancias"].round(4)},
               {}, ["Breiman (2001)"])


def ag_055(ctx):
    d = _merged(ctx)
    feats = _features(d)
    if d.empty or len(feats) < 2:
        return _req("Se requieren predictores.", ["Breiman (2001)"])
    r = ml.random_forest_regressor(d, feats, "mp_total")
    if "importancias" not in r:
        return _req(r.get("error", "Modelo no ajustable."))
    fig = figures.bars(r["importancias"].head(8).rename("importancias")
                       .rename_axis("variable").reset_index(),
                       "variable", "importancias", horizontal=True,
                       title="Importancia por permutación", ytitle="Δ R²")
    return _ok("Ranking de importancia por permutación (30 repeticiones): variables top "
               f"{', '.join(r['importancias'].head(3).index)}",
               {"importancias": r["importancias"].round(4)}, {"barras": fig},
               ["Breiman (2001)"])


def ag_056(ctx):
    d = _merged(ctx)
    feats = _features(d)
    if d.empty or len(feats) < 2:
        return _req("Se requieren predictores.", ["Cortes & Vapnik (1995)"])
    d = d.assign(impacto=(d["zona"] == "impacto").astype(int))
    r = ml.svm_classifier(d, feats, "impacto")
    if "error" in r:
        return _req(r["error"])
    return _ok(r["interpretacion"], {"metricas": pd.DataFrame([r["metrics"]]),
                                     "confusion": r["matriz_confusion"]},
               {}, ["Cortes & Vapnik (1995) Mach. Learn. 20:273-297"])


def ag_057(ctx):
    d = _merged(ctx)
    feats = _features(d)
    if d.empty or len(feats) < 2:
        return _req("Se requieren variables ambientales.", ["Rousseeuw (1987)"])
    site_mean = d.groupby("sitio")[feats].mean().dropna()
    r = multivariate.kmeans_clustering(site_mean, k_min=2, k_max=min(6, len(site_mean)))
    fig = figures.bars(r["tabla"], "k", "silueta", title="Silueta por k", ytitle="coeficiente")
    return _ok(f"k óptimo por silueta: {r['tabla'].sort_values('silueta', ascending=False).iloc[0]['k']}. "
               "Los clusters son hipótesis de zonas homogéneas de gestión.",
               {"scan": r["tabla"].round(3)}, {"barras": fig},
               ["Rousseeuw (1987) J. Comput. Appl. Math. 20:53-65"])


def ag_058(ctx):
    d = _merged(ctx)
    feats = _features(d)
    if d.empty or len(feats) < 2:
        return _req("Se requieren variables.", ["Ester et al. (1996)"])
    r = ml.dbscan_outliers(d.reset_index(drop=True), feats, eps=1.8)
    return _ok(r["interpretacion"],
               {"estaciones": r["tabla"][["sitio", "zona", "cluster_dbscan"]].drop_duplicates()},
               {}, ["Ester et al. (1996) KDD 96:226-231"])


def ag_059(ctx):
    d = _merged(ctx)
    feats = _features(d)
    if d.empty or len(feats) < 3:
        return _req("Se requieren ≥3 variables.", ["Chandola et al. (2009)"])
    r = ml.pca_anomaly(d.reset_index(drop=True), feats)
    anom = r["tabla"][r["tabla"]["anomalia"]]
    return _ok(r["interpretacion"] + ". Las muestras atípicas pueden ser eventos de descarga no documentados.",
               {"anomalias": anom[["codigo", "sitio", "zona", "mp_total"]].head(15)},
               {}, ["Chandola et al. (2009) ACM Comput. Surv. 41:1-58"])


def ag_060(ctx):
    d = _merged(ctx)
    feats = _features(d)
    if d.empty or len(feats) < 2:
        return _req("Se requieren predictores.", ["Friedman (2001)"])
    d = d.assign(riesgo_alto=(d["mp_total"] > d["mp_total"].quantile(.75)).astype(int))
    r = ml.gradient_boosting(d, feats, "riesgo_alto", kind="clasificacion")
    return _ok(f"Gradient Boosting: {r['metrics']}. Predice la condición de carga alta "
               "(P75) a partir del ambiente; útil como alerta temprana.",
               {"importancias": r.get("importancias", pd.Series(dtype=float)).round(4)},
               {}, ["Friedman (2001) Ann. Statist. 29:1189-1232"])


def ag_061(ctx):
    d = _merged(ctx)
    if d.empty:
        return _req("Se requieren datos de MP.", ["Efron (1979)"])
    rows = []
    for sitio, g in d.groupby("sitio"):
        b = multilevel.bootstrap_stat(g["mp_total"], stat="median", n_boot=3000)
        if "ic95" in b:
            rows.append({"sitio": sitio, "mediana": b["estimado"],
                         "ic_inf": b["ic95"][0], "ic_sup": b["ic95"][1]})
    tab = pd.DataFrame(rows)
    fig = figures.bars(tab, "sitio", "mediana", title="Mediana de MP con IC95 bootstrap",
                       ytitle="part/m³")
    fig.add_scatter(x=tab["sitio"], y=tab["ic_sup"], mode="markers", name="IC sup",
                    marker=dict(symbol="line-ns", color="#0E2433"))
    return _ok("Intervalos bootstrap (B=3000) de la mediana por sitio: "
               + "; ".join(f"{r.sitio} {r.ic_inf:.0f}-{r.ic_sup:.0f}" for r in tab.head(4).itertuples()),
               {"bootstrap": tab.round(1)}, {"barras": fig}, ["Efron (1979) Ann. Statist. 7:1-26"])


def ag_062(ctx):
    d = _merged(ctx)
    feats = _features(d)
    if d.empty or len(feats) < 2:
        return _req("Se requieren predictores.", ["Hanley & McNeil (1982)"])
    d = d.assign(riesgo_alto=(d["mp_total"] > d["mp_total"].median()).astype(int))
    r = ml.roc_data(d, feats, "riesgo_alto")
    if "error" in r:
        return _req(r["error"])
    fig = figures.line(pd.DataFrame({"fpr": r["fpr"], "tpr": r["tpr"]}), x="fpr", y="tpr",
                       title=f"Curva ROC (AUC={r['auc']:.3f})", xtitle="1 - especificidad",
                       ytitle="sensibilidad")
    return _ok(f"AUC={r['auc']:.3f}; {r['interpretacion']}. "
               "OR>1 del predictor = factor de riesgo de carga alta.",
               {"odds_ratios": r["odds_ratios"].round(3)}, {"roc": fig},
               ["Hanley & McNeil (1982) Radiology 143:29-36"])


HANDLERS_B = {f"AG-{i:03d}": f for i, f in enumerate([
    ag_033, ag_034, ag_035, ag_036, ag_037, ag_038, ag_039, ag_040, ag_041, ag_042,
    ag_043, ag_044, ag_045, ag_046, ag_047, ag_048, ag_049, ag_050, ag_051, ag_052,
    ag_053, ag_054, ag_055, ag_056, ag_057, ag_058, ag_059, ag_060, ag_061, ag_062],
    start=33)}
