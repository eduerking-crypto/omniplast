"""Handlers de agentes A (fisicoquimica), B (sedimentos) y C (nucleo MP): AG-001..032."""
import numpy as np
import pandas as pd

from ..core import indices, sources, uncertainty
from ..core.geo import spatial, trajectory
from ..core.statistics import inferential, multilevel, multivariate, timeseries
from ..reporting import figures
from ..config import DEMO_HULL, MP_TIERS, POLYMER_DENSITY  # noqa: F401

REFS_QAQC = ["Koelmans et al. (2019) Environ. Sci. Technol. 53:9615-9623"]


def _req(msg, refs=None, pasos=None):
    return {"estado": "requiere", "resumen": msg, "tablas": {}, "figuras": {},
            "refs": refs or [], "pasos": pasos or []}


def _ok(resumen, tablas=None, figuras=None, refs=None, pasos=None):
    return {"estado": "ok", "resumen": resumen, "tablas": tablas or {},
            "figuras": figuras or {}, "refs": refs or [], "pasos": pasos or []}


def _agua(ctx):
    env = ctx["env"]
    return env[env["matriz"] == "agua"].copy() if not env.empty else env


def _sed(ctx):
    env = ctx["env"]
    return env[env["matriz"] == "sedimento"].copy() if not env.empty else env


def ag_001(ctx):
    agua = _agua(ctx)
    if agua.empty or "temperatura_c" not in agua:
        return _req("Se requieren mediciones de temperatura en agua.", ["Hartmann et al. (2019)"])
    res = inferential.kruskal_wallis(agua, "temperatura_c", "zona")
    prof = agua[agua["profundidad_m"] > 0]
    fig = figures.box(agua, x="zona", y="temperatura_c", title="Temperatura por zona",
                      ytitle="°C", color="temporada")
    estrat = ""
    if not prof.empty and agua["profundidad_m"].nunique() > 1:
        corr = inferential.correlation_test(agua["profundidad_m"], agua["temperatura_c"], method="spearman")
        estrat = f" Perfil vertical: r_s={corr.get('r', np.nan):.2f} (p={corr.get('p', np.nan):.3g})."
    return _ok(f"Mediana general {agua['temperatura_c'].median():.1f} °C. {res['interpretacion']}.{estrat}",
               {"por_zona": agua.groupby("zona")["temperatura_c"].describe().round(2)},
               {"cajas": fig}, ["Hartmann et al. (2019)"])


def ag_002(ctx):
    agua = _agua(ctx)
    if agua.empty or "od_mg_l" not in agua:
        return _req("Se requieren mediciones de oxígeno disuelto (od_mg_l).", ["US EPA (1986)"])
    T = agua["temperatura_c"].fillna(28)
    S = agua["salinidad_psu"].fillna(0)
    do_sat = (14.652 - 0.41022 * T + 0.007991 * T ** 2 - 7.7774e-5 * T ** 3) * (1 - 0.0097 * S)
    agua = agua.assign(saturacion_pct=agua["od_mg_l"] / do_sat.clip(lower=1) * 100)
    hipoxia = (agua["od_mg_l"] < 2).mean() * 100
    res = inferential.kruskal_wallis(agua, "od_mg_l", "zona")
    fig = figures.box(agua, x="zona", y="saturacion_pct", title="Saturación de O₂ por zona",
                      ytitle="% saturación")
    return _ok(f"{hipoxia:.0f}% de muestras en hipoxia (<2 mg/L). {res['interpretacion']}",
               {"por_zona": agua.groupby("zona")["saturacion_pct"].describe().round(1)},
               {"cajas": fig},
               ["Benson & Krause (1984) USGS", "US EPA (1986) Quality Criteria for Water"],
               ["Sitios con saturación <60% requieren acción de restauración"])


def ag_003(ctx):
    agua = _agua(ctx)
    if agua.empty or "dbo_mg_l" not in agua:
        return _req("Se requieren DBO y/o DQO en agua.", ["Sawyer et al. (2003)"])
    agua = agua.assign(ratio=agua["dqo_mg_l"] / agua["dbo_mg_l"].replace(0, np.nan))
    res = inferential.kruskal_wallis(agua, "dbo_mg_l", "zona")
    fig = figures.box(agua, x="zona", y="dbo_mg_l", title="DBO₅ por zona", ytitle="mg/L")
    return _ok(f"DBO mediana={agua['dbo_mg_l'].median():.0f} mg/L; razón DQO/DBO mediana="
               f"{agua['ratio'].median():.1f} (>2.5 sugiere fracción recalcitrante). {res['interpretacion']}",
               {"por_zona": agua.groupby("zona")[["dbo_mg_l", "dqo_mg_l", "ratio"]].describe().round(1)},
               {"cajas": fig}, ["Sawyer et al. (2003) Chemistry for Environmental Engineering"],
               ["Razón DQO/DBO > 2.5 indica compuestos poco biodegradables (posible aditivos de plástico)"])


def ag_004(ctx):
    agua = _agua(ctx)
    for c in ("nitrato_mg_l", "nitrito_mg_l", "amonio_mg_l", "fosfato_mg_l"):
        if c not in agua:
            return _req(f"Falta la variable {c} para calcular DIN y N:P.",
                        ["Redfield et al. (1963)"])
    agua = agua.assign(din=agua[["nitrato_mg_l", "nitrito_mg_l", "amonio_mg_l"]].sum(axis=1, min_count=1))
    agua = agua.assign(np_ratio=(agua["din"] * 14.007 / 62) / (agua["fosfato_mg_l"].replace(0, np.nan) / 3.066),
                       din_umol=agua["din"] * 16.13, dip_umol=agua["fosfato_mg_l"] * 32.06)
    res = inferential.kruskal_wallis(agua.dropna(subset=["din"]), "din", "zona")
    fig = figures.box(agua, x="zona", y="np_ratio", title="Ratio N:P (molar) vs Redfield 16:1",
                      ytitle="N:P")
    fig.add_hline(y=16, line=dict(color="#FF6B6B", dash="dot"), annotation_text="Redfield 16:1")
    limit = agua["np_ratio"].median()
    estado = "limitación por P" if limit > 20 else "limitación por N" if limit < 10 else "balanceado"
    return _ok(f"DIN mediana={agua['din'].median():.2f} mg/L; N:P mediano={limit:.1f} → {estado}. "
               f"{res['interpretacion']}",
               {"por_zona": agua.groupby("zona")[["din", "fosfato_mg_l", "np_ratio"]].median().round(2)},
               {"cajas": fig}, ["Redfield et al. (1963)"],
               ["N:P > 20 favorece cianobacterias; revisar fuentes de fósforo urbano"])


def ag_005(ctx):
    agua = _agua(ctx)
    if agua.empty or "salinidad_psu" not in agua:
        return _req("Se requieren mediciones de salinidad.", ["Bowden (1980)"])
    res = inferential.kruskal_wallis(agua, "salinidad_psu", "zona")
    med = agua.groupby("zona")["salinidad_psu"].median()
    tipo = "estratificado/hipersalino según estación" if med.max() - med.min() > 15 else "mixohalino bien mezclado"
    fig = figures.box(agua, x="zona", y="salinidad_psu", title="Salinidad por zona", ytitle="PSU")
    return _ok(f"Rango salino {agua['salinidad_psu'].min():.0f}-{agua['salinidad_psu'].max():.0f} PSU "
               f"→ sistema {tipo}. {res['interpretacion']}",
               {"por_zona": agua.groupby("zona")["salinidad_psu"].describe().round(1)},
               {"cajas": fig}, ["Bowden (1980) Estuarine Oceanography"])


def ag_006(ctx):
    agua = _agua(ctx)
    if agua.empty or "turbidez_ntu" not in agua:
        return _req("Se requieren turbidez y SST.", ["Kirk (2011)"])
    corr = inferential.correlation_test(agua["turbidez_ntu"], agua.get("tss_mg_l"), method="spearman")
    kd = 1.7 / (agua["turbidez_ntu"] + 1)
    zeu = np.log10(0.01) / (-kd.clip(lower=.1))
    zeu = zeu.clip(0.05, 30)
    fig = figures.scatter(agua, x="turbidez_ntu", y="tss_mg_l", trendline=True,
                          title="Turbidez vs SST", xtitle="NTU", ytitle="mg/L")
    return _ok(f"Kd≈{kd.median():.1f} m⁻¹; zona eufólica estimada ≈{zeu.median():.1f} m "
               f"(fotoinhibición superficial probable si <1.5 m). Correlación turbidez-SST: "
               f"{corr.get('interpretacion', '')}",
               {"por_zona": agua.groupby("zona")[["turbidez_ntu", "tss_mg_l"]].median().round(1)},
               {"dispersion": fig}, ["Kirk (2011) Light and Photosynthesis in Aquatic Ecosystems"],
               ["Ze champagne: fórmula empírica local calibrable"])


def ag_007(ctx):
    agua = _agua(ctx)
    if agua.empty:
        return _req("Se requieren parámetros fisicoquímicos en agua.", ["Brown et al. (1970)"])
    rows = []
    for sitio, g in agua.groupby("sitio"):
        means = g.select_dtypes("number").mean(numeric_only=True).to_dict()
        w = indices.wqi_nsf(means)
        rows.append({"sitio": sitio, "wqi": w["wqi"], "clase": w["clase"],
                     "n_parametros": w.get("n_parametros", 0)})
    tab = pd.DataFrame(rows).dropna(subset=["wqi"]).sort_values("wqi")
    fig = figures.bars(tab, "sitio", "wqi", title="WQI-NSF por estación", ytitle="WQI (0-100)")
    prom = tab["wqi"].mean()
    return _ok(f"WQI promedio del sistema: {prom:.1f} → clase '{tab['clase'].mode().iloc[0] if len(tab) else '—'}'. "
               f"Estación peor: {tab.iloc[0]['sitio']} ({tab.iloc[0]['wqi']:.1f})",
               {"wqi_por_sitio": tab}, {"barras": fig},
               ["Brown et al. (1970) Water Quality Index, NSF"],
               ["Si faltan >3 parámetros del WQI, reportar el índice con la lista de ausentes"])


def ag_008(ctx):
    agua = _agua(ctx)
    if agua.empty or "clorofila_a_ug_l" not in agua:
        return _req("Se requiere clorofila a (y fósforo total si existe).", ["Carlson (1977)"])
    tp_ug = agua.get("fosfato_mg_l", pd.Series(np.nan, index=agua.index)) * 3.066 * 1000
    rows = []
    for sitio, g in agua.groupby("sitio"):
        t = indices.tsi_carlson(chla=g["clorofila_a_ug_l"].median(),
                                 tp_ug=tp_ug.loc[g.index].median())
        rows.append({"sitio": sitio, "tsi": t["tsi"], "clase": t["clase"]})
    tab = pd.DataFrame(rows).dropna().sort_values("tsi")
    fig = figures.bars(tab, "sitio", "tsi", title="Estado trófico TSI por estación", ytitle="TSI")
    return _ok(f"TSI medio={tab['tsi'].mean():.0f} → '{tab['clase'].mode().iloc[0] if len(tab) else '—'}'",
               {"tsi_por_sitio": tab}, {"barras": fig},
               ["Carlson (1977) Limnol. Oceanogr. 22:361-369"])


def ag_009(ctx):
    agua = _agua(ctx)
    if agua.empty or "e_coli_nmp" not in agua:
        return _req("Se requieren coliformes/E. coli.", ["US EPA (2012)"])
    exc = agua.assign(exc_ecoli=agua["e_coli_nmp"] > 410,
                      exc_colif=agua.get("coliformes_totales_nmp", np.nan) > 1000)
    pct_ecoli = exc["exc_ecoli"].mean() * 100
    pct_colif = exc["exc_colif"].mean() * 100 if exc["exc_colif"].notna().any() else np.nan
    fig = figures.box(agua, x="zona", y="e_coli_nmp", title="E. coli por zona",
                      ytitle="NMP/100 mL")
    return _ok(f"{pct_ecoli:.0f}% de muestras superan el STV de E. coli (410 NMP/100 mL, contacto primario). "
               f"{pct_colif:.0f}% superan 1000 NMP de coliformes totales (referencial).",
               {"por_zona": exc.groupby("zona")[["e_coli_nmp", "coliformes_totales_nmp"]].agg(["median", "max"]).round(0)},
               {"cajas": fig},
               ["US EPA (2012) Recreational Water Quality Criteria", "OMS (2003) Guías de agua recreativa"],
               ["Superar STV implica riesgo gastroenteritis >36/1000 bañistas"])


def ag_010(ctx):
    agua = _agua(ctx)
    if agua.empty:
        return _req("Se requieren datos fisicoquímicos.")
    cortes = {"od_mg_l": ("<", 4), "dbo_mg_l": (">", 12), "turbidez_ntu": (">", 50),
              "fosfato_mg_l": (">", 0.5), "nitrato_mg_l": (">", 2.0), "e_coli_nmp": (">", 410)}
    rows = []
    for param, (op, ref) in cortes.items():
        if param not in agua:
            continue
        if op == "<":
            exc = (agua[param] < ref).mean() * 100
        else:
            exc = (agua[param] > ref).mean() * 100
        rows.append({"parametro": param, "pct_excedencia": round(exc, 1), "corte": f"{op} {ref}"})
    tab = pd.DataFrame(rows)
    score = 100 - tab["pct_excedencia"].mean() if len(tab) else np.nan
    fig = figures.radar(list(tab["parametro"]), list(100 - tab["pct_excedencia"]),
                        title="Cumplimiento por parámetro (%)")
    return _ok(f"Scorecard fisicoquímico global: {score:.0f}/100 de cumplimiento",
               {"excedencias": tab.sort_values("pct_excedencia", ascending=False)},
               {"radar": fig}, ["Elaboración OMNIPLAST (cortes referenciales verificables)"])


def ag_011(ctx):
    sed, mp = _sed(ctx), ctx["mp"]
    if sed.empty or mp[mp["matriz"] == "sedimento"].empty or "arena_pct" not in sed:
        return _req("Se requieren granulometría y MP en sedimento.", ["Besseling et al. (2017)"])
    mps = mp[mp["matriz"] == "sedimento"].groupby("sample_id")["abundancia"].sum()
    sed = sed.merge(mps.rename("mp_total"), on="sample_id", how="left").dropna(subset=["mp_total"])
    corr = inferential.correlation_test(sed["arena_pct"], sed["mp_total"], method="spearman")
    fig = figures.scatter(sed, x="arena_pct", y="mp_total", trendline=True,
                          title="Arena vs MP en sedimento", xtitle="% arena", ytitle="part/kg")
    return _ok(f"Correlación arena-MP: {corr.get('interpretacion', '')}. "
               "Los MP se acumulan preferentemente en finos (mayor superficie y MO).",
               {"por_sitio": sed.groupby("sitio")[["arena_pct", "mp_total"]].median().round(1)},
               {"dispersion": fig}, ["Besseling et al. (2017) Environ. Pollut. 218:839-845"])


def ag_012(ctx):
    sed, mp = _sed(ctx), ctx["mp"]
    if sed.empty or mp[mp["matriz"] == "sedimento"].empty or "toc_pct" not in sed:
        return _req("Se requieren TOC y MP en sedimento.", ["Besseling et al. (2017)"])
    mps = mp[mp["matriz"] == "sedimento"].groupby("sample_id")["abundancia"].sum()
    sed = sed.merge(mps.rename("mp_total"), on="sample_id", how="left").dropna(subset=["mp_total"])
    corr = inferential.correlation_test(sed["toc_pct"], sed["mp_total"], method="spearman")
    reg = multilevel.ols_regression(sed, "mp_total", ["toc_pct"], robust=True)
    fig = figures.scatter(sed, x="toc_pct", y="mp_total", trendline=True,
                          title="Materia orgánica vs MP", xtitle="COT %", ytitle="part/kg")
    return _ok(f"Correlación TOC-MP: {corr.get('interpretacion', '')}. OLS: {reg.get('interpretacion', '')}",
               {"por_sitio": sed.groupby("sitio")[["toc_pct", "mp_total"]].median().round(1)},
               {"dispersion": fig}, ["Besseling et al. (2017)"])


def ag_013(ctx):
    mp = ctx["mp"]
    mpa = mp[(mp["matriz"] == "agua") & mp["profundidad_m"].notna()]
    if mpa.empty:
        return _req("Se requieren perfiles con profundidad (0/1.5/3 m).", ["Kowalski et al. (2016)"])
    per_sample = mpa.groupby(["sample_id", "profundidad_m", "sitio"], observed=True)["abundancia"].sum().reset_index()
    corr = inferential.correlation_test(per_sample["profundidad_m"], per_sample["abundancia"], method="spearman")
    fig = figures.box(per_sample, x="profundidad_m", y="abundancia", title="MP por profundidad",
                      ytitle="part/m³")
    dens = mpa.assign(densidad=mpa["polimero"].map(POLYMER_DENSITY)).dropna(subset=["densidad"])
    flot = dens[dens["densidad"] < 1000].groupby("profundidad_m")["abundancia"].mean()
    return _ok(f"Perfil vertical: {corr.get('interpretacion', '')}. "
               f"Polímeros flotantes en superficie vs 3 m: {flot.get(0.0, np.nan):.0f} vs {flot.get(3.0, np.nan):.0f} part/m³.",
               {"perfil": per_sample.groupby("profundidad_m")["abundancia"].describe().round(0)},
               {"cajas": fig}, ["Kowalski et al. (2016) Environ. Sci. Technol. 50:6825-6834"])


def ag_014(ctx):
    mp = ctx["mp"]
    mps = mp[mp["matriz"] == "sedimento"].groupby("sample_id")["abundancia"].sum()
    if mps.empty:
        return _req("Se requieren MP en sedimento.", ["Elaboración OMNIPLAST"])
    p = ctx["params"]
    area_m2 = p.get("area_km2", 20) * 1e6
    rate_cm_yr = p.get("sed_rate_cm_yr", 0.5)
    def flux(conc, rate, bulk=1500.0):
        inventory = conc * bulk * 0.05
        return inventory * (rate / 5) * area_m2 / 365
    mc = uncertainty.monte_carlo(lambda conc, rate: flux(conc, rate),
                                {"conc": ("lognormal", float(np.log(mps.mean())), float(np.log(mps).std())),
                                 "rate": ("triangular", rate_cm_yr * .6, rate_cm_yr, rate_cm_yr * 1.4)},
                                n=8000)
    fig = figures.bars(pd.DataFrame({"escenario": ["P2.5", "mediana", "P97.5"],
                                     "flujo_part_dia": [mc["salida"]["p2.5"], mc["salida"]["p50"],
                                                        mc["salida"]["p97.5"]]}),
                        "escenario", "flujo_part_dia", title="Flujo de MP al sedimento",
                        ytitle="partículas/día")
    return _ok(f"Flujo mediano estimado: {mc['salida']['p50']:.3g} partículas/día al sedimento "
               f"(P2.5-P97.5: {mc['salida']['p2.5']:.3g} - {mc['salida']['p97.5']:.3g}). "
               "SUPUESTOS: densidad aparente 1.5 g/cm³, capa activa 5 cm, "
               f"área {p.get('area_km2', 20)} km² y tasa {rate_cm_yr} cm/año.",
               {"mc_percentiles": pd.DataFrame([mc["salida"]])}, {"barras": fig},
               ["Elaboración OMNIPLAST (inventario + tasa de sedimentación)"],
               ["Calibrar con ²¹⁰Pb/137Cs para tasa de sedimentación real"])


def ag_015(ctx):
    mp = ctx["mp"]
    mps = mp[mp["matriz"] == "sedimento"]
    if mps.empty:
        return _req("Se requieren MP en sedimento.", ["Hakanson (1980)"])
    por_sitio = mps.groupby("sitio")["abundancia"].sum()
    bg = por_sitio.min()
    rows = []
    for sitio, total in por_sitio.items():
        cf = indices.contamination_factor(total, bg)
        rows.append({"sitio": sitio, "mp_total": total, "CF": cf, "clase_CF": indices.cf_classify(cf),
                     "PLI": indices.pli_site([cf])})
    tab = pd.DataFrame(rows)
    fig = figures.bars(tab, "sitio", "CF", title="Factor de contaminación por sitio", ytitle="CF")
    return _ok(f"Fondo de referencia (sitio mínimo): {bg:.0f} part/kg. PLI del sistema: "
               f"{np.exp(np.mean(np.log(tab['PLI'].clip(lower=1e-9)))):.2f} "
               f"({'deterioro de calidad' if np.exp(np.mean(np.log(tab['PLI'].clip(lower=1e-9)))) > 1 else 'sin deterioro'}).",
               {"pli_por_sitio": tab}, {"barras": fig},
               ["Hakanson (1980) Hydrobiologia 91:832-848", "Tomlinson et al. (1980)"])


def ag_016(ctx):
    sed = _sed(ctx)
    if sed.empty or "toc_pct" not in sed:
        return _req("Se requieren TOC y textura de sedimento.", ["Berner (1981)"])
    sed = sed.assign(redox_proxy=np.clip(sed["toc_pct"] * (1 + sed.get("limo_arcilla_pct", 50) / 200), 0, 100))
    fig = figures.bars(sed.groupby("sitio")["redox_proxy"].mean().reset_index(),
                       "sitio", "redox_proxy", title="Proxy de demanda de O₂ del sedimento",
                       ytitle="proxy (a.u.)")
    return _ok("Proxy redox: MO×finos. Valores altos implican hipoxia bentónica y sulfuros "
               "(H2S) en summer stagnation; medir Eh directamente para confirmar.",
               {"por_sitio": sed.groupby("sitio")[["toc_pct", "redox_proxy"]].mean().round(2)},
               {"barras": fig}, ["Berner (1981) Early Diagenesis"])


def ag_017(ctx):
    sed = _sed(ctx)
    if sed.empty or "arena_pct" not in sed:
        return _req("Se requieren porcentajes de arena y finos.", ["Folk (1954)"])
    sed = sed.assign(textura=np.where(sed["arena_pct"] >= 90, "Arena",
                              np.where(sed["arena_pct"] >= 50, "Arena fangosa",
                                       np.where(sed["arena_pct"] >= 10, "Fangoso arenoso", "Fango (mudstone)"))))
    counts = sed["textura"].value_counts()
    fig = figures.pie(list(counts.index), list(counts.values), title="Clases texturales (Folk)")
    return _ok(f"Distribución textural: {counts.to_dict()}",
               {"por_sitio": sed.groupby("sitio")[["arena_pct", "limo_arcilla_pct", "textura"]].first()},
               {"torta": fig}, ["Folk (1954) Petrology of Sedimentary Rocks"])


def ag_018(ctx):
    sed, mp = _sed(ctx), ctx["mp"]
    mps = mp[mp["matriz"] == "sedimento"].groupby("sample_id")["abundancia"].sum()
    if sed.empty or mps.empty:
        return _req("Se requieren variables de sedimento y MP.")
    sed = sed.merge(mps.rename("mp_total"), on="sample_id", how="left").dropna(subset=["mp_total"])
    num = sed.select_dtypes("number")
    cols = [c for c in num.columns if c not in ("sample_id", "lat", "lon", "mp_total")]
    rows = []
    for c in cols:
        if sed[c].notna().sum() < 5:
            continue
        r = inferential.correlation_test(sed[c], sed["mp_total"], method="spearman")
        rows.append({"variable": c, "r": r.get("r"), "p": r.get("p")})
    tab = pd.DataFrame(rows).sort_values("p")
    if not tab.empty:
        r_mat, p_mat = inferential.correlation_matrix(sed[["mp_total"] + tab.head(4)["variable"].tolist()],
                                                      method="spearman")
        fig = figures.heatmap(r_mat, title="Correlaciones MP×sedimento (Spearman)",
                              labels=("variable", "variable"))
    else:
        fig = None
    figs = {"heatmap": fig} if fig is not None else {}
    return _ok(f"Asociaciones MP-sedimento: {len(tab)} variables probadas; "
               f"top: {tab.iloc[0]['variable']} (p={tab.iloc[0]['p']:.3g})" if len(tab) else "Sin variables suficientes",
               {"correlaciones": tab.round(4)}, figs)


def ag_019(ctx):
    mp = ctx["mp"]
    mpa = mp[mp["matriz"] == "agua"].groupby(["sample_id", "sitio", "zona"], observed=True)["abundancia"].sum().reset_index()
    if mpa.empty:
        return _req("Se requieren registros de MP en agua.", ["Hartmann et al. (2019)"])
    res = inferential.kruskal_wallis(mpa, "abundancia", "zona")
    tab = mpa.groupby("zona")["abundancia"].describe()[["count", "mean", "50%", "max"]].round(0)
    fig = figures.box(mpa, x="zona", y="abundancia", title="MP en agua por zona", ytitle="part/m³")
    tier = indices.mp_tier(mpa["abundancia"].median(), "agua")
    return _ok(f"Mediana del sistema: {mpa['abundancia'].median():.0f} part/m³ (carga {tier}). "
               f"{res['interpretacion']}",
               {"por_zona": tab}, {"cajas": fig},
               ["Hartmann et al. (2019) ES&T 53:1039-1047", "Cortes de MP_TIERS en config.py (referenciales)"])


def ag_020(ctx):
    mp = ctx["mp"]
    if mp.empty:
        return _req("Se requieren registros con banda de tamaño.", ["Kooi et al. (2016)"])
    bands = mp.assign(banda=pd.cut(mp["tamano_max_um"], [0, 100, 300, 500, 1000, 5000],
                                  labels=["1-100", "100-300", "300-500", "500-1000", "1000-5000"]))
    counts = bands.groupby("banda", observed=True)["abundancia"].sum().reset_index()
    total = counts["abundancia"].sum()
    counts["pct"] = counts["abundancia"] / total * 100
    frac300 = counts.loc[counts["banda"].astype(str).isin(["1-100", "100-300"]), "pct"].sum()
    fig = figures.bars(counts, "banda", "abundancia", title="Distribución de tamaños",
                       ytitle="partículas", xtitle="µm")
    return _ok(f"Fracción ≤300 µm: {frac300:.1f}% del total → rango biodisponible para zooplancton "
               "y filtradores (riesgo de transferencia trófica).",
               {"bandas": counts.round(1)}, {"barras": fig},
               ["Kooi et al. (2016) Environ. Res. Lett. 11:104013"],
               ["Complementar con conteo por microscopía estereoscópica en malla <300 µm"])


def ag_021(ctx):
    mp = ctx["mp"]
    if mp.empty:
        return _req("Se requieren registros de MP con polímero.", ["Erni-Cassola et al. (2019)"])
    mat = mp.groupby(["matriz", "polimero"], observed=True)["abundancia"].sum().reset_index()
    mat["pct"] = mat.groupby("matriz")["abundancia"].transform(lambda s: s / s.sum() * 100)
    fig = figures.stacked(mat, "matriz", "pct", "polimero", title="Polímeros por matriz",
                          ytitle="% del total")
    return _ok(f"Dominio: {mat.sort_values('abundancia', ascending=False).iloc[0]['polimero']} en "
               f"{mat.sort_values('abundancia', ascending=False).iloc[0]['matriz']}. "
               "PE+PP suelen superar 50% (envases); PET alto indica textiles.",
               {"composicion": mat.round(1)}, {"apilado": fig},
               ["Erni-Cassola et al. (2019) Environ. Sci. Technol. 53:6044-6052"])


def ag_022(ctx):
    mp = ctx["mp"]
    if mp.empty:
        return _req("Se requieren registros de MP con morfología.", ["Cole et al. (2011)"])
    morph = mp.groupby(["zona", "morfologia"], observed=True)["abundancia"].sum().reset_index()
    morph["pct"] = morph.groupby("zona")["abundancia"].transform(lambda s: s / s.sum() * 100)
    fibra = morph[morph["morfologia"] == "fibra"].set_index("zona")["pct"]
    fig = figures.stacked(morph, "zona", "pct", "morfologia", title="Morfologías por zona", ytitle="%")
    return _ok(f"Participación de fibras por zona: {fibra.round(1).to_dict()}. "
               ">40% de fibras es firma típica de vertimientos domésticos/aguas grises.",
               {"morfologia": morph.round(1)}, {"apilado": fig}, ["Cole et al. (2011) Environ. Sci. Technol."])


def ag_023(ctx):
    mp = ctx["mp"]
    if mp.empty or "color" not in mp:
        return _req("Se requieren registros con color.", ["Enders et al. (2015)"])
    colors = mp.groupby("color", observed=True)["abundancia"].sum().reset_index()
    colors["pct"] = colors["abundancia"] / colors["abundancia"].sum() * 100
    calidos = colors.loc[colors["color"].isin(["rojo", "naranja", "amarillo", "rosa"]), "pct"].sum()
    fig = figures.pie(list(colors["color"]), list(colors["abundancia"]), title="Espectro de colores")
    return _ok(f"Colores cálidos (foto-degradación/aging): {calidos:.1f}%. "
               "Alta fracción transparente dificulta la ingestión visual por depredadores.",
               {"colores": colors.round(1)}, {"torta": fig}, ["Enders et al. (2015)"])


def ag_024(ctx):
    mp = ctx["mp"]
    mpa = mp[mp["matriz"] == "agua"]
    if mpa.empty or mpa["profundidad_m"].nunique() < 2:
        return _req("Se requieren perfiles verticales (≥2 profundidades).", ["Kowalski et al. (2016)"])
    per = mpa.groupby(["sample_id", "profundidad_m"], observed=True)["abundancia"].sum().reset_index()
    glm = multilevel.glm_counts(per.assign(profundidad=per["profundidad_m"]),
                                "abundancia", ["profundidad"], family="negativebinomial")
    fig = figures.scatter(per, x="profundidad_m", y="abundancia", trendline=True,
                          title="MP vs profundidad (GLM NB)", xtitle="m", ytitle="part/m³")
    return _ok(f"GLM NB: {glm.get('interpretacion', glm.get('error', ''))}. "
               "Descenso con profundidad esperado por flotabilidad y resuspensión costera.",
               {"perfil": per}, {"dispersion": fig}, ["Kowalski et al. (2016)"])


def ag_025(ctx):
    mp = ctx["mp"]
    if mp.empty:
        return _req("Se requieren campañas con temporada.", ["GESAMP (2019)"])
    per = mp.groupby(["matriz", "temporada", "sitio"], observed=True)["abundancia"].sum().reset_index()
    rows = []
    for matriz in per["matriz"].unique():
        sub = per[per["matriz"] == matriz]
        piv = sub.pivot_table(index="sitio", columns="temporada", values="abundancia").dropna()
        if piv.shape[1] >= 2 and len(piv) >= 3:
            w = inferential.wilcoxon_paired(piv.reset_index(), piv.columns[1], "sitio",
                                            piv.columns[0], piv.columns[1])
            ratio = piv[piv.columns[1]].mean() / piv[piv.columns[0]].mean()
            rows.append({"matriz": matriz, "ratio_lluvias_secas": round(float(ratio), 2),
                         "p": w.get("p"), "test": w.get("test", "")})
    fig = figures.box(mp.groupby(["sample_id", "temporada"], observed=True)["abundancia"].sum().reset_index(),
                      x="temporada", y="abundancia", title="MP por temporada", color="temporada")
    tab = pd.DataFrame(rows)
    return _ok(f"Comparación estacional por sitio (pareada): "
               f"{tab.to_dict('records') if len(tab) else 'insuficiente solapamiento sitio×temporada'}. "
               "Escorrentía pluvial eleva cargas en lluvias (esperado en sistemas urbanos).",
               {"pareado": tab}, {"cajas": fig}, ["GESAMP (2019) Guidelines for monitoring plastic litter"])


def ag_026(ctx):
    mp = ctx["mp"]
    mpa = mp[mp["matriz"] == "agua"].groupby(["sitio", "lat", "lon"], observed=True)["abundancia"].sum().reset_index()
    if len(mpa) < 5:
        return _req("Se requieren ≥5 estaciones georreferenciadas.", ["Getis & Ord (1992)"])
    gi = spatial.getis_ord_gi_star(mpa["lat"], mpa["lon"], mpa["abundancia"], k=6)
    gi["sitio"] = mpa["sitio"].values
    hot = gi[gi["tipo"] == "hotspot"]
    fig = figures.bars(gi.sort_values("gi_z", ascending=False), "sitio", "gi_z",
                       title="Hotspots Gi* de MP en agua", ytitle="z Gi*")
    return _ok(f"Hotspots 90%+: {', '.join(hot['sitio']) if len(hot) else 'ninguno significativo'} "
               f"({len(hot)} estaciones sobre P90).",
               {"getis_ord": gi.round(3)}, {"barras": fig},
               ["Getis & Ord (1992) Geogr. Anal. 24:189-206"])


def ag_027(ctx):
    mp = ctx["mp"]
    per = mp.groupby(["matriz", "sitio"], observed=True)["abundancia"].sum().unstack("matriz")
    if per.shape[1] < 2:
        return _req("Se requieren ≥2 matrices muestreadas.", ["Hartmann et al. (2019)"])
    rows = []
    for m1 in per.columns:
        for m2 in per.columns:
            if m1 >= m2:
                continue
            r = inferential.correlation_test(per[m1], per[m2], method="spearman")
            rows.append({"matriz_1": m1, "matriz_2": m2, "r": r.get("r"), "p": r.get("p")})
    tab = pd.DataFrame(rows)
    fig = figures.scatter(per.reset_index(), x=per.columns[0], y=per.columns[1],
                          title="Co-variación entre matrices", xtitle=per.columns[0], ytitle=per.columns[1],
                          text="sitio")
    return _ok(f"Correlaciones entre matrices: {tab.to_dict('records')}. "
               "Agua-sedimento altas indican intercambio activo (resuspensión).",
               {"correlaciones": tab.round(4)}, {"dispersion": fig}, ["Hartmann et al. (2019)"])


def ag_028(ctx):
    mp, bio = ctx["mp"], ctx["bio"]
    mpb = mp[mp["matriz"] == "biota"]
    if mpb.empty:
        return _req("Se requieren registros de MP en biota.", ["Lusher et al. (2013)"])
    from ..core.eco import biomarkers
    res = biomarkers.ingestion_metrics(mpb.assign(sitio=mpb["sitio"]))
    fig = figures.box(res["por_muestra"], x="sitio", y="carga_total",
                      title="MP ingeridos por sitio", ytitle="partículas/individuo")
    return _ok(f"Frecuencia de ingestión: {res['frecuencia_ingestion_pct']:.0f}%. "
               f"{res['interpretacion']}",
               {"resumen": res["resumen"]}, {"cajas": fig},
               ["Lusher et al. (2013) Mar. Pollut. Bull. 77:1-8"],
               ["Reportar SIEMPRE % de individuos con 0 partículas (cero inflado informativo)"])


def ag_029(ctx):
    mp = ctx["mp"]
    mpa = mp[mp["matriz"] == "agua"].groupby("sample_id")["abundancia"].sum()
    if mpa.empty:
        return _req("Se requieren concentraciones de MP en agua.")
    p = ctx["params"]
    volumen = p.get("volumen_m3", 3e7)
    mc = uncertainty.monte_carlo(lambda conc, vol: conc * vol,
                                 {"conc": ("lognormal", float(np.log(mpa.mean())), float(np.log(mpa).std())),
                                  "vol": ("uniform", volumen * .7, volumen * 1.3)}, n=10000)
    return _ok(f"Stock lagunar estimado (agua): {mc['salida']['p50']:.3g} partículas "
               f"({mc['salida']['p50'] / 1e12:.2f} billones; P97.5 {mc['salida']['p97.5']:.3g}) "
               f"con volumen supuesto {volumen:.0e} m³ (editable).",
               {"mc": pd.DataFrame([mc["salida"]]), "sensibilidad": mc["sensibilidad"]},
               {}, ["Elaboración OMNIPLAST (Monte Carlo)"],
               ["El stock depende del supuesto de mezcla completa: usar solo como orden de magnitud"])


def ag_030(ctx):
    mp, samples = ctx["mp"], ctx["samples"]
    if mp.empty:
        return _req("Se requieren datos de laboratorio.", REFS_QAQC)
    id_col = "sample_id" if "sample_id" in samples.columns else "id"
    qc_col = "qa_qc" if "qa_qc" in samples.columns else None
    if qc_col is None:
        return _req("La tabla de muestras no registra controles QA/QC.", REFS_QAQC)
    qc = samples.drop_duplicates(id_col).set_index(id_col)[qc_col].apply(
        lambda d: d if isinstance(d, dict) else {})
    with_blank = sum(1 for d in qc if d.get("blanco_proceso"))
    with_dup = sum(1 for d in qc if d.get("duplicado"))
    unknown_pct = (mp["polimero"] == "UNKNOWN").sum() / len(mp) * 100 if len(mp) else np.nan
    visual = mp["metodo_id"].str.contains("visual|Nile", case=False, na=False).mean() * 100
    issues = []
    if with_blank == 0:
        issues.append("Sin blancos de proceso registrados: riesgo de sobreestimación por contaminación aérea.")
    if unknown_pct > 20:
        issues.append(f"{unknown_pct:.0f}% de polímeros sin identificar: identifica con FTIR/Raman antes de PHI/riesgo.")
    if visual > 0:
        issues.append(f"{visual:.0f}% de partículas identificadas solo visualmente (no espectroscópico).")
    fig = figures.pie(["Blanco de proceso", "Duplicado de campo", "Sin control"],
                      [with_blank, with_dup, max(0, len(qc) - with_blank - with_dup)],
                      title="Cobertura de controles QA/QC")
    return _ok("Auditoría QA/QC: " + ("; ".join(issues) if issues else "sin alertas automáticas."),
               {"cobertura": pd.DataFrame({"control": ["blancos", "duplicados", "total_muestras"],
                                           "n": [with_blank, with_dup, len(qc)]})},
               {"torta": fig}, REFS_QAQC,
               ["Meta: <10% de UNKNOWN, blancos con <5% del total y recuperación 70-130% en spiked samples"])


def ag_031(ctx):
    mp = ctx["mp"]
    if mp.empty:
        return _req("Se requieren tamaños por registro.", ["Elaboración OMNIPLAST"])
    mp = mp.assign(tamano=(mp["tamano_min_um"] + mp["tamano_max_um"]) / 2)
    per_poly = mp.groupby("polimero", observed=True)["tamano"].median().reset_index()
    res = inferential.kruskal_wallis(mp, "tamano", "polimero")
    fig = figures.bars(per_poly.sort_values("tamano"), "polimero", "tamano",
                       title="Tamaño mediano por polímero", ytitle="µm")
    return _ok(f"{res.get('interpretacion', '')}. "
               "Filamentos PET/PA finos explican su biodisponibilidad elevada.",
               {"mediana_tamano": per_poly.round(0)}, {"barras": fig})


def ag_032(ctx):
    mp, env = ctx["mp"], ctx["env"]
    mpa = mp[mp["matriz"] == "agua"]
    if mpa.empty or env.empty:
        return _req("Se requieren MP en agua y DBO.", ["Browne et al. (2011)"])
    fibra = mpa.assign(es_fibra=mpa["morfologia"] == "fibra")
    ratio = fibra.groupby("sitio")["es_fibra"].mean()
    dbo = env[env["matriz"] == "agua"].groupby("sitio")["dbo_mg_l"].mean()
    df = pd.DataFrame({"ratio_fibra": ratio, "dbo": dbo}).dropna()
    if df.empty:
        return _req("Sin solapamiento sitio×DBO.")
    corr = inferential.correlation_test(df["ratio_fibra"], df["dbo"], method="spearman")
    fig = figures.scatter(df.reset_index(), x="dbo", y="ratio_fibra", trendline=True,
                          text="sitio", title="Fibras vs DBO por sitio", xtitle="DBO mg/L",
                          ytitle="proporción de fibras")
    return _ok(f"Correlación fibra-DBO: {corr.get('interpretacion', '')}. "
               "Las fibras (PET/PA) trazan vertimientos de aguas residuales domésticas.",
               {"por_sitio": df.round(3)}, {"dispersion": fig}, ["Browne et al. (2011) ES&T"])


HANDLERS_A = {f"AG-{i:03d}": f for i, f in enumerate([
    ag_001, ag_002, ag_003, ag_004, ag_005, ag_006, ag_007, ag_008, ag_009, ag_010,
    ag_011, ag_012, ag_013, ag_014, ag_015, ag_016, ag_017, ag_018,
    ag_019, ag_020, ag_021, ag_022, ag_023, ag_024, ag_025, ag_026, ag_027, ag_028,
    ag_029, ag_030, ag_031, ag_032], start=1)}
