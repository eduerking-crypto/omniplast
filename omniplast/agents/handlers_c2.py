"""Handlers G: plastisfera (AG-075..082) y H: biota/ecotoxicologia (AG-083..087)."""
import numpy as np
import pandas as pd

from ..core.eco import biomarkers
from ..core.genetics import community
from ..core.statistics import inferential
from ..reporting import figures
from ..utils import benjamini_hochberg
from .handlers_a import _ok, _req
from .handlers_c1 import _taxa


def ag_075(ctx):
    mb, mp = ctx["mb"], ctx["mp"]
    if mb is None or mb.empty:
        return _req("Se requieren conteos 16S.", ["Zettler et al. (2013)"])
    mps = mp[mp["matriz"] == "sedimento"].groupby("sample_id")["abundancia"].sum()
    mb2 = mb.merge(mps.rename("mp_total"), on="sample_id", how="inner")
    if mb2.empty or mb2["mp_total"].nunique() < 2:
        return _req("Sin solapamiento microbioma×MP para el contraste alto/bajo.")
    umbral = mb2["mp_total"].median()
    alta = mb2["mp_total"] > umbral
    taxa = _taxa(mb)
    rel = mb2[taxa].div(mb2[taxa].sum(axis=1).replace(0, 1), axis=0)
    rows = []
    for taxon in taxa:
        a, b = rel.loc[alta, taxon], rel.loc[~alta, taxon]
        if a.sum() + b.sum() == 0:
            continue
        from scipy import stats
        try:
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        except ValueError:
            continue
        pooled_sd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2) or 1e-12
        rows.append({"taxon": taxon, "rel_alta": a.mean(), "rel_baja": b.mean(),
                     "log2fc": np.log2((a.mean() + 1e-6) / (b.mean() + 1e-6)),
                     "d_cohen": (a.mean() - b.mean()) / pooled_sd, "p": p})
    tab = pd.DataFrame(rows)
    tab["p_adj"] = benjamini_hochberg(tab["p"].fillna(1).values)
    sig = tab[tab["p_adj"] < 0.05].sort_values("log2fc", ascending=False)
    fig = figures.bars(sig.head(10), "taxon", "log2fc",
                       title="Taxa enriquecidos en sitios con alta carga de MP",
                       ytitle="log2FC")
    return _ok(f"{len(sig)} taxa con enriquecimiento significativo (FDR<0.05). "
               + (f"Dominan: {', '.join(sig['taxon'].head(5))}" if len(sig) else
                  "Sin señal fuerte de enriquecimiento en este contraste."),
               {"enriquecimiento": sig.round(4)},
               {"barras": fig} if len(sig) else {},
               ["Zettler, Mincer & Amaral-Zettler (2013) Environ. Sci. Technol. 47:7137-7146"])


def ag_076(ctx):
    mb = ctx["mb"]
    if mb is None or mb.empty:
        return _req("Se requieren conteos 16S.", ["Shade & Handelsman (2012)"])
    r = community.core_microbiome(mb)
    return _ok(r["interpretacion"],
               {"core": r["core"].to_frame("prevalencia").reset_index(names="taxon"),
                "prevalencias": r["prevalencias"].to_frame("prev").reset_index(names="taxon")},
               {}, ["Shade & Handelsman (2012) FEMS Microbiol. Rev. 36:734-750"])


def ag_077(ctx):
    p = ctx["params"]
    if "tiempos" in p and "od" in p:
        t = np.asarray(p["tiempos"], float)
        od = np.asarray(p["od"], float)
        nota = "Series DEL USUARIO"
    else:
        t = np.array([0, 6, 12, 18, 24, 30, 36, 42, 48, 60, 72], float)
        lag = 12
        A, mu = 0.85, 0.08
        od = A * np.exp(-np.exp(mu * np.e / A * (lag - t) + 1))
        od = np.clip(od + np.random.default_rng(3).normal(0, .01, len(t)), 0, None)
        nota = "EJEMPLO SINTÉTICO (curva de crecimiento sobre película de plástico)"
    from scipy.optimize import curve_fit

    def gompertz(tt, A, mu, lam):
        return A * np.exp(-np.exp(mu * np.e / max(A, 1e-9) * (lam - tt) + 1))
    try:
        popt, _ = curve_fit(gompertz, t, od, p0=[np.max(od), .08, 10], maxfev=20000)
        A, mu, lam = popt
        pred = gompertz(t, *popt)
        r2 = 1 - ((od - pred) ** 2).sum() / ((od - od.mean()) ** 2).sum()
        fig = figures.scatter(pd.DataFrame({"tiempo_h": t, "od600": od}), x="tiempo_h",
                              y="od600", title=f"Gompertz ajustada — {nota}",
                              xtitle="horas", ytitle="OD600")
        fig.add_scatter(x=t, y=pred, mode="lines", name="Gompertz", line=dict(color="#FF6B6B"))
        return _ok(f"{nota}: A={A:.3f}, µ_max={mu:.4f} h⁻¹, λ={lam:.1f} h, R²={r2:.3f}. "
                   "µ máxima sobre plástico como sustrato ÚNICO es el indicador directo de "
                   "capacidad biodegradadora del consorcio.",
                   {"parametros": pd.DataFrame([{"A": A, "mu_max_h": mu, "lag_h": lam,
                                                 "r2": r2}])},
                   {"curva": fig},
                   ["Zwietering et al. (1990) Appl. Environ. Microbiol. 56:1875-1881"],
                   ["Complementar con ISO 14852 (mineralización acuosa) y pérdida de masa gravimétrica"])
    except Exception as e:
        return _req(f"Ajuste fallido: {e}. Verifique la serie (≥8 puntos, monótona creciente).",
                    ["Zwietering et al. (1990)"])


def ag_078(ctx):
    q = ctx["qpcr"]
    if q.empty:
        return _req("Se requieren datos qPCR del panel de genes.", ["Danso et al. (2019)"])
    weights = {"petase": 1.0, "mhetase": 0.8, "alkb": 0.6, "laccasa": 0.6, "phac": 0.3}
    rows = []
    for sitio, g in q.groupby("sitio", observed=True):
        score = 0.0
        parts = []
        for gene, w in weights.items():
            sub = g[g["target"] == gene]
            if len(sub) and sub["detectado"].any():
                copies = sub.loc[sub["detectado"], "valor"].median()
                contribution = w * float(np.log10(max(copies, 1)))
                score += contribution
                parts.append(f"{gene}={contribution:.1f}")
        rows.append({"sitio": sitio, "potencial_biodegradacion": round(score, 2),
                     "detalle": "; ".join(parts)})
    tab = pd.DataFrame(rows).sort_values("potencial_biodegradacion", ascending=False)
    fig = figures.bars(tab, "sitio", "potencial_biodegradacion", horizontal=True,
                       title="Potencial de biodegradación por sitio (índice génico ponderado)")
    return _ok(f"Sitio con mayor potencial: {tab.iloc[0]['sitio']} "
               f"({tab.iloc[0]['potencial_biodegradacion']}). "
               "Índice construido con pesos por papel del gen (petase>mhetase>alkB≈laccasa>phaC); "
               "es comparativo entre sitios, NO una tasa de degradación.",
               {"ranking": tab}, {"barras": fig},
               ["Danso, Chow & Streit (2019) AEM 85:e01095-19",
                "Yoshida et al. (2016) Science 351:1196-1203"])


def ag_079(ctx):
    q = ctx["qpcr"]
    if q.empty:
        return _req("Se requieren datos qPCR del panel.", ["Yoshida et al. (2016)"])
    genes = [g for g in ("petase", "mhetase", "laccasa", "alkb", "phac")
             if g in set(q["target"])]
    if len(genes) < 2:
        return _req("Se requieren ≥2 genes del panel en los datos.")
    wide = q.pivot_table(index="sample_id", columns="target", values="valor", aggfunc="mean")
    lg = wide[genes].apply(lambda s: np.log10(s.clip(lower=1)))
    rows = []
    for i, g1 in enumerate(genes):
        for g2 in genes[i + 1:]:
            r = inferential.correlation_test(lg[g1], lg[g2], method="spearman")
            rows.append({"gen_1": g1, "gen_2": g2, "r": r.get("r"), "p": r.get("p")})
    tab = pd.DataFrame(rows)
    mat = tab.pivot(index="gen_1", columns="gen_2", values="r")
    full = mat.combine_first(mat.T)
    arr = np.array(full.values, dtype=float, copy=True)
    np.fill_diagonal(arr, 1.0)
    full = pd.DataFrame(arr, index=full.index, columns=full.columns)
    fig = figures.heatmap(full, title="Correlación entre genes del panel (log₁₀ copias)")
    return _ok(f"Correlaciones gen-gen: {tab.dropna().sort_values('p').iloc[0]['gen_1']}×"
               f"{tab.dropna().sort_values('p').iloc[0]['gen_2']} más asociados "
               "(co-regulación sugiere consorcios degradadores con rutas complementarias).",
               {"correlaciones": tab.round(4)}, {"heatmap": fig},
               ["Yoshida et al. (2016) Science 351:1196-1203"])


def ag_080(ctx):
    mb = ctx["mb"]
    if mb is None or mb.empty:
        return _req("Se requieren conteos 16S.", ["Segata et al. (2011)"])
    r = community.lefse_effect(mb, mb["zona"])
    if "error" in r:
        return _req(r.get("error", "LEfSe no computable."), ["Segata et al. (2011)"])
    sig = r["tabla"][r["tabla"]["kw_p_adj"] < 0.05]
    return _ok(f"{r['interpretacion']}. "
               + (f"Biomarcadores: {', '.join(sig['taxon'].head(6))}" if len(sig) else ""),
               {"lefse": r["tabla"].head(15).round(4)}, {},
               ["Segata et al. (2011) Genome Biol. 12:R50 (lógica LEfSe adaptada)"])


def ag_081(ctx):
    q = ctx["qpcr"]
    targets = set(q.get("target", [])) if not q.empty else set()
    if "intI1" not in targets:
        return _req("Agregue el gen intI1 (integrasa de integrones clase 1) al panel qPCR "
                    "para cuantificar el co-marcador de resistencia antimicrobiana (AMR) "
                    "en la plastísfera. Protocolo: qPCR con cebadores publicados para intI1 "
                    "(p. ej. los de Gillings et al. 2015), incluyendo NTC y curva estándar.",
                    ["Gillings, M. R. et al. (2015) Integrons: past, present, and future",
                     "Revisión: plásticos como vector de AMR (Arias-Andres et al. 2018)"],
                    ["Contraste intI1 (copias/g) contra MP_total y genes de degradación con AG-067"])
    wide = q.pivot_table(index="sample_id", columns="target", values="valor", aggfunc="mean")
    mp = ctx["mp"]
    mps = mp[mp["matriz"] == "sedimento"].groupby("sample_id")["abundancia"].sum().rename("mp_total")
    d = wide.join(mps, how="inner").dropna(subset=["intI1"])
    if len(d) < 5:
        return _req("≥5 muestras con intI1 y MP medidos.")
    r = inferential.correlation_test(np.log10(d["intI1"].clip(lower=1)), d["mp_total"],
                                     method="spearman")
    return _ok(f"Correlación intI1×MP: {r.get('interpretacion', '')}. "
               "Co-variación positiva respalda la hipótesis del plástico como vector de AMR.",
               {"datos": d.round(2)}, {}, ["Gillings et al. (2015)"])


BIOFILM_TAXA = {"Vibrio": 1.0, "Pseudomonas": 1.0, "Bacillus": .9, "Marinobacter": .9,
                "Shewanella": .9, "Flavobacterium": .8, "Acinetobacter": .8,
                "Rhodococcus": .8, "Alcanivorax": .8, "Alteromonas": .8,
                "Tenacibaculum": .7, "Oleispira": .7, "Thalassospira": .7}


def ag_082(ctx):
    mb = ctx["mb"]
    if mb is None or mb.empty:
        return _req("Se requieren conteos 16S.", ["Oberbeckmann et al. (2015)"])
    taxa = _taxa(mb)
    rel = mb[taxa].div(mb[taxa].sum(axis=1).replace(0, 1), axis=0)
    score = sum(rel[t] * w for t, w in BIOFILM_TAXA.items() if t in rel.columns)
    df = pd.DataFrame({"sample_id": mb["sample_id"].values,
                       "sitio": mb["sitio"].values if "sitio" in mb else "",
                       "score_biofilm": score.values})
    resumen = df.groupby("sitio", observed=True)["score_biofilm"].mean().sort_values(ascending=False)
    fig = figures.bars(resumen.reset_index(), "sitio", "score_biofilm",
                       title="Score de formadores de biofilm sobre MP (ponderado)")
    return _ok(f"Score medio del sistema: {df['score_biofilm'].mean():.3f} (máx teórico 1.0). "
               "Taxa formadores de biofilm dominantes en la plastísfera suelen ser "
               "oprtunistas (Vibrio) — relevantes para sanidad de pesca y acuicultura.",
               {"score_por_sitio": resumen.round(3).reset_index()}, {"barras": fig},
               ["Oberbeckmann et al. (2015) Mar. Environ. Res. 120:220-232",
                "Zettler et al. (2013) ES&T 47:7137-7146"])


def _mpb(ctx):
    mp = ctx["mp"]
    return mp[mp["matriz"] == "biota"].copy()


def ag_083(ctx):
    mpb = _mpb(ctx)
    if mpb.empty:
        return _req("Se requieren registros de MP en biota.", ["Lusher et al. (2013)"])
    per = mpb.groupby(["sample_id", "subtipo"], observed=True)["abundancia"].sum().reset_index()
    especies = per["subtipo"].nunique()
    if especies < 2:
        return _req("Se requieren ≥2 especies para comparar.",
                    ["Lusher et al. (2013)"])
    e1, e2 = sorted(per["subtipo"].unique())[:2]
    a = per.loc[per["subtipo"] == e1, "abundancia"]
    b = per.loc[per["subtipo"] == e2, "abundancia"]
    from scipy import stats
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    fig = figures.box(per, x="subtipo", y="abundancia",
                      title="MP ingeridos por especie", ytitle="partículas/individuo")
    return _ok(f"{e1}: mediana {a.median():.0f} vs {e2}: mediana {b.median():.0f} "
               f"(Mann-Whitney p={p:.3g}). "
               "Los bivalvos filtran continuamente y suelen superar a peces en carga por individuo.",
               {"por_especie": per.groupby("subtipo")["abundancia"].describe().round(2)},
               {"cajas": fig}, ["Lusher et al. (2013) Mar. Pollut. Bull. 77:1-8"])


def ag_084(ctx):
    mpb = _mpb(ctx)
    if mpb.empty:
        return _req("Se requieren registros de MP en biota.", ["Lusher et al. (2013)"])
    per = mpb.groupby(["sample_id", "subtipo"], observed=True)["abundancia"].sum().reset_index()
    pct5 = (per["abundancia"] >= 5).mean() * 100
    fig = figures.box(per, x="subtipo", y="abundancia",
                      title="Intensidad de contaminación en tracto digestivo",
                      ytitle="partículas/individuo")
    return _ok(f"Carga media={per['abundancia'].mean():.1f} partículas/individuo; "
               f"{pct5:.0f}% de las muestras con ≥5 partículas (criterio de 'ingestión "
               "elevada' referencial en ecotoxicología de MP).",
               {"distribucion": per.round(1)}, {"cajas": fig},
               ["Lusher et al. (2013)"],
               ["Reportar siempre n individuos por pool y variabilidad individual (SD)"])


def ag_085(ctx):
    mpb = _mpb(ctx)
    if mpb.empty or mpb["subtipo"].nunique() < 2:
        return _req("Se requieren ≥2 especies con morfologías registradas.",
                    ["Browne et al. (2008)"])
    ct = pd.crosstab(mpb["subtipo"], mpb["morfologia"], values=mpb["abundancia"],
                     aggfunc="sum")
    ct = ct.fillna(0).round(0)
    if ct.shape[0] < 2 or ct.shape[1] < 2:
        return _req("Tabla de contingencia insuficiente (≥2×2).")
    r = inferential.chi_square_test(ct.values)
    fig = figures.heatmap(ct, title="Morfología ingerida × especie (partículas)")
    return _ok(f"{r['interpretacion']}. La selectividad de fibras vs fragmentos difiere "
               "entre grupos tróficos: filtradores retienen más fibras.",
               {"contingencia": ct}, {"heatmap": fig},
               ["Browne et al. (2008) Environ. Sci. Technol. 42:5026-5031"])


def ag_086(ctx):
    bio = ctx["bio"]
    markers = ["sod_u_mg_prot", "cat_u_mg_prot", "gst_u_mg_prot", "lpo_nmol_tba_mg"]
    sub = bio[bio["parametro"].isin(markers)] if not bio.empty else bio
    if sub.empty:
        return _req("Se requieren biomarcadores (SOD/CAT/GST/LPO) en biota.",
                    ["Regoli & Giuliani (2014)"])
    pivot = sub.pivot_table(index=["sample_id", "codigo", "zona", "subtipo"],
                            columns="parametro", values="valor", aggfunc="mean").reset_index()
    rows = []
    for m in markers:
        if m in pivot.columns and pivot["zona"].nunique() >= 2:
            kw = inferential.kruskal_wallis(pivot, m, "zona")
            rows.append({"biomarcador": m, "KW_H": kw.get("estadistico"),
                         "p": kw.get("p")})
    tab = pd.DataFrame(rows)
    fig = figures.box(pivot, x="zona", y="lpo_nmol_tba_mg",
                      title="LPO (peroxidación lipídica) por zona",
                      ytitle="nmol TBARS/mg") if "lpo_nmol_tba_mg" in pivot.columns else None
    figs = {"cajas": fig} if fig is not None else {}
    return _ok("LPO elevado con antioxidantes alterados configura el síndrome de estrés "
               "oxidativo compatible con exposición a MP + COVs adsorbidos. "
               + "; ".join(f"{r.biomarcador}: p={r.p:.3g}" for r in tab.itertuples()),
               {"pruebas": tab.round(4), "pivot": pivot.round(2)}, figs,
               ["Regoli & Giuliani (2014) Mar. Environ. Res. 96:62-77"],
               ["Interpretar por tejido: hepatopáncreas ≠ branquia"])


def ag_087(ctx):
    bio = ctx["bio"]
    markers = ["sod_u_mg_prot", "cat_u_mg_prot", "gst_u_mg_prot", "lpo_nmol_tba_mg"]
    sub = bio[bio["parametro"].isin(markers)] if not bio.empty else bio
    if sub.empty:
        return _req("Se requieren biomarcadores para el IBR.",
                    ["Beliaeff & Burgeot (2002)"])
    pivot = sub.pivot_table(index=["sample_id", "sitio", "zona"],
                            columns="parametro", values="valor", aggfunc="mean").reset_index()
    zonas = pivot["zona"].unique()
    if len(zonas) < 2:
        return _req("Se requieren ≥2 zonas (incluida una de referencia).",
                    ["Beliaeff & Burgeot (2002)"])
    ref_zona = zonas[0]
    ref = pivot[pivot["zona"] == ref_zona]
    rows, perfiles = [], {}
    for zona in zonas:
        if zona == ref_zona:
            continue
        tr = pivot[pivot["zona"] == zona]
        medias_t = {m: float(tr[m].mean()) for m in markers if m in tr.columns}
        medias_r = {m: float(ref[m].mean()) for m in markers if m in ref.columns}
        sds_r = {m: float(ref[m].std(ddof=1)) if len(ref) > 1 else 1.0
                 for m in markers if m in ref.columns}
        ibr = biomarkers.ibr_beliaeff(medias_t, medias_r, sds_r)
        if "ibr" not in ibr:
            continue
        rows.append({"zona": zona, "ibr": ibr["ibr"], "clase": ibr["interpretacion"].split("→")[-1].strip()})
        perfiles[zona] = ibr["perfil"]
    tab = pd.DataFrame(rows)
    if not tab.empty and perfiles:
        primera = list(perfiles.values())[0]
        fig = figures.radar(list(primera.keys()),
                            [abs(v) * 10 for v in primera.values()],
                            title="Perfil IBR (|z| ×10) — " + list(perfiles)[0])
        fig2 = figures.bars(tab, "zona", "ibr", title="IBR/2 por zona", ytitle="IBR/2")
        figs = {"radar": fig, "barras": fig2}
    else:
        figs = {}
    ref_nombre = ref_zona
    return _ok(f"IBR/2 por zona (referencia={ref_nombre}): "
               + "; ".join(f"{r.zona}: {r.ibr:.1f}" for r in tab.itertuples())
               if len(tab) else "IBR no computable con estos datos.",
               {"ibr": tab.round(2)}, figs,
               ["Beliaeff & Burgeot (2002) Biomarkers 7:15-34",
                "Lam & Gray (2003) Mar. Ecol. Prog. Ser. 255:183-196"])


HANDLERS_C2 = {f"AG-{i:03d}": f for i, f in enumerate([
    ag_075, ag_076, ag_077, ag_078, ag_079, ag_080, ag_081, ag_082,
    ag_083, ag_084, ag_085, ag_086, ag_087], start=75)}
