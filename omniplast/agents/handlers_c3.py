"""Handlers H-biota restantes, I-geo/teledeteccion y J-riesgo (AG-088..100)."""
import numpy as np
import pandas as pd

from ..config import DEMO_HULL
from ..core import indices, sources
from ..core.eco import biomarkers, onehealth, toxicology
from ..core.connectors import connectors
from ..core.geo import remote_sensing, spatial, trajectory
from ..core.statistics import inferential
from ..reporting import figures
from .handlers_a import _ok, _req


def ag_088(ctx):
    bio = ctx["bio"]
    sub = bio[bio["parametro"].isin(["factor_condicion", "hsi_pct", "gsi_pct"])] if not bio.empty else bio
    if sub.empty:
        return _req("Se requieren CF/HSI/GSI de peces.", ["Bagenal (1978)"])
    pivot = sub.pivot_table(index=["sample_id", "codigo", "zona", "subtipo"],
                            columns="parametro", values="valor", aggfunc="mean").reset_index()
    rows = []
    for m in ("factor_condicion", "hsi_pct", "gsi_pct"):
        if m in pivot.columns and pivot["zona"].nunique() >= 2:
            kw = inferential.kruskal_wallis(pivot, m, "zona")
            rows.append({"indice": m, "KW_H": kw.get("estadistico"), "p": kw.get("p")})
    tab = pd.DataFrame(rows)
    fig = figures.box(pivot, x="zona", y="factor_condicion",
                      title="Factor de condición por zona") if "factor_condicion" in pivot.columns else None
    figs = {"cajas": fig} if fig is not None else {}
    return _ok("Índices de condición energética: "
               + "; ".join(f"{r.indice}: p={r.p:.3g}" for r in tab.itertuples())
               + ". CF reducido con HSI elevado sugiere costo energético de desintoxicación.",
               {"pruebas": tab.round(4), "pivot": pivot.round(3)}, figs,
               ["Bagenal (1978) Methods for Assessment of Fish Production"])


def ag_089(ctx):
    p = ctx["params"]
    if "lesiones" in p and len(p["lesiones"]) >= 3:
        df = pd.DataFrame(p["lesiones"])
        nota = "Datos DEL USUARIO"
    else:
        df = pd.DataFrame([
            {"lesion": "Vacuolización lipídica hepatocitaria", "prevalencia_pct": 70,
             "severidad": 2, "factor_importancia": 3},
            {"lesion": "Hiperplasia epitelial branquial", "prevalencia_pct": 45,
             "severidad": 1, "factor_importancia": 2},
            {"lesion": "Infiltración leucocitaria focal", "prevalencia_pct": 60,
             "severidad": 2, "factor_importancia": 2},
            {"lesion": "Necrosis focal", "prevalencia_pct": 15, "severidad": 4,
             "factor_importancia": 4},
        ])
        nota = "EJEMPLO SINTÉTICO (estructura requerida)"
    r = biomarkers.histopathology_index(df)
    if "error" in r:
        return _req(r["error"], ["Bernet et al. (1999)"])
    return _ok(f"{nota}: {r['interpretacion']}",
               {"indice": r["tabla"].round(2)}, {},
               ["Bernet et al. (1999) J. Fish Dis. 22:23-34"],
               ["Claves: 0=sin alteración a 5=masiva; importancia 1=cicatricial a 4=progresiva"])


def ag_090(ctx):
    mp = ctx["mp"]
    if mp.empty:
        return _req("Se requieren registros de MP.", ["Roldán (2003); Borja et al. (2000)"])
    mpa = mp[mp["matriz"] == "agua"]
    shannon_rows = []
    for sitio, g in mpa.groupby("sitio", observed=True):
        morfotipos = g.groupby("morfologia", observed=True)["abundancia"].sum().values
        shannon_rows.append({"sitio": sitio,
                             "shannon_morfotipos": indices.shannon_index(morfotipos)})
    tab = pd.DataFrame(shannon_rows).sort_values("shannon_morfotipos")
    fig = figures.bars(tab, "sitio", "shannon_morfotipos",
                       title="Diversidad de morfotipos de MP (Shannon)", ytitle="H'")
    return _ok("La diversidad de morfotipos proxy la mezcla de fuentes: H' alta = "
               "múltiples fuentes (urbano+pesca+turismo); H' baja = fuente dominante única. "
               "PARA AMBI/BMWP-Col: requiere conteos de macroinvertebrados bentónicos con "
               "grupos ecológicos AZTI (AMBI) o familias y puntajes (BMWP-Col, Roldán 2003: "
               ">120 muy buena, 101-120 buena, 61-100 aceptable, 36-60 dudosa, 16-35 crítica).",
               {"shannon": tab.round(3)}, {"barras": fig},
               ["Roldán (2003) Análisis y teoría de las comunidades dulceacuícolas",
                "Borja et al. (2000) Mar. Pollut. Bull. 40:1100-1108"])


def _bands(ctx):
    p = ctx["params"]
    if "bands" in p and isinstance(p.get("bands"), dict):
        return p["bands"], p.get("lagoon_mask"), "RASTER DEL USUARIO"
    escena = remote_sensing.synthetic_lagoon_bands()
    return escena["bands"], escena["lagoon_mask"], "ESCENA SINTÉTICA DE DEMOSTRACIÓN"


def ag_091(ctx):
    bands, mask, nota = _bands(ctx)
    try:
        ndwi = remote_sensing.compute_index(bands, "NDWI")
        ndti = remote_sensing.compute_index(bands, "NDTI")
        ndci = remote_sensing.compute_index(bands, "NDCI")
        mndwi = remote_sensing.compute_index(bands, "MNDWI")
    except KeyError as e:
        return _req(f"El índice requiere la banda {e}. Use bandas Sentinel-2 (B3,B4,B5,B8,B11).")
    vals = lambda a: float(np.nanmean(a[mask])) if mask is not None else float(np.nanmean(a))
    tab = pd.DataFrame({"indice": ["NDWI", "MNDWI", "NDTI", "NDCI"],
                        "media": [vals(ndwi), vals(mndwi), vals(ndti), vals(ndci)],
                        "descripcion": ["agua (McFeeters)", "agua SWIR (Xu)",
                                        "turbidez (Lacaux)", "clorofila (Mishra)"]})
    fig = figures.heatmap(ndwi, title=f"NDWI — {nota}", colorscale="Teal")
    fig2 = figures.heatmap(ndti, title=f"NDTI (turbidez) — {nota}", colorscale="Oranges")
    return _ok(f"{nota}: NDWI medio={vals(ndwi):.3f} (>0.2 = agua), NDTI medio={vals(ndti):.3f}. "
               "La pluma de turbidez derivada del canal urbano es visible en NDTI.",
               {"indices": tab.round(3)}, {"ndwi": fig, "ndti": fig2},
               ["McFeeters (1996) IJRS 17:1425-1432", "Xu (2006) IJRS 27:3025-3033",
                "Lacaux et al. (2007) IJRS 28:2877-2905"])


def ag_092(ctx):
    bands, mask, nota = _bands(ctx)
    try:
        fdi = remote_sensing.compute_index(bands, "FDI")
    except KeyError as e:
        return _req(f"FDI requiere banda {e} (B4, B8, B11).")
    mask_fdi, thr = remote_sensing.classify_floating_plastics(fdi, water_mask=mask)
    agua_px = int(mask.sum()) if mask is not None else fdi.size
    cand_pct = float(mask_fdi.sum()) / agua_px * 100 if agua_px else 0
    fig = figures.heatmap(fdi, title=f"FDI (plásticos flotantes candidatos) — {nota}",
                          colorscale="Viridis")
    return _ok(f"{nota}: umbral FDI=P98.5={thr:.4f}; {cand_pct:.2f}% del área de agua con "
               "firma de escombros flotantes. VALIDACIÓN OBLIGATORIA en campo: FDI detecta "
               "'floating debris' genérico (espumas, madera, Sargassum), no polímeros puros "
               "(Biermann et al. 2020: precisión 86% en escenarios controlados).",
               {"fdi_stats": pd.DataFrame([{"umbral": thr, "candidatos_pct": cand_pct,
                                            "pixeles_agua": agua_px}])},
               {"fdi": fig},
               ["Biermann et al. (2020) Sci. Rep. 10:19836", "Hu (2009) RSE 113:398-409"])


def ag_093(ctx):
    mp = ctx["mp"]
    mpa = mp[mp["matriz"] == "agua"].groupby(["sitio", "lat", "lon"], observed=True)["abundancia"].mean().reset_index()
    if len(mpa) < 5:
        return _req("≥5 estaciones georreferenciadas.", ["Isaaks & Srivastava (1989)"])
    grid_lon, grid_lat = spatial.build_grid(DEMO_HULL, n=40)
    grid = spatial.idw_interpolate(mpa["lat"], mpa["lon"], mpa["abundancia"],
                                   grid_lon, grid_lat, power=2)
    kr, fit = spatial.ordinary_kriging(mpa["lat"], mpa["lon"], mpa["abundancia"],
                                       grid_lon, grid_lat)
    fig = figures.heatmap(grid, title="Superficie IDW de MP en agua (part/m³)")
    variograma = spatial.variogram_empirical(mpa["lat"], mpa["lon"], mpa["abundancia"])
    kriging_nota = ""
    if fit:
        kriging_nota = (f" Variograma: nugget={fit['nugget']:.3g}, meseta={fit['meseta']:.3g}, "
                        f"alcance={fit['alcance_km']:.2f} km (R²={fit['r2']:.2f}).")
    return _ok(f"Mapa de concentración interpolada (IDW, potencia 2) sobre {len(mpa)} estaciones."
               + kriging_nota,
               {"estaciones": mpa.round(2)}, {"superficie": fig},
               ["Isaaks & Srivastava (1989) An Introduction to Applied Geostatistics"],
               ["Para publicación, reportar validación cruzada leave-one-out del kriging"])


def ag_094(ctx):
    mp = ctx["mp"]
    mpa = mp[mp["matriz"] == "agua"].groupby(["sitio", "lat", "lon"], observed=True)["abundancia"].mean().reset_index()
    if len(mpa) < 5:
        return _req("≥5 estaciones.", ["Moran (1950)"])
    mi = spatial.morans_i(mpa["lat"], mpa["lon"], mpa["abundancia"], k=8)
    lisa = spatial.local_moran(mpa["lat"], mpa["lon"], mpa["abundancia"], k=8)
    lisa.insert(0, "sitio", mpa["sitio"].values)
    hh = lisa[lisa["clase"].str.startswith("HH")]
    return _ok(f"{mi['interpretacion']}. Clusters HH (hotspots locales): "
               f"{', '.join(hh['sitio']) if len(hh) else 'ninguno'}. "
               "La autocorrelación espacial viola el supuesto de independencia: usar modelos "
               "mixtos/GEE (AG-046/047) en vez de pruebas clásicas.",
               {"moran_global": pd.DataFrame([mi]), "lisa": lisa.round(3)},
               {}, ["Moran (1950) Biometrika 37:17-23", "Anselin (1995) Geogr. Anal. 27:93-115"])


def ag_095(ctx):
    mp = ctx["mp"]
    mpa = mp[mp["matriz"] == "agua"].groupby(["sitio", "lat", "lon"], observed=True)["abundancia"].mean().reset_index()
    origen = ctx["params"].get("origen")
    if origen is None:
        if mpa.empty:
            return _req("Sin estaciones: fije 'origen' en parámetros.")
        src = mpa.sort_values("abundancia", ascending=False).iloc[0]
        origen = (float(src["lat"]), float(src["lon"]))
        nota = f"Origen automático = sitio más contaminado ({src['sitio']})"
    else:
        nota = "Origen fijado por el usuario"
    r = trajectory.simulate_drift(origen[0], origen[1], hours=24, n_particles=1500,
                                  bbox=DEMO_HULL,
                                  wind_speed_ms=ctx["params"].get("viento_ms", 3.5),
                                  wind_direction_deg=ctx["params"].get("viento_dir", 45),
                                  current_speed_ms=ctx["params"].get("corriente_ms", 0.12),
                                  current_direction_deg=ctx["params"].get("corriente_dir", 100),
                                  diameter_mm=ctx["params"].get("diametro_mm", 0.5),
                                  polymer_density=ctx["params"].get("densidad", 960))
    pos = r["posiciones_finales"]
    fig = figures.scatter(pos, x="lon", y="lat", color="estado",
                          title=f"Dispersión 24 h — {nota}", xtitle="longitud", ytitle="latitud")
    fig2 = figures.pie(list(r["destinos"].keys()), list(r["destinos"].values()),
                       title="Destino de partículas simuladas")
    return _ok(f"{nota}. {r['interpretacion']}",
               {"destinos": pd.DataFrame([r["destinos"]]),
                "desplazamiento": pd.DataFrame([r["desplazamiento_m"]])},
               {"trayectorias": fig, "destinos": fig2},
               ["Kowalski et al. (2016) ES&T 50:6825-6834",
                "Waldschläger & Schüttrumpf (2019) Water Res. 162:369-382"])


def ag_096(ctx):
    res = remote_sensing.fetch_sentinel2_items(
        DEMO_HULL, date_range="2024-01-01/2025-12-31", cloud_max=20, limit=8,
        api_key=ctx["params"].get("pc_key"))
    if res["ok"]:
        items = pd.DataFrame(res["items"])
        return _ok(f"{len(items)} escenas Sentinel-2 L2A disponibles en el área "
                   f"(2024-2025, nubes<20%). La descarga de píxeles requiere token de "
                   "Planetary Computer (registro gratuito) o uso directo desde Python "
                   "(pystac-client + planetary-computer).",
                   {"escenas": items}, {},
                   ["Microsoft Planetary Computer (https://planetarycomputer.microsoft.com)"])
    return _req(f"Sin acceso al catálogo STAC: {res['error']}. {res.get('hint', '')}",
                ["Microsoft Planetary Computer"],
                ["Use AG-091/092 con la escena sintética mientras configura el token"])


def ag_097(ctx):
    mp = ctx["mp"]
    mpa = mp[mp["matriz"] == "agua"]
    if mpa.empty:
        return _req("Se requieren MP en agua con polímero identificado.",
                    ["Frias & Nash (2019)"])
    pec = mpa.groupby("polimero", observed=True)["abundancia"].mean() / 1000.0
    r = toxicology.hazard_quotient_table({k: float(v) for k, v in pec.items()})
    tab = r["tabla"]
    fig = figures.bars(tab.dropna(subset=["hq"]), "polimero", "hq",
                       title="Cociente de riesgo HQ por polímero (PEC/PNEC)", ytitle="HQ")
    return _ok(f"{r['interpretacion']}. {r['nota']}",
               {"hq": tab.round(3), "pnec_usados": pd.DataFrame(
                   {"polimero": list(r["pnec_default"].keys()),
                    "pnec_p_L": list(r["pnec_default"].values())})},
               {"barras": fig},
               ["EU Technical Guidance Document R.10 (2003)",
                "Frias & Nash (2019) Mar. Pollut. Bull. 138:316-322"],
               ["PNECs referenciales: NO son valores regulatorios oficiales para MP"])


def ag_098(ctx):
    mp = ctx["mp"]
    mpb = mp[mp["matriz"] == "biota"]
    if mpb.empty:
        return _req("Se requieren MP en biota (partículas por individuo).",
                    ["WHO (2019)"])
    per = mpb.groupby(["sample_id", "subtipo"], observed=True)["abundancia"].sum().reset_index()
    masa_g = ctx["params"].get("masa_tejido_g", {"ostra": 15.0, "pez": 120.0})
    conc = []
    for row in per.itertuples():
        clave = "ostra" if "ostra" in str(row.subtipo).lower() else "pez"
        g = masa_g.get(clave, 100.0)
        conc.append(max(0.0, row.abundancia / max(g, 1e-9)))
    r = onehealth.exposure_seafood(np.array(conc))
    if "error" in r:
        return _req(r["error"], ["WHO (2019)"])
    fig = figures.bars(pd.DataFrame({"percentil": list(r["p_kg_bw_dia"].keys()),
                                     "valor": list(r["p_kg_bw_dia"].values())}),
                       "percentil", "valor", title="Exposición (partículas/kg pc/día)",
                       ytitle="p/kg/día")
    return _ok(f"{r['interpretacion']}",
               {"exposicion": pd.DataFrame([r["p_kg_bw_dia"]]),
                "sensibilidad": r["sensibilidad"].round(3)},
               {"barras": fig},
               ["WHO (2019) Microplastics in drinking-water",
                "Nota: sin umbral toxicológico regulatorio; comparar escenarios entre sí"],
               ["Supuestos: porción de {masa_g} g y frecuencia semanal editable en parámetros"])


def ag_099(ctx):
    p = ctx["params"]
    if "encuesta" in p and len(p["encuesta"]) >= 5:
        df = pd.DataFrame(p["encuesta"])
        nota = "Encuesta DEL USUARIO"
    else:
        rng = np.random.default_rng(11)
        items = [f"p{i}" for i in range(1, 6)]
        data = {"p1": rng.integers(1, 6, 24), "p2": rng.integers(1, 6, 24),
                "p3": rng.integers(1, 6, 24), "p4": rng.integers(1, 6, 24),
                "p5": rng.integers(1, 6, 24),
                "estrato": rng.choice(["pescador", "turista", "residente"], 24)}
        df = pd.DataFrame(data)
        nota = "EJEMPLO SINTÉTICO (escala Likert 1-5, 24 encuestados)"
    likert = [c for c in df.columns if c.startswith("p") and df[c].dtype.kind in "iu"]
    if len(likert) < 3:
        return _req("Se requieren ≥3 ítems Likert (columnas p1..pn, escala 1-5).",
                    ["Likert (1932)"])
    r = onehealth.social_survey_analysis(df, likert, group_col="estrato" if "estrato" in df.columns else None)
    return _ok(f"{nota}: {r['interpretacion']}",
               {"resumen": r["resumen"], "por_grupo": r.get("por_grupo", pd.DataFrame())},
               {}, ["Likert (1932) Arch. Psychol. 140:1-55"],
               ["Sugerencia: incluya preguntas sobre percepción de riesgo, disposición a pagar "
                "por productos sin plástico y conocimiento del problema"])


def ag_100(ctx):
    mp = ctx["mp"]
    env = ctx["env"]
    q = ctx["qpcr"]
    hallazgos = []

    mps = mp[mp["matriz"] == "sedimento"].groupby("sitio")["abundancia"].sum()
    if not mps.empty:
        bg = mps.min()
        pli = indices.pli_site([indices.contamination_factor(mps.max(), bg)])
        hallazgos.append({"prioridad": "Alta" if pli and pli > 1 else "Media",
                          "hallazgo": f"PLI máx={pli:.2f}" if pli == pli else "PLI n/d",
                          "accion": "Programa de reducción en fuentes terrestres con "
                                    "focalización en hotspots (canales urbanos)",
                          "evidencia": f"{mps.max():.0f} part/kg vs fondo {bg:.0f}"})
    phi = indices.polymer_hazard_index(mp)
    if phi["phi"] == phi["phi"]:
        hallazgos.append({"prioridad": "Alta" if phi["phi"] >= 100 else "Media",
                          "hallazgo": f"PHI={phi['phi']:.0f} ({phi['clase']})",
                          "accion": "Sustitución regulatoria de polímeros de alto peligro "
                                    "(PVC/ABS/PA) en usos de un solo uso",
                          "evidencia": f"{phi['desconocidos_pct']:.0f}% sin identificar"})
    agua = env[env["matriz"] == "agua"] if not env.empty else env
    if not agua.empty:
        wqi_vals = []
        for _, g in agua.groupby("sitio"):
            wqi_vals.append(indices.wqi_nsf(g.iloc[0].to_dict())["wqi"])
        wqi_vals = [v for v in wqi_vals if v == v]
        if wqi_vals:
            wqi_media = float(np.mean(wqi_vals))
            hallazgos.append({"prioridad": "Alta" if wqi_media < 50 else "Media",
                              "hallazgo": f"WQI medio={wqi_media:.0f}/100",
                              "accion": "Saneamiento del sistema hidráulico (PTAR y "
                                        "desviación de vertimientos crudos)",
                              "evidencia": f"mínimo {min(wqi_vals):.0f} en estación peor"})
        if "e_coli_nmp" in agua:
            exc = (agua["e_coli_nmp"] > 410).mean() * 100
            hallazgos.append({"prioridad": "Alta" if exc > 30 else "Media",
                              "hallazgo": f"{exc:.0f}% muestras > STV E. coli (410 NMP)",
                              "accion": "Alertas sanitarias para pesca artesanal y contacto "
                                        "primario; vigilancia epidemiológica",
                              "evidencia": f"mediana {agua['e_coli_nmp'].median():.0f} NMP/100mL"})
    if not mp.empty:
        mpa = mp[mp["matriz"] == "agua"]
        fibra_pct = (mpa["morfologia"] == "fibra").sum() / len(mpa) * 100
        hallazgos.append({"prioridad": "Media" if fibra_pct > 40 else "Baja",
                          "hallazgo": f"{fibra_pct:.0f}% de registros son fibras",
                          "accion": "Filtros de microfibra en lavanderías domésticas/industriales "
                                    "y rejas finas en PTAR",
                          "evidencia": "morfología dominante"})
        carga_mediana = mpa.groupby("sample_id")["abundancia"].sum().median()
        hallazgos.append({"prioridad": "Media" if carga_mediana > 1000 else "Baja",
                          "hallazgo": f"Carga mediana agua={carga_mediana:.0f} part/m³",
                          "accion": "Barreras flotantes y limpieza focalizada en desembocaduras "
                                    "de canales con monitoreo Gi* trimestral",
                          "evidencia": indices.mp_tier(carga_mediana, "agua")})
    if not q.empty and "petase" in set(q["target"]):
        det = q[q["target"] == "petase"]["detectado"].mean() * 100
        hallazgos.append({"prioridad": "Media",
                          "hallazgo": f"petase detectado en {det:.0f}% de sedimentos",
                          "accion": "Investigación aplicada en bio-remediación: enriquecer "
                                    "consorcios degradadores nativos (línea doctoral)",
                          "evidencia": "qPCR panel de degradación"})
    orden = {"Alta": 0, "Media": 1, "Baja": 2}
    tab = pd.DataFrame(hallazgos)
    tab = tab.iloc[tab["prioridad"].map(orden).sort_values().index].reset_index(drop=True)
    return _ok(f"Síntesis trazable: {len(tab)} recomendaciones priorizadas por umbrales "
              "incumplidos. Cada acción es auditable contra su indicador numérico.",
               {"recomendaciones": tab}, {},
               ["Umbrales citados en config.py: Hakanson 1980; Lithner 2011; Brown 1970; "
                "US EPA 2012; Kruskal/Getis en la literatura interna del proyecto"],
               ["Revisar semestralmente: los umbrales referenciales evolucionan con la literatura"])


HANDLERS_C3 = {f"AG-{i:03d}": f for i, f in enumerate([
    ag_088, ag_089, ag_090, ag_091, ag_092, ag_093, ag_094, ag_095, ag_096,
    ag_097, ag_098, ag_099, ag_100], start=88)}
