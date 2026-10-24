"""Indices de calidad por microplasticos y del ecosistema.

PLI: Tomlinson et al. (1980). CF: Hakanson (1980). PHI (Polymer Hazard Index):
puntajes de riesgo quimico de Lithner et al. (2011) ES&T 45:4727-4735 combinados
con la abundancia porcentual (como en Turner & Holmes 2019 y estudios derivados).
WQI-NSF: Brown et al. (1970) Water Pollution Control Federation.
TSI: Carlson (1977) Limnol. Oceanogr. 22:361-369.
IICAM (indice propio de OMNIPLAST, propuesto): compone exposicion (carga MP),
peligro polimerico, estado fisicoquimico (WQI), respuesta biologica y toxicologica
en una escala 0-100 de calidad del ecosistema frente a MP, con pesos configurables
y cuantificacion de incertidumbre por Monte Carlo.
"""
import numpy as np
import pandas as pd

from ..config import (AMBI_CLASSES, BMWP_COL_CLASSES, CF_CLASSES, CHART_SEQUENCE,
                      HQ_CLASSES, MP_TIERS, MP_TIERS_NOTE, PHI_CLASSES, POLYMER_HAZARD,
                      QUALITY_CLASSES, TSI_CLASSES, WQI_CLASSES, WQI_WEIGHTS)
from ..utils import classify_value


def contamination_factor(mp_series, background):
    """CF = C_muestra / C_referencia (fondo). Requiere > 0."""
    bg = float(background) if background and background > 0 else np.nan
    if np.isnan(bg):
        return np.nan
    return float(mp_series) / bg


def pli_site(cf_series):
    """PLI = (prod CF_i)^(1/n); usa CFs de multiples polimeros/muestras del sitio."""
    cfs = np.asarray([c for c in cf_series if c is not None and c == c and c > 0], float)
    if len(cfs) == 0:
        return np.nan
    return float(np.exp(np.mean(np.log(cfs))))


def pli_classify(pli):
    if pli != pli:
        return "no computable"
    if pli < 1:
        return "No contaminado (PLI<1)"
    if pli == 1:
        return "Nivel basal (PLI=1)"
    return f"Deterioro de la calidad (PLI>1 → {pli:.2f})"


def cf_classify(cf):
    return classify_value(cf, CF_CLASSES)


def polymer_hazard_index(df_records, polymer_col="polimero", count_col="abundancia"):
    """PHI = sum(H_p x %abundancia_p) sobre polimeros con puntaje conocido."""
    sub = df_records.copy()
    sub[count_col] = pd.to_numeric(sub[count_col], errors="coerce")
    sub = sub.dropna(subset=[count_col])
    if sub.empty:
        return {"phi": np.nan, "clase": "sin datos", "desconocidos_pct": np.nan}
    total = sub[count_col].sum()
    scored = 0.0
    unknown = 0.0
    for poly, g in sub.groupby(polymer_col, observed=True):
        pct = g[count_col].sum() / total * 100
        hazard = POLYMER_HAZARD.get(str(poly).upper().strip(), None)
        if hazard is None:
            unknown += pct
        else:
            scored += hazard * pct
    phi = scored / 100 if total > 0 else np.nan
    clase = classify_value(phi, PHI_CLASSES)
    return {"phi": float(phi), "clase": clase, "desconocidos_pct": float(unknown),
            "nota": "Polimeros sin puntaje publico se excluyen del PHI y se reportan como % desconocido"}


def mp_tier(value, matrix="agua"):
    bands = MP_TIERS.get(matrix, MP_TIERS["agua"])
    return classify_value(float(value) if value == value else np.nan, bands)


def _wqi_rating(param, value):
    """Curvas de valoracion NSF aprox. por interpolacion de puntos publicados."""
    curves = {
        "ph": [(0, 0), (3, 10), (5, 40), (6.5, 70), (7, 85), (7.5, 92), (8, 88),
               (8.5, 75), (9, 50), (11, 10), (14, 0)],
        "temperatura_c": [(0, 20), (10, 35), (15, 55), (20, 75), (25, 88), (28, 90),
                          (30, 80), (32, 60), (35, 40), (40, 10)],
        "turbidez_ntu": [(0, 95), (5, 80), (10, 70), (25, 55), (50, 40), (100, 25),
                         (200, 10), (500, 5)],
        "dbo_mg_l": [(0, 100), (1, 95), (2, 90), (3, 85), (5, 75), (8, 60), (10, 50),
                     (15, 25), (20, 15), (30, 5), (50, 0)],
        "tss_mg_l": [(0, 95), (10, 85), (25, 70), (50, 55), (100, 40), (200, 20),
                     (400, 10), (1000, 5)],
        "nitrato_mg_l": [(0, 95), (2, 90), (5, 80), (10, 65), (20, 45), (50, 20), (100, 10)],
        "nitrito_mg_l": [(0, 95), (.01, 90), (.05, 80), (.1, 70), (.5, 50), (1, 35),
                         (2, 15), (5, 5)],
        "fosfato_mg_l": [(0, 95), (.02, 90), (.05, 80), (.1, 70), (.3, 50), (.6, 35),
                         (1, 15), (2, 5)],
        "coliformes_totales_nmp": [(1, 100), (10, 90), (100, 80), (200, 70), (500, 60),
                                  (1000, 45), (2000, 30), (5000, 15), (10000, 5), (1e6, 0)],
        "od_mg_l": [(0, 0), (2, 13), (4, 27), (6, 50), (8, 75), (10, 90), (12, 98),
                    (14, 100), (16, 95), (20, 85)],
    }
    curve = curves.get(param)
    if curve is None:
        return None
    xs = [c[0] for c in curve]
    ys = [c[1] for c in curve]
    return float(np.interp(value, xs, ys))


def wqi_nsf(sample_dict):
    """WQI NSF 9 parametros: 100 x sum(w_i x q_i) / sum(w_i)."""
    num, den = 0.0, 0.0
    ratings = {}
    for param, weight in WQI_WEIGHTS.items():
        val = sample_dict.get(param)
        if val is None or val != val:
            continue
        q = _wqi_rating(param, float(val))
        if q is None:
            continue
        num += weight * q
        den += weight
        ratings[param] = {"q": round(q, 1), "w": weight}
    if den == 0:
        return {"wqi": np.nan, "clase": "sin parámetros disponibles", "ratings": ratings}
    wqi = num / den
    return {"wqi": float(wqi), "clase": classify_value(wqi, WQI_CLASSES),
            "ratings": ratings, "n_parametros": len(ratings)}


def tsi_carlson(chla=None, tp_ug=None, sd_m=None):
    """TSI de Carlson; tp debe venir en ug/L."""
    tsi_vals = {}
    if chla is not None and chla > 0:
        tsi_vals["clorofila_a"] = 30.6 + 9.81 * np.log(chla)
    if tp_ug is not None and tp_ug > 0:
        tsi_vals["fósforo_total"] = 14.42 * np.log(tp_ug) + 4.15
    if sd_m is not None and sd_m > 0:
        tsi_vals["disco_secchi"] = 60 - 14.43 * np.log(sd_m)
    if not tsi_vals:
        return {"tsi": np.nan, "clase": "sin datos"}
    tsi = float(np.mean(list(tsi_vals.values())))
    return {"tsi": tsi, "por_variable": {k: float(v) for k, v in tsi_vals.items()},
            "clase": classify_value(tsi, TSI_CLASSES)}


def hazard_quotient(pec, pnec):
    if pnec is None or pnec == 0:
        return {"hq": np.nan, "clase": "PNEC no definido"}
    hq = pec / pnec
    return {"hq": float(hq), "clase": classify_value(hq, HQ_CLASSES),
            "interpretacion": ("HQ>1: riesgo potencial de efecto; HQ<1: riesgo aceptable "
                               "(marco PEC/PNEC, ecotoxicologia regulatoria EU R.10)")}


def _score_from_percentile(value, series, invert=True):
    """Score 0-100: percentil del valor en la distribucion observada."""
    s = pd.to_numeric(pd.Series(series), errors="coerce").dropna()
    if len(s) < 3 or value != value:
        return 50.0
    pct = float((s < value).mean() * 100)
    return pct if not invert else 100 - pct


def iicam(site_row, mp_all_sites, phi_site=None, wqi_site=None, bio_shannon=None,
          biomarker_score=None, weights=None, mc_runs=2000, seed=42):
    """Indice Integrado de Calidad Ambiental por Microplasticos (propuesto por OMNIPLAST).

    site_row: dict con 'mp_total' (abundancia total en la matriz), 'matriz'.
    mp_all_sites: serie con la abundancia total de todos los sitios (para percentiles).
    phi_site / wqi_site / bio_shannon / biomarker_score: subindices normalizados 0-100
    (si son None, se re-normalizan los pesos restantes).
    weights: dict con claves mp, phi, wqi, bio, tox (suman 1 idealmente).
    Retorna calidad 0-100 (100 = estado excelente) + clase + banda de incertidumbre MC.
    """
    rng = np.random.default_rng(seed)
    w = weights or {"mp": .30, "phi": .20, "wqi": .20, "bio": .15, "tox": .15}
    mp_val = site_row.get("mp_total")
    mp_score = (None if mp_val is None or mp_val != mp_val
                else _score_from_percentile(mp_val, mp_all_sites, invert=True))
    phi_score = float(phi_site) if phi_site is not None else None
    if phi_score is not None:
        phi_score = float(np.clip(100 - np.log10(max(phi_score, 1)) * 20, 0, 100))
    wqi_score = wqi_site
    bio_score = bio_shannon
    tox_score = biomarker_score
    subs = {"mp": mp_score, "phi": phi_score, "wqi": wqi_score, "bio": bio_score, "tox": tox_score}
    active = {k: v for k, v in subs.items() if v is not None and v == v}
    if not active:
        return {"iicam": np.nan, "clase": "sin información", "subindices": subs}
    w_active = {k: w[k] for k in active}
    wsum = sum(w_active.values())
    score = sum(active[k] * (w_active[k] / wsum) for k in active)
    mc = []
    for _ in range(mc_runs):
        perturbed = {k: np.clip(v + rng.normal(0, 5), 0, 100) for k, v in active.items()}
        mc.append(sum(perturbed[k] * (w_active[k] / wsum) for k in perturbed))
    lo, hi = np.percentile(mc, [2.5, 97.5])
    clase = classify_value(score, QUALITY_CLASSES)
    return {"iicam": float(score), "clase": clase,
            "ic95_mc": (float(lo), float(hi)),
            "subindices": {k: (float(v) if v is not None else None) for k, v in subs.items()},
            "pesos_usados": w_active,
            "componentes_activos": list(active.keys()),
            "nota": "IICAM: índice compuesto propuesto (documentado en docs/METODOLOGIA.md); "
                    "validar con panel de expertos antes de uso regulatorio"}


def ambi(species_table):
    """species_table: DataFrame con columnas [especie, abundancia, grupo_ecologico (I-V)]."""
    df = species_table.copy()
    if "grupo_ecologico" not in df.columns:
        return {"ambi": np.nan, "clase": "sin grupos ecológicos asignados (requiere lista AZTI)"}
    df["abundancia"] = pd.to_numeric(df["abundancia"], errors="coerce").fillna(0)
    total = df["abundancia"].sum()
    if total <= 0:
        return {"ambi": np.nan, "clase": "sin individuos"}
    pcts = {}
    for grp in ["I", "II", "III", "IV", "V"]:
        pcts[grp] = df.loc[df["grupo_ecologico"] == grp, "abundancia"].sum() / total * 100
    coef = {"I": 0, "II": 1.5, "III": 3, "IV": 4.5, "V": 6}
    ambi_val = sum(coef[g] * pcts[g] for g in pcts) / 10
    return {"ambi": float(ambi_val),
            "clase": classify_value(ambi_val, AMBI_CLASSES),
            "composicion_pct": {k: round(v, 2) for k, v in pcts.items()},
            "formula": "AMBI = (0·%I + 1.5·%II + 3·%III + 4.5·%IV + 6·%V)/10 (Borja et al. 2000)"}


def bmwp_col(families_scores):
    """families_scores: dict familia -> puntaje de sensibilidad (1-10) presente (abundancia>0)."""
    total = sum(families_scores.values())
    n_fam = len(families_scores)
    aspt = total / n_fam if n_fam else np.nan
    clase = classify_value(total, BMWP_COL_CLASSES)
    return {"bmwp_col": float(total), "aspt": float(aspt), "n_familias": n_fam,
            "clase": clase,
            "nota": "BMWP-Col: adaptación de Roldán (2003) para Colombia"}


def shannon_index(counts):
    counts = np.asarray(counts, float)
    counts = counts[counts > 0]
    if len(counts) == 0:
        return np.nan
    p = counts / counts.sum()
    return float(-(p * np.log(p)).sum())


def richness_hill(counts, order=0):
    """Numero de Hill de orden q (q=0 riqueza, q=1 Shannon exp, q=2 Simpson inverso)."""
    counts = np.asarray(counts, float)
    counts = counts[counts > 0]
    if len(counts) == 0:
        return np.nan
    p = counts / counts.sum()
    if order == 0:
        return float(len(counts))
    if abs(order - 1) < 1e-9:
        return float(np.exp(-(p * np.log(p)).sum()))
    return float((p ** order).sum() ** (1 / (1 - order)))
