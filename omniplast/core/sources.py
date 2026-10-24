"""Apportionment de fuentes de MP por huella polimero-morfologia (enfoque probabilistico).

Huellas construidas a partir de la literatura de fuentes dominantes en humedales
costeros urbanos (revisar: Cole et al. 2011; Browne 2015; Bertoldi et al. 2020;
Gago et al. 2020). Los pesos son editables por el investigador.
"""
import numpy as np
import pandas as pd

SOURCE_FINGERPRINTS = {
    "Textiles y lavado doméstico": {("PET", "fibra"): .5, ("PA", "fibra"): .28,
                                    ("PP", "fibra"): .1, ("PE", "fibra"): .12},
    "Envases y embalaje urbano": {("PE", "fragmento"): .3, ("PP", "fragmento"): .24,
                                  ("PE", "film"): .2, ("PP", "film"): .14,
                                  ("PET", "fragmento"): .12},
    "Pesca y acuicultura": {("PE", "fragmento"): .2, ("PP", "fragmento"): .2,
                            ("PA", "fibra"): .2, ("PS", "espuma"): .15,
                            ("PVC", "fragmento"): .1, ("PE", "film"): .15},
    "Espumas de construcción y desechables": {("EPS", "espuma"): .6, ("PS", "espuma"): .3,
                                              ("PS", "fragmento"): .1},
    "Agricultura (mulch, invernadero)": {("PE", "film"): .55, ("PP", "film"): .2,
                                         ("PVC", "fragmento"): .1, ("PE", "fragmento"): .15},
    "Tráfico urbano (TRWP y polvo)": {("UNKNOWN", "fragmento"): .5, ("PP", "fragmento"): .2,
                                      ("PE", "fragmento"): .15, ("PA", "fibra"): .15},
    "Vertimientos de aguas residuales": {("PET", "fibra"): .2, ("PP", "fibra"): .12,
                                         ("PE", "film"): .18, ("PP", "fragmento"): .15,
                                         ("PE", "fragmento"): .15, ("PS", "fragmento"): .1,
                                         ("PA", "fibra"): .1},
    "Turismo y recreación": {("PP", "fragmento"): .25, ("PS", "espuma"): .2,
                             ("PE", "fragmento"): .2, ("PET", "fragmento"): .2,
                             ("PA", "fragmento"): .15},
}


def apportion_sources(mp_records, polymer_col="polimero", morph_col="morfologia",
                      count_col="abundancia", temperature=0.05):
    """Asignacion probabilistica: similitud coseno huella observada vs fuentes + softmax."""
    sub = mp_records.copy()
    sub[count_col] = pd.to_numeric(sub[count_col], errors="coerce")
    sub = sub.dropna(subset=[count_col])
    if sub.empty:
        return {"test": "Fuentes de MP", "error": "Sin registros"}
    observed = sub.groupby([polymer_col, morph_col], observed=True)[count_col].sum()
    total = observed.sum()
    obs_vec = observed / total if total > 0 else observed
    keys = set()
    for fp in SOURCE_FINGERPRINTS.values():
        keys.update(fp.keys())
    keys.update(obs_vec.index)
    keys = sorted(keys)
    obs = np.array([obs_vec.get(k, 0) for k in keys])
    scores = {}
    for src, fp in SOURCE_FINGERPRINTS.items():
        vec = np.array([fp.get(k, 0) for k in keys])
        cos = float(np.dot(obs, vec) / (np.linalg.norm(obs) * np.linalg.norm(vec) + 1e-12))
        scores[src] = cos
    if max(scores.values()) <= 0:
        return {"test": "Fuentes de MP", "error": "Composición no comparable con huellas"}
    exp = np.array(list(scores.values())) / temperature
    w = np.exp(exp - exp.max())
    probs = w / w.sum()
    order = np.argsort(probs)[::-1]
    table = pd.DataFrame({"fuente": [list(scores)[i] for i in order],
                          "aporte_pct": [float(probs[i] * 100) for i in order]}).round(1)
    norm_probs = probs / probs.sum()
    entropy = float(-(norm_probs * np.log(norm_probs + 1e-12)).sum())
    hmax = float(np.log(len(probs)))
    confidence = ("alta" if entropy / hmax < .55 else "media" if entropy / hmax < .8 else "baja")
    top = table.iloc[0]
    return {"test": "Apportionment de fuentes (huella polímero×morfología)",
            "tabla": table, "huella_observada": obs_vec.round(3).to_dict(),
            "entropia": entropy, "confianza": confidence,
            "interpretacion": f"Fuente dominante probable: {top['fuente']} ({top['aporte_pct']}%); "
                              f"confianza {confidence}. Las atribuciones son hipótesis contrastables "
                              "con campañas de caracterización local (peso editable)"}


def tire_wear_indicator(mp_records, unknown_fragment_share=None):
    """Proxy TRWP: proporcion de fragmentos negros UNKNOWN/PP de tamano pequeño."""
    sub = mp_records
    if sub.empty:
        return {"test": "TRWP", "error": "Sin datos"}
    small = pd.to_numeric(sub.get("tamano_max_um"), errors="coerce") <= 500
    dark = sub.get("color", pd.Series(index=sub.index, dtype=object)).astype(str).str.lower().isin(["negro", "gris"])
    frag = sub.get("morfologia", pd.Series(index=sub.index, dtype=object)).astype(str) == "fragmento"
    trwp_like = (small & frag).mean() * 100
    trwp_dark = (small & frag & dark).mean() * 100
    return {"test": "Indicador proxy de TRWP (desgaste de neumáticos)",
            "fragmentos_pequenos_pct": float(trwp_like), "fragmentos_negros_pct": float(trwp_dark),
            "interpretacion": (f"{trwp_dark:.1f}% de las partículas son fragmentos negros ≤500 µm, "
                              "fracción consistente con TRWP urbano; confirmar con FTIR/microscopía "
                              "de RRWP y análisis de negros de carbón")}
