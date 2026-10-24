"""Evaluacion de riesgo ecotoxicologico: HQ (PEC/PNEC), SSD, PAF."""
import numpy as np
import pandas as pd
from scipy import stats

from ...config import HQ_CLASSES


def hazard_quotient_table(conc_by_polymer, pnec_by_polymer=None):
    """HQ por polímero: PEC mediana vs PNEC (valores de literatura, editables)."""
    default_pnec_water = {"PE": 1.9, "PP": 1.1, "PS": 8.0, "PET": 5.0, "PVC": 50.0,
                          "PA": 3.0, "PU": 100.0}
    pnec = pnec_by_polymer or default_pnec_water
    rows = []
    for poly, conc in conc_by_polymer.items():
        ref = pnec.get(poly)
        if ref is None:
            rows.append({"polimero": poly, "pec": float(conc), "pnec": np.nan,
                         "hq": np.nan, "clase": "PNEC no disponible (definir por literatura)"})
            continue
        hq = float(conc) / float(ref)
        clase = [c for c, lo, hi in HQ_CLASSES if lo <= hq < hi][-1]
        rows.append({"polimero": poly, "pec": float(conc), "pnec": float(ref), "hq": hq,
                     "clase": clase})
    tab = pd.DataFrame(rows).sort_values("hq", ascending=False, na_position="last")
    worst = tab.iloc[0] if len(tab) else None
    return {"test": "Cocientes de riesgo HQ (PEC/PNEC)", "tabla": tab,
            "pnec_default": default_pnec_water,
            "interpretacion": (f"Mayor HQ: {worst['polimero']} = {worst['hq']:.2g} ({worst['clase']})"
                               if worst is not None and worst["hq"] == worst["hq"]
                               else "Definir PNECs por polímero para calcular HQ"),
            "nota": "PNEC por defecto referencial de literatura ecotoxicológica acuática; "
                    "ajustar por tamaño de partícula y matriz antes de uso formal"}


def ssd_fit(species_conc, endpoint="NOEC/EC10"):
    """Distribucion de sensibilidad de especies (lognormal) + HC5 + PAF.

    species_conc: serie/dict especie -> concentración critica (misma unidad que el PEC).
    """
    data = pd.Series(species_conc).dropna()
    data = data[data > 0]
    if len(data) < 5:
        return {"test": "SSD", "error": "≥5 especies para ajustar un SSD "
                "(requiere compilar literatura o bioensayos propios)"}
    ln = np.log(data.values)
    mu, sd = ln.mean(), ln.std(ddof=1)
    hc5 = float(np.exp(mu + stats.norm.ppf(0.05) * sd))
    def paf(pec):
        return float(stats.norm.cdf((np.log(pec) - mu) / sd) * 100)
    curve_x = np.exp(np.linspace(ln.min() - 1.5, ln.max() + 1.5, 200))
    curve_y = stats.norm.cdf((np.log(curve_x) - mu) / sd) * 100
    return {"test": f"SSD lognormal ({endpoint}, {len(data)} especies)",
            "hc5": hc5, "mu_log": float(mu), "sd_log": float(sd),
            "paf_func": paf, "curve_x": curve_x, "curve_y": curve_y,
            "especies": data.to_dict(),
            "interpretacion": f"HC5={hc5:.3g} (concentración que protege al 95% de las especies); "
                              "usar PAF(x) para fracción potencialmente afectada en cada PEC"}


def ms_paf(pafs):
    """msPAF por respuesta aditiva (comun en mezclas: suma de PAF de toxico uniforme)."""
    vals = [v for v in pafs if v is not None and v == v]
    if not vals:
        return {"test": "msPAF", "error": "Sin PAFs disponibles"}
    total = float(sum(vals))
    return {"test": "msPAF (suma de PAF, tóxicos uniformes)",
            "mspaf_pct": total,
            "interpretacion": f"msPAF={total:.1f}% de las especies potencialmente afectadas "
                              f"({'riesgo ecológico alto' if total > 5 else 'riesgo bajo-moderado'})"}
