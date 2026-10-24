"""Analisis qPCR conforme a MIQE (Bustin et al. 2009): dCt, ddCt, eficiencia, fusión.

Convenciones:
- dCt = Ct_objetivo - Ct_referencia (16S rRNA recomendado como normalizador).
- ddCt = dCt_muestra - dCt_control.
- Cambio relativo de expresion = 2^(-ddCt) (livak) o con eficiencia: (1+E)^(-ddCt).
- Ct > 38 se considera no deteccion practica (Bustin 2009); 35 < Ct <= 38 zona gris.
"""
import numpy as np
import pandas as pd
from scipy import stats


def efficiency_from_slope(slope):
    """E = 10^(-1/pendiente) - 1; ideal 90% a 110% (pendiente -3.58 a -3.1)."""
    if not slope or slope == 0:
        return np.nan
    e = 10 ** (-1 / slope) - 1
    return {"eficiencia_pct": float(e * 100), "valida_miqe": bool(90 <= e * 100 <= 110),
            "pendiente": float(slope),
            "nota": "MIQE: pendiente entre -3.58 y -3.10 ≈ eficiencia 90-110%"}


def standard_curve(dilutions, ct_values):
    """Regresion Ct ~ log10(copia relativa): pendiente, intercepto, R2, eficiencia."""
    data = pd.DataFrame({"log_conc": np.log10(pd.to_numeric(dilutions, errors="coerce")),
                         "ct": pd.to_numeric(ct_values, errors="coerce")}).dropna()
    if len(data) < 3:
        return {"test": "Curva estándar", "error": "≥3 puntos de calibración"}
    res = stats.linregress(data["log_conc"], data["ct"])
    eff = efficiency_from_slope(res.slope)
    r2 = float(res.rvalue ** 2)
    return {"test": "Curva estándar qPCR", "pendiente": float(res.slope),
            "intercepto": float(res.intercept), "r2": r2,
            "eficiencia": eff, "r2_valido_miqe": bool(r2 >= 0.99),
            "interpretacion": f"E={eff['eficiencia_pct']:.1f}%, R²={r2:.4f}, "
                              f"pendiente={res.slope:.3f} "
                              f"({'conforme MIQE' if eff['valida_miqe'] and r2 >= 0.99 else 'REVISAR: fuera de rango MIQE'})"}


def ddct_table(df, target, ref="16s", sample_col="codigo", group_col="zona",
               control_group=None, eff_target=1.0, eff_ref=1.0):
    """df: columnas [sample_col, group_col, target..., ref...] con Ct en filas largas o anchas."""
    t = df.pivot_table(index=["sample_id", sample_col, group_col] if "sample_id" in df.columns
                       else [sample_col, group_col],
                       columns="target", values="ct", aggfunc="mean")
    if target not in t.columns or ref not in t.columns:
        return {"test": "ΔΔCt", "error": f"Se requieren targets '{target}' y '{ref}' con Ct"}
    t = t.dropna(subset=[ref])
    t["dct"] = t[target] - t[ref]
    if control_group:
        ref_dct = t.loc[t.index.get_level_values(group_col) == control_group, "dct"]
        ctrl = float(np.nanmean(ref_dct)) if len(ref_dct) else float(np.nanmean(t["dct"]))
    else:
        ctrl = float(np.nanmin(t["dct"]))
    t["ddct"] = t["dct"] - ctrl
    factor_t = (1 + eff_target) ** (-t["dct"])
    factor_ctrl = (1 + eff_target) ** (-ctrl)
    t["fc"] = factor_t / factor_ctrl
    t["fc"] = t["fc"].clip(upper=1e12)
    return {"test": f"ΔΔCt ({target} vs {ref})", "tabla": t.reset_index(),
            "control_dct": ctrl,
            "no_detectados": int(t[target].isna().sum()),
            "interpretacion": f"FC mediana={np.nanmedian(t['fc']):.2g}; control dCt={ctrl:.2f}"}


def ddct_group_test(ddct_result, group_col="zona"):
    t = ddct_result["tabla"]
    if group_col not in t.columns or t[group_col].nunique() < 2:
        return {"test": "t sobre ΔCt", "error": "Se requieren ≥2 grupos"}
    groups = [g["dct"].dropna().values for _, g in t.groupby(group_col, observed=True)
              if g["dct"].notna().sum() >= 3]
    groups = [x for x in groups if len(x) >= 3]
    if len(groups) < 2:
        return {"test": "t sobre ΔCt", "error": "≥3 muestras con detección por grupo"}
    h, p = stats.kruskal(*groups)
    medias = t.groupby(group_col, observed=True)["dct"].agg(["count", "median", "mean"])
    return {"test": "Kruskal-Wallis sobre ΔCt (método recomendado Yuan et al. 2006)",
            "estadistico": float(h), "p": float(p),
            "resumen_grupos": medias,
            "interpretacion": f"H={h:.3f}, p={p:.4g} → {'diferencias de expresión entre grupos' if p < .05 else 'sin diferencias detectables'}"}


def melt_curve_peaks(temps, fluorescence, height=0.1, prominence=0.05):
    """Detección de Tm: picos en -dF/dT."""
    temps = np.asarray(temps, float)
    f = np.asarray(fluorescence, float)
    if len(temps) < 10:
        return {"test": "Curva de fusión", "error": "Serie insuficiente"}
    neg_deriv = -np.gradient(f, temps)
    peaks, props = scipy_find_peaks(neg_deriv, height=height, prominence=prominence * neg_deriv.max())
    tms = temps[peaks]
    return {"test": "Curva de fusión (Tm calling)", "tms": [float(x) for x in tms],
            "n_picos": int(len(tms)),
            "interpretacion": (f"1 pico (Tm≈{tms[0]:.1f}°C): producto específico"
                                if len(tms) == 1 else
                                f"{len(tms)} picos {np.round(tms, 1).tolist()}: revisar especificidad (dímeros/primer no específico)"
                                if len(tms) > 1 else "Sin picos: falla de amplificación")}


def scipy_find_peaks(x, height=0.1, prominence=None):
    from scipy.signal import find_peaks
    return find_peaks(x, height=height, prominence=prominence)


def miqe_audit(qpcr_df, ntc_ct=None):
    """Auditoría MIQE básica sobre la tabla qPCR."""
    issues = []
    df = qpcr_df
    if "ct" in df.columns:
        ct = pd.to_numeric(df["ct"], errors="coerce")
        n_gris = int(((ct > 35) & (ct <= 38)).sum())
        n_alto = int((ct > 38).sum())
        if n_gris:
            issues.append(f"{n_gris} Ct en zona gris (35-38): reportar como detectables con cautela")
        if n_alto:
            issues.append(f"{n_alto} Ct > 38: se recomienda tratar como no detectados")
        rep_sd = df.groupby(["sample_id", "target"])["ct"].std()
        if (rep_sd > 0.5).any():
            issues.append(f"{int((rep_sd > 0.5).sum())} pocillos con SD de réplicas > 0.5 Ct: repetir")
    if ntc_ct is not None and ntc_ct is not np.nan and ntc_ct < 35:
        issues.append(f"NTC con Ct={ntc_ct}: contaminación de reactivos, invalidar lote")
    checks = [
        ("Curva estándar con ≥5 puntos y R²≥0.99", "Pendiente/eficiencia registradas en datos del run"),
        ("Réplicas técnicas ≥2", "Conteo replicado"),
        ("Controles NTC y no-RT incluidos", ntc_ct is not None),
        ("Eficiencia objetivo y referencia comparables", "Registrar E de ambos"),
    ]
    return {"test": "Auditoría MIQE (Bustin et al. 2009)", "alertas": issues,
            "checklist": [{"criterio": c, "estado": "revisar manualmente" if i < 3 else
                           ("cumplido" if ok else "revisar")} for i, (c, ok) in enumerate(checks)],
            "interpretacion": "; ".join(issues) if issues else "Sin alertas automáticas"}


def gene_environment_association(qpcr_wide, env_wide, gene, env_vars, method="spearman"):
    """Asociación copies/g (log10) de un gen con variables ambientales y carga de MP."""
    from ..statistics import inferential
    merged = qpcr_wide.merge(env_wide, on="sample_id", how="inner", suffixes=("", "_env"))
    merged["log_copies"] = np.log10(merged[gene].clip(lower=1e-3))
    rows = []
    for var in env_vars:
        if var not in merged.columns:
            continue
        res = inferential.correlation_test(merged["log_copies"], merged[var], method=method)
        rows.append({"variable": var, "r": res.get("r"), "p": res.get("p"),
                     "interpretacion": res.get("interpretacion", "")})
    tab = pd.DataFrame(rows).sort_values("p")
    if not tab.empty and tab.iloc[0]["p"] < 0.05:
        top = tab.iloc[0]
        concl = f"Asociación principal: {top['variable']} (r={top['r']:.2f}, p={top['p']:.3g})"
    else:
        concl = "Sin asociaciones significativas tras corrección implícita (revisar Holm)"
    return {"test": f"Asociación gen-ambiente ({gene})", "tabla": tab, "interpretacion": concl}
