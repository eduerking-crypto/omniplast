"""Biomarcadores de estres oxidativo e IBR integrado (Beliaeff & Burgeot 2002).

IBR (Integrated Biomarker Response): enfoque de desviacion estandar generalizado;
se implementa la variante IBR/2 = sum|SI|/2 usada ampliamente en ecotoxicologia
acuícola (Lam & Gray 2003; Sanchez et al. 2013).
"""
import numpy as np
import pandas as pd


def oxidative_stress_table(bio_df, group_col="zona"):
    """Pivot de biomarcadores por grupo con estadística descriptiva."""
    markers = ["sod_u_mg_prot", "cat_u_mg_prot", "gst_u_mg_prot", "lpo_nmol_tba_mg"]
    sub = bio_df[bio_df["parametro"].isin(markers)]
    if sub.empty:
        return {"test": "Biomarcadores", "error": "Sin datos de biomarcadores (SOD/CAT/GST/LPO)"}
    pivot = sub.pivot_table(index="codigo", columns="parametro", values="valor", aggfunc="mean")
    pivot = pivot.join(sub.drop_duplicates("codigo").set_index("codigo")[group_col])
    summary = pivot.groupby(group_col, observed=True).agg(["count", "mean", "std"]).round(2)
    return {"test": "Panel de estrés oxidativo", "pivot": pivot.reset_index(),
            "resumen": summary, "marcadores": markers,
            "interpretacion": "LPO elevado + antioxidantes alterados indican daño oxidativo "
                              "asociado a contaminación (incl. exposición a MP y aditivos)"}


def ibr_beliaeff(marker_means_treated, marker_means_reference, marker_sds_reference,
                 inhibition_markers=()):
    """IBR/2 a partir de medias del sitio impactado, del control y SD del control.

    inhibition_markers: marcadores cuya DISMINUCION indica dano (ej. antioxidantes deprimidos);
    por defecto se asume que LPO elevado es dano y SOD/CAT/GST elevadas son respuesta compensatoria.
    """
    keys = list(marker_means_reference.keys())
    if not keys:
        return {"test": "IBR", "error": "Sin biomarcadores de referencia"}
    z = []
    for k in keys:
        m_t = marker_means_treated.get(k)
        m_r = marker_means_reference.get(k)
        sd_r = marker_sds_reference.get(k, 1)
        if m_t is None or m_r is None or not sd_r or sd_r <= 0:
            continue
        val = (m_t - m_r) / sd_r
        if k in inhibition_markers:
            val = -val
        z.append((k, val))
    if not z:
        return {"test": "IBR", "error": "Datos incompletos de medias/SD de referencia"}
    ibr = float(sum(abs(v) for _, v in z) / 2)
    quality = ("sin estrés destacable" if ibr < 10 else "estrés leve" if ibr < 30
               else "estrés moderado" if ibr < 60 else "estrés alto / integridad comprometida")
    return {"test": "IBR/2 (Beliaeff & Burgeot 2002)", "ibr": ibr, "perfil": dict(z),
            "interpretacion": f"IBR/2={ibr:.1f} → {quality}",
            "nota": "Clases referenciales; calibrar con controles locales de la especie objetivo"}


def condition_indices(bio_df, group_col="zona"):
    sub = bio_df[bio_df["parametro"].isin(["factor_condicion", "hsi_pct", "gsi_pct"])]
    if sub.empty:
        return {"test": "Índices de condición", "error": "Sin CF/HSI/GSI (requiere datos en peces)"}
    summary = sub.pivot_table(index=group_col, columns="parametro", values="valor",
                              aggfunc=["mean", "std"]).round(3)
    return {"test": "Índices de condición (CF, HSI, GSI)", "resumen": summary,
            "interpretacion": "HSI elevado + CF reducido sugieren代价 energético por exposición "
                              "a contaminantes (revisar con estadística inferencial por sitio)"}


def histopathology_index(lesions_df):
    """lesions_df: columnas [lesion, prevalencia_pct, severidad_0_5, factor_importancia 1-4].

    Indice de alteracion organo-especifico (estilo Bernet et al. 1999):
    IAL = sum(prev_i x sever_i x w_i) / (10 x sum(w_i)) x 10 → 0-10.
    """
    df = lesions_df.copy()
    required = {"lesion", "prevalencia_pct", "severidad"}
    if not required.issubset(df.columns):
        return {"test": "Histopatología", "error": f"Se requieren columnas {required}"}
    if "factor_importancia" not in df.columns:
        df["factor_importancia"] = 2
    df["severidad"] = pd.to_numeric(df["severidad"], errors="coerce").clip(0, 5)
    df["prevalencia_pct"] = pd.to_numeric(df["prevalencia_pct"], errors="coerce").clip(0, 100)
    df["factor_importancia"] = pd.to_numeric(df["factor_importancia"], errors="coerce").fillna(2)
    total_w = df["factor_importancia"].sum() or 1
    score = float((df["prevalencia_pct"] * df["severidad"] * df["factor_importancia"]).sum()
                  / (100 * total_w) * 2)
    score = float(np.clip(score, 0, 10))
    clase = ("sin alteraciones relevantes" if score < 2 else
             "alteraciones leves" if score < 4 else
             "alteraciones moderadas" if score < 6 else
             "alteraciones graves" if score < 8 else "alteraciones muy graves")
    return {"test": "Índice histopatológico (Bernet et al. 1999, adaptado)",
            "indice": score, "tabla": df, "clase": clase,
            "interpretacion": f"Índice={score:.1f}/10 → {clase}"}


def ingestion_metrics(mp_biota_df, group_col="sitio", per_individual_col="abundancia"):
    """Frecuencia de ingestión y carga por individuo a partir de registros MP en biota."""
    if mp_biota_df.empty:
        return {"test": "Ingestión en biota", "error": "Sin registros MP en biota"}
    per_sample = mp_biota_df.groupby(["sample_id", group_col], observed=True)[
        per_individual_col].sum().reset_index()
    per_sample.columns = ["sample_id", group_col, "carga_total"]
    per_sample["ingestion_frecuente"] = per_sample["carga_total"] > 0
    freq = per_sample["ingestion_frecuente"].mean() * 100
    resumen = per_sample.groupby(group_col, observed=True)["carga_total"].agg(
        ["count", "mean", "std", "median"]).round(2)
    return {"test": "Ingestión de MP en biota", "por_muestra": per_sample,
            "frecuencia_ingestion_pct": float(freq), "resumen": resumen,
            "interpretacion": f"{freq:.0f}% de las muestras con MP detectado; carga media="
                              f"{per_sample['carga_total'].mean():.1f} partículas/individuo"}
