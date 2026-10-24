"""Asesor metodologico: arbol de decision estadistico y diseno experimental."""
from ..core.statistics.bayes import power_analysis
from ..core.statistics.inferential import cohen_d


def recommend_test(objetivo="comparar", tipo="continua", n_grupos=2, pareado=False,
                   normalidad=None, heterocedasticidad=None, conteos=False,
                   multivariado=False, repetido=False):
    """Arbol de decision estadistico con justificacion y alternativa robusta."""
    if conteos:
        return {"prueba": "GLM Binomial Negativa (o Poisson si α≈0)",
                "justificacion": "Los conteos de MP son sobredispersos y no normales por "
                                 "naturaleza; el ANOVA/t sobre conteos violan supuestos.",
                "alternativa": "GEE por sitio o GLMM si hay agrupación espacial/temporal",
                "supuestos": "independencia de residuos; verificar α de sobredispersión",
                "ref": "Cameron & Trivedi (2013); Zuur et al. (2009)"}
    if multivariado:
        if normalidad:
            return {"prueba": "MANOVA / RDA (si hay gradiente ambiental)",
                    "justificacion": "Respuestas multivariadas normales con homocedasticidad "
                                     "multivariada (Box M) y linealidad.",
                    "alternativa": "PERMANOVA (más robusto, sin supuestos distribucionales)",
                    "supuestos": "normalidad multivariada, Box M (sensible)",
                    "ref": "Anderson (2001); Legendre & Legendre (2012)"}
        return {"prueba": "PERMANOVA (Bray-Curtis) + betadisper",
                "justificacion": "Libre de distribución para datos de comunidad; betadisper "
                                  "verifica que la diferencia no sea solo dispersión.",
                "alternativa": "ANOSIM como complemento descriptivo de separación",
                "supuestos": "dispersión similar entre grupos (verificar con betadisper)",
                "ref": "Anderson (2001)"}
    if n_grupos == 1:
        return {"prueba": "Prueba t de una muestra o IC de referencia",
                "justificacion": "Comparar contra un valor de referencia publicado (p. ej. fondo).",
                "alternativa": "Bootstrap del IC si n<10 o no normalidad",
                "supuestos": "normalidad (Shapiro) o n grande",
                "ref": "Zar (2010)"}
    if n_grupos == 2:
        if pareado:
            return {"prueba": "t pareada (normal) o Wilcoxon (no normal)",
                    "justificacion": "Las mediciones pareadas por sitio/remuestreo controlan "
                                     "la variabilidad entre unidades.",
                    "alternativa": "Friedman si >2 condiciones pareadas",
                    "supuestos": "normalidad de las diferencias (para t)",
                    "ref": "Zar (2010)"}
        if normalidad is None:
            return {"prueba": "Primero diagnostique: Shapiro-Wilk + Levene (use AG-033)",
                    "justificacion": "La elección t vs Mann-Whitney depende del diagnóstico.",
                    "alternativa": "Con dudosa normalidad y n<30, Mann-Whitney + Cliff's δ",
                    "supuestos": "—",
                    "ref": "Zar (2010)"}
        if normalidad:
            if heterocedasticidad:
                return {"prueba": "t de Welch (no asume varianzas iguales)",
                        "justificacion": "Robusta ante desigualdad de varianzas; siempre "
                                         "preferible a la t clásica.",
                        "alternativa": "Bootstrap de la diferencia de medias",
                        "supuestos": "normalidad por grupo; independencia",
                        "ref": "Welch (1947)"}
            return {"prueba": "t de Student + d de Cohen",
                    "justificacion": "Normalidad y homocedasticidad verificadas.",
                    "alternativa": "t de Welch (más conservadora)",
                    "supuestos": "normalidad, homocedasticidad",
                    "ref": "Zar (2010)"}
        return {"prueba": "Mann-Whitney U + δ de Cliff",
                "justificacion": "Sin normalidad, compara medianas con potencia decente "
                                 "y reporta tamaño de efecto no paramétrico.",
                "alternativa": "t de Welch con bootstrap si n≥30 (CLT)",
                "supuestos": "independencia; formas similares de distribución",
                "ref": "Zar (2010)"}
    if n_grupos > 2:
        if repetido:
            return {"prueba": "ANOVA de medidas repetidas / Friedman",
                    "justificacion": "Múltiples mediciones de la misma unidad (sitio×campaña).",
                    "alternativa": "Modelo mixto (intercepto+pendiente aleatoria) — superior",
                    "supuestos": "esfericidad (Mauchly) para RM-ANOVA",
                    "ref": "Zuur et al. (2009)"}
        if normalidad:
            if heterocedasticidad:
                return {"prueba": "ANOVA-Welch + Games-Howell post-hoc",
                        "justificacion": "Post-hoc válido con varianzas desiguales y n "
                                         "desiguales entre zonas.",
                        "alternativa": "Kruskal-Wallis + Dunn",
                        "supuestos": "normalidad por grupo",
                        "ref": "Games & Howell (1976)"}
            return {"prueba": "ANOVA + Tukey HSD + η²/ω²",
                    "justificacion": "Comparación omnibus clásica con post-hoc equilibrado.",
                    "alternativa": "Kruskal-Wallis + Dunn (Holm) como robusta",
                    "supuestos": "normalidad, homocedasticidad",
                    "ref": "Fisher (1925); Tukey (1949)"}
        return {"prueba": "Kruskal-Wallis + Dunn (corrección Holm)",
                "justificacion": "Omnibus no paramétrico con post-hoc corregido por "
                                 "comparaciones múltiples.",
                "alternativa": "PERMUTA: ANOVA por permutaciones si hay vinculación con covariables",
                "supuestos": "independencia",
                "ref": "Kruskal & Wallis (1952); Dunn (1964)"}
    return {"prueba": "Describa el objetivo con más precisión en el chat.",
            "justificacion": "", "alternativa": "", "supuestos": "", "ref": ""}


def sampling_design(d_esperado=0.8, potencia=0.8, estaciones=6, replicas_min=3):
    n = power_analysis(effect_d=d_esperado, power=potencia)
    n_grupo = n.get("n_por_grupo", 26)
    return {"n_por_grupo": n_grupo,
            "campanas_recomendadas": 4,
            "estaciones_recomendadas": estaciones,
            "replicas_por_estacion": max(replicas_min, 3),
            "total_muestras_aprox": n_grupo * 2 + estaciones * replicas_min,
            "estructura": "Gradiente impacto-referencia (o BACI si existe línea base 'antes')",
            "checklist": ["Estación de referencia externa", "Blancos de proceso por campaña",
                          "Duplicados de campo ≥10%", "Registro de volumen filtrado",
                          "Confirmación espectral ≥10% de partículas"],
            "justificacion": f"Con d={d_esperado} y potencia {potencia:.0%} se requieren "
                             f"≈{n_grupo} observaciones por grupo (aprox. normal, α=0.05)."}


def effect_size_interpretation(d):
    d = abs(d)
    if d < 0.2:
        return "trivial (<0.2)"
    if d < 0.5:
        return "pequeño (0.2-0.5)"
    if d < 0.8:
        return "mediano (0.5-0.8)"
    return "grande (≥0.8)"
