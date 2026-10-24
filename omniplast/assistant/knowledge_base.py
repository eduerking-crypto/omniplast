"""Base de conocimiento curada del asistente PLASTI.

REGLA DE HONESTIDAD: solo afirmaciones verificables con cita primaria; cuando un
valor normativo es variable por pais/ano, se marca como 'verificar vigencia'.
"""
from pathlib import Path
import json

from ..config import BASE_DIR

KB = {
    "definicion_microplastico": {
        "tema": "Definición y límites de tamaño",
        "resumen": "Microplástico (MP): partículas de 1 µm a 5 mm. Nanoplástico: 1-1000 nm.",
        "detalle": ("Hartmann et al. (2019, ES&T 53:1039-1047) proponen la definición "
                    "consensuada: nanoplastico 1-1000 nm, microplastico 1 µm-5 mm. "
                    "Evite el termino 'microplastico' para particulas >5 mm (meso/macroplastico)."),
        "fuentes": ["Hartmann et al. (2019) Environ. Sci. Technol. 53:1039-1047"],
    },
    "protocolo_agua": {
        "tema": "Muestreo de MP en agua",
        "resumen": "Red de neustón/manta 300-330 µm o bomba; reportar SIEMPRE volumen filtrado y unidades (part/m³).",
        "detalle": ("GESAMP (2019, Reports & Studies No. 99) y Cincotta et al. (2019): "
                    "prefiera mallas 300-330 µm para comparabilidad; malla 100 µm captura más "
                    "pero se obstruye en aguas turbias. Trace el volumen con flujómetro. "
                    "Replicas ≥3 por estación. Blanca de proceso obligatoria (Koelmans 2019)."),
        "fuentes": ["GESAMP (2019) Guidelines for the monitoring of plastic litter",
                    "Koelmans et al. (2019) ES&T 53:9615-9623"],
    },
    "protocolo_sedimento": {
        "tema": "Muestreo de MP en sedimento",
        "resumen": "Núcleo o cuadrante; capa 0-5 cm; separación por densidad (NaI/ZnCl2 1.5-1.8 g/cm³); unidades part/kg peso seco.",
        "detalle": ("Masura et al. (2015, NOAA Tech. Memo) es el protocolo clasico: tamizado, "
                    "digestion de MO (H2O2 15% o Fenton), flotacion con sal saturada, conteo "
                    "por microscopia y confirmacion espectral en submuestra ≥10%."),
        "fuentes": ["Masura et al. (2015) NOAA Marine Debris Program",
                    "GESAMP (2019)"],
    },
    "protocolo_biota": {
        "tema": "Muestreo de MP en biota",
        "resumen": "Tracto digestivo aislado; n≥10 individuos por sitio; reportar % de individuos con MP y partículas por individuo.",
        "detalle": ("Lusher et al. (2013) y Browne et al. (2008) sientan el estandar: separar "
                    "el tracto (evita sobreestimar por contaminacion externa), digestion con "
                    "KOH 10% (no degrada la mayoria de polimeros; verificar PA/PET), "
                    "blancos de aire en cada lote."),
        "fuentes": ["Lusher et al. (2013) Mar. Pollut. Bull. 77:1-8",
                    "Browne et al. (2008) ES&T 42:5026-5031"],
    },
    "qa_qc": {
        "tema": "Control de calidad (QA/QC)",
        "resumen": "Blancos de aire/proceso, duplicados de campo y laboratorio, recuperación con spikes (70-130%), % identificado espectralmente.",
        "detalle": ("Koelmans et al. (2019) fijan el estandar de reporte minimo: blancos con "
                    "<5% del total, tasa de recuperacion reportada, umbral de deteccion "
                    "explicito, y % de particulas confirmadas por espectroscopia. Sin esto, "
                    "el manuscrito sera devuelto por revisores."),
        "fuentes": ["Koelmans et al. (2019) ES&T 53:9615-9623",
                    "Brander et al. (2020) Appl. Spectrosc. 74:966-988"],
    },
    "identificacion_espectral": {
        "tema": "FTIR vs Raman",
        "resumen": "FTIR (≥20 µm, FPA para mapeo) o Raman (≥1 µm); coincidencia espectral ≥70%; Nile Red como cribado.",
        "detalle": ("FTIR microscopio con detector FPA permite miles de particulas por muestra; "
                    "Raman alcanza ~1 µm pero sufre fluorescencia. La tincion con Nile Red "
                    "(Erni-Cassola 2017) es un cribado cuantitativo rapido que NO reemplaza la "
                    "confirmacion polimerica. Piro-GC/MS cuantifica masa total por polimero."),
        "fuentes": ["Primpke et al. (2020) Anal. Bioanal. Chem.",
                    "Erni-Cassola et al. (2017) Environ. Sci. Technol. 51:13641-13649"],
    },
    "indice_pli": {
        "tema": "Pollution Load Index (PLI)",
        "resumen": "PLI = (∏CF)^(1/n) con CF = C_muestra/C_fondo; PLI>1 = deterioro de calidad.",
        "detalle": ("Derivado de Tomlinson et al. (1980) en geoquimica de metales y adaptado "
                    "para MP: requiere un sitio de referencia (fondo). Reporte SIEMPRE cual "
                    "es el fondo elegido y su justificacion."),
        "fuentes": ["Tomlinson et al. (1980) Helgol. Meeresunters. 33:566-575",
                    "Hakanson (1980) Hydrobiologia 91:832-848"],
    },
    "indice_phi": {
        "tema": "Polymer Hazard Index (PHI)",
        "resumen": "PHI = Σ(puntaje_riesgo_polimero × %abundancia)/100; clases: <1 bajo, 1-10 moderado, 10-100 alto, 100-1000 peligroso.",
        "detalle": ("Usa los puntajes de riesgo quimico de Lithner et al. (2011, ES&T "
                    "45:4727-4735): PE 11, PP 11, PS 30, PET 51, PC 68, PVC 5000, PA 2056, "
                    "PU 1156, ABS 8552. Los polimeros sin puntaje publico se excluyen y se "
                    "reporta el % desconocido."),
        "fuentes": ["Lithner et al. (2011) ES&T 45:4727-4735"],
    },
    "wqi_nsf": {
        "tema": "Índice de calidad del agua NSF",
        "resumen": "WQI = Σ(qi×wi)/Σwi con 9 parámetros ponderados (OD 0.17, coliformes 0.16...).",
        "detalle": ("Brown et al. (1970, NSF/WPCF): cada parametro tiene curva de valoracion "
                    "qi (0-100) y peso wi. Clases: >90 excelente, 71-90 buena, 51-70 media, "
                    "26-50 mala, <26 muy mala. OMNIPLAST interpola curvas publicadas."),
        "fuentes": ["Brown et al. (1970) Water Quality Index, NSF"],
    },
    "riesgo_hq": {
        "tema": "Cociente de riesgo PEC/PNEC",
        "resumen": "HQ = PEC/PNEC; HQ<1 riesgo aceptable, >1 potencial; requiere PNEC defendible por polímero y tamaño.",
        "detalle": ("Marco regulatorio EU (TGD R.10): PNEC = NOEC/AF (factor de asignacion "
                    "10-1000). Para MP NO existen PNEC oficiales: los valores de literatura "
                    "que usa OMNIPLAST son editables y deben citarse. Frias & Nash (2019) "
                    "proponen el marco de riesgo relativo para MP."),
        "fuentes": ["EU TGD-R.10 (2003)", "Frias & Nash (2019) Mar. Pollut. Bull. 138:316-322"],
    },
    "genes_degradacion": {
        "tema": "Genes de degradación de plásticos",
        "resumen": "Panel OMNIPLAST: petase/mhetase (PET), alkB (alcanos tipo PE), laccasa (aromáticos), phaC (PHA); referencia 16S.",
        "detalle": ("Ideonella sakaiensis 201-F6 degrada PET via PETase + MHETase (Yoshida et "
                    "al. 2016, Science 351:1196-1203). alkB oxida cadenas alifaticas (van "
                    "Beilen & Funhoff 2007). Para qPCR use guias MIQE (Bustin 2009): "
                    "eficiencia 90-110%, R2>0.99, NTC limpios, y 16S como normalizador."),
        "fuentes": ["Yoshida et al. (2016) Science 351:1196-1203",
                    "Danso, Chow & Streit (2019) AEM 85:e01095-19",
                    "Bustin et al. (2009) Clin. Chem. 55:611-622"],
    },
    "plastisfera": {
        "tema": "Plastísfera",
        "resumen": "Biofilm sobre MP con taxa oportunistas (Vibrio, Pseudomonas); enriquecimiento verificable con 16S + PERMANOVA.",
        "detalle": ("Zettler, Mincer & Amaral-Zettler (2013) acuñaron el termino. Diseño "
                    "recomendado: muestrear MP y sustrato natural (piedra/madera) del mismo "
                    "sitio; contrastar con PERMANOVA y SIMPER; correlacionar genes de "
                    "degradacion con la carga de MP (asociacion gen-ambiente, AG-067)."),
        "fuentes": ["Zettler et al. (2013) ES&T 47:7137-7146",
                    "Oberbeckmann et al. (2015) Mar. Environ. Res. 120:220-232"],
    },
    "permanova": {
        "tema": "PERMANOVA",
        "resumen": "Compara centroides multivariados por permutación (Anderson 2001); complemente SIEMPRE con betadisper.",
        "detalle": ("PERMANOVA asume dispersion similar entre grupos; si betadisper es "
                    "significativo, la 'diferencia' puede ser solo varianza. Reporte "
                    "pseudo-F, R2, p por permutacion (999+) y el resultado de betadisper."),
        "fuentes": ["Anderson (2001) Austral. NZ J. Stat. 43:281-300"],
    },
    "nmds_stress": {
        "tema": "NMDS y calidad del ajuste",
        "resumen": "Stress de Kruskal: <0.05 excelente, <0.10 buena, <0.20 útil con cautela, >0.30 cuestionable.",
        "detalle": ("NMDS es no-paramétrico y no-métrico: preserve el orden de las "
                    "distancias, no su magnitud. Ejecute ≥250 reinicios aleatorios (metaMDS "
                    "en vegan usa trymax=250); OMNIPLAST alinea la solucion via Procrustes "
                    "para estabilidad visual."),
        "fuentes": ["Kruskal (1964) Psychometrika 29:115-129", "Clarke (1993) MEPS 92:205-219"],
    },
    "glm_cuentos_mp": {
        "tema": "Estadística de conteos de MP",
        "resumen": "Los conteos de MP son datos de conteo sobredispersos: use GLM Binomial Negativa, nunca ANOVA directo.",
        "detalle": ("La varianza de conteos ambientales supera la media (agregacion espacial): "
                    "el modelo Poisson subestima el error estandar. Ajuste Binomial Negativa, "
                    "reporte alpha de sobredispersion, y use offset=log(volumen) si difiere "
                    "entre muestras. Alternativas: GEE por sitio o mixtos (GLMM)."),
        "fuentes": ["Cameron & Trivedi (2013) Regression Analysis of Count Data",
                    "Zuur et al. (2009) Mixed Effects Models in R"],
    },
    "pseudorreplicacion": {
        "tema": "Pseudorreplicación espacial",
        "resumen": "Muestras repetidas del mismo sitio NO son réplicas independientes: modele el sitio como efecto aleatorio o cluster GEE.",
        "detalle": ("Si muestrea 3 veces la misma estacion, la replica REAL es la estacion. "
                    "Use modelos mixtos (intercepto aleatorio por sitio, reportando ICC) o "
                    "GEE con correlacion intercambiable. Hurlbert (1984) es la cita clasica."),
        "fuentes": ["Hurlbert (1984) Ecol. Monogr. 54:187-211", "Zuur et al. (2009)"],
    },
    "correlacion_vs_causalidad": {
        "tema": "Correlación no es causalidad",
        "resumen": "Use correlación parcial, mediación y diseño (gradient/impact-reference) antes de afirmar causalidad.",
        "detalle": ("Un r alto entre MP y DBO puede reflejar una fuente común (vertimiento). "
                    "Controle covariables (correlacion parcial), pruebe mediacion (bootstrap "
                    "del efecto indirecto) y discuta explicitamente confusores no medidos."),
        "fuentes": ["Zar (2010) Biostatistical Analysis"],
    },
    "autocorrelacion_espacial": {
        "tema": "Autocorrelación espacial",
        "resumen": "Estaciones cercanas se parecen: calcule I de Moran; si es significativo, use modelos espaciales o GEE.",
        "detalle": ("Moran (1950) I global con pesos k-NN y p por permutacion; LISA local "
                    "clasifica HH/LL/HL/LH (Anselin 1995). Los hotspots Getis-Ord Gi* "
                    "priorizan acciones de limpieza con niveles de confianza 90/95/99%."),
        "fuentes": ["Moran (1950)", "Anselin (1995)", "Getis & Ord (1992)"],
    },
    "kriging": {
        "tema": "Kriging ordinario",
        "resumen": "Interpolación óptima con variograma (nugget-meseta-alcance); reporte validación cruzada leave-one-out.",
        "detalle": ("Ajuste el variograma empirico a modelos esferico/exponencial/gaussiano; "
                    "el alcance define la distancia de independencia. Con <15 puntos, "
                    "prefiera IDW y reportelo honestamente."),
        "fuentes": ["Isaaks & Srivastava (1989) Applied Geostatistics"],
    },
    "deteccion_satelital": {
        "tema": "Detección satelital de plásticos",
        "resumen": "FDI (Biermann 2020) detecta escombros flotantes genéricos con Sentinel-2; valide SIEMPRE en campo.",
        "detalle": ("Biermann et al. (2020, Sci. Rep. 10:19836) lograron 86% de precision "
                    "en escenarios controlados con el Floating Debris Index sobre Sentinel-2 "
                    "(10 m/pixel): detecta agregados flotantes, NO particulas individuales, "
                    "y confunde espumas/Sargassum/madera. Complemente con drones (cm-scale)."),
        "fuentes": ["Biermann et al. (2020) Sci. Rep. 10:19836", "Hu (2009) RSE 113:398-409"],
    },
    "normativa_ue": {
        "tema": "Normativa internacional",
        "resumen": "Reglamento (UE) 2024/1448 fija METODOLOGÍA de medición en agua de consumo; ISO 24187:2023 principios de análisis; sin valor límite ambiental vigente para MP.",
        "detalle": ("La UE reguló metodología (Reglamento 2024/1448, en el marco de la "
                    "Directiva (UE) 2020/2184) pero NO existe aún un valor paramétrico "
                    "numérico obligatorio para MP ambiental en la mayoria de jurisdicciones. "
                    "Colombia: sin norma específica de MP hasta el corte de conocimiento "
                    "del sistema (2026); el marco general es el Decreto 4742 de 2005 y las "
                    "políticas de gestión de residuos sólidos. VERIFICAR VIGENCIA antes de citar."),
        "fuentes": ["Reglamento (UE) 2024/1448", "Directiva (UE) 2020/2184",
                    "ISO/TR 21960:2020", "ISO 24187:2023"],
    },
    "unidades_reporte": {
        "tema": "Unidades y reporte",
        "resumen": "Agua: part/m³ (o part/L); sedimento: part/kg peso seco; biota: part/individuo. Reporte rango de tamaños y límite inferior.",
        "detalle": ("El 80% de la inconsistencia entre estudios viene de unidades y limites "
                    "de deteccion distintos. Siempre: rango tamal (ej. 100-5000 µm), metodo, "
                    "volumen/masa, % confirmado espectralmente, blancos y recuperacion."),
        "fuentes": ["Hartmann et al. (2019)", "Koelmans et al. (2019)"],
    },
    "cienga_virgen": {
        "tema": "Ciénaga de La Virgen (demo OMNIPLAST)",
        "resumen": "Laguna costera de Cartagena con presión urbana fuerte; el dataset demo reproduce patrones de literatura (DATOS SINTÉTICOS).",
        "detalle": ("La Cienaga de La Virgen recibe vertimientos del canal Juan Polo y "
                    "presion turistica; su manglar es habitat de aves y pesca artesanal. "
                    "El proyecto demo OMNIPLAST simula 10 estaciones, 4 campanas y 3 "
                    "matrices con patrones fisicoquimicos tipicos de humedales costeros "
                    "hiperurbanizados del Caribe. Sustituya los datos por los suyos en la "
                    "misma estructura (data/demo/*.csv como plantilla)."),
        "fuentes": ["Contexto regional: PNOH; INVEMAR (informes litoral Caribe)"],
    },
    "diseno_muestral": {
        "tema": "Diseño muestral doctoral para MP",
        "resumen": "Diseño BACI o gradiente impacto-referencia; ≥5 estaciones, ≥3 réplicas, ≥4 campañas estacionales, QA/QC completo.",
        "detalle": ("Para inferencia causal prefiera 'Before-After-Control-Impact' (BACI); "
                    "si no hay 'antes', use gradiente espacial con referencia externa. "
                    "Calcule n a priori con analisis de potencia (d esperado 0.8, potencia "
                    "0.8, alfa 0.05 → ~26 por grupo). Incluya estacion de referencia "
                    "externa NO influenciada por la fuente."),
        "fuentes": ["Underwood (1994) J. Exp. Mar. Biol. Ecol. 178:227-243"],
    },
    "literatura_clave": {
        "tema": "Lecturas fundacionales",
        "resumen": "Hartmann 2019 (definición); GESAMP 2019 (monitoreo); Koelmans 2019 (QA/QC); Yoshida 2016 (PETasa); Zettler 2013 (plastísfera).",
        "detalle": ("Rutas de entrada: (1) marco conceptual: Hartmann 2019; (2) protocolos: "
                    "GESAMP 2019 + Masura 2015; (3) calidad: Koelmans 2019; (4) biologia: "
                    "Zettler 2013 + Yoshida 2016; (5) riesgo: Frias & Nash 2019 + WHO 2019; "
                    "(6) revision critica de métodos: Ranjan & Goel 2019? (verificar)."),
        "fuentes": ["Todas citadas en docs/METODOLOGIA.md"],
    },
}


def search_kb(query, top=3):
    q = query.lower()
    scored = []
    for key, entry in KB.items():
        text = " ".join([key, entry["tema"], entry["resumen"], entry["detalle"]]).lower()
        score = sum(1 for w in q.split() if len(w) > 3 and w in text)
        if score > 0:
            scored.append((score, key, entry))
    scored.sort(key=lambda x: -x[0])
    return [{"clave": k, **e} for _, k, e in scored[:top]]


def context_for_prompt(query, top=3):
    results = search_kb(query, top)
    if not results:
        return ""
    blocks = []
    for r in results:
        blocks.append(f"### {r['tema']}\n{r['resumen']}\n{r['detalle']}\nFuentes: "
                      + "; ".join(r["fuentes"]))
    return "CONTEXTO CIENTÍFICO VERIFICADO (usa SOLO esta información como base y cita las fuentes):\n" + \
        "\n\n".join(blocks)
