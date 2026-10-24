"""Registro de los 100 agentes especialistas de OMNIPLAST.

Cada agente es un analista virtual con dominio propio que ejecuta un flujo de
trabajo completo sobre los datos del proyecto activo (o entrega la guia y plantilla
exacta cuando falta un insumo especifico). Niveles: Basico / Intermedio / Doctoral.
"""

CATEGORIES = {
    "A": ("Físico-química del agua", "#05BFDB"),
    "B": ("Sedimentos", "#7C9EB2"),
    "C": ("Microplásticos — núcleo", "#00A896"),
    "D": ("Estadística inferencial", "#0A2540"),
    "E": "Machine Learning",
    "F": ("Genética & molecular", "#E9B44C"),
    "G": ("Plastísfera & biodegradación", "#3E9B4F"),
    "H": ("Biota & ecotoxicología", "#FF6B6B"),
    "I": ("Teledetección & GIS", "#4CC9F0"),
    "J": ("Riesgo & One Health", "#F07CAD"),
}
CATEGORIES["E"] = ("Machine Learning & IA", "#9B5DE5")

_SPEC = """
AG-001|A|Perfil térmico y estratificación|Intermedio|Temperatura por zona, campaña y profundidad, con comparación inferencial y detección de estratificación.|Hartmann et al. (2019)|
AG-002|A|Déficit de oxígeno disuelto|Intermedio|Saturación de O2 (fórmula de Benson-Krause), hipoxia por sitio y zonas críticas.|US EPA (1986); Benson & Krause (1984)|
AG-003|A|Carga orgánica DBO/DQO|Intermedio|DBO5, DQO y su relación (biodegradabilidad) por zona y temporada.|Sawyer et al. (2003)|
AG-004|A|Estado nutricional y ratio N:P|Intermedio|DIN (NO3+NO2+NH4), DIP y desviación del ratio de Redfield.|Redfield et al. (1963)|
AG-005|A|Salinidad y conectividad estuarina|Intermedio|Gradiente salino, mezcla y tipo de estuario (SIMPER de zonas).|Bowden (1980)|
AG-006|A|Turbidez, SST y atenuación de luz|Intermedio|Relación turbidez-SST y estimación empírica de la zona eufólica (Kd).|Kirk (2011)|
AG-007|A|Índice WQI-NSF por estación|Doctoral|WQI de 9 parámetros con curvas de valoración NSF y clasificación por estación.|Brown et al. (1970)|
AG-008|A|Estado trófico TSI (Carlson)|Doctoral|TSI por clorofila y fósforo total; clasificación trófica del sistema.|Carlson (1977)|
AG-009|A|Calidad microbiológica sanitaria|Intermedio|Coliformes y E. coli vs umbrales de contacto primario; % de excedencias.|US EPA RWQC (2012); OMS (2003)|
AG-010|A|Scorecard fisicoquímico integrado|Doctoral|Tasa de excedencia por parámetro + radar de estado por sitio (scoring multicriterio).|Elaboración OMNIPLAST|
AG-011|B|Granulometría vs microplásticos|Intermedio|Relación arena% y finos con la carga de MP en sedimento.|Besseling et al. (2017)|
AG-012|B|Materia orgánica vs MP|Intermedio|TOC como predictor de MP en sedimento (correlación + regresión).|Besseling et al. (2017)|
AG-013|B|Perfiles verticales en la columna|Intermedio|MP a 0/1.5/3 m; hundimiento diferencial por polimero (densidad).|Kowalski et al. (2016)|
AG-014|B|Flujo de acumulación sedimentaria|Doctoral|Estimación de flujo de MP (part/kg/día) con supuestos explícitos de área y tasa de sedimentación.|Elaboración OMNIPLAST|
AG-015|B|PLI sedimentario por sitio|Doctoral|Factor de contaminación y Pollution Load Index por estación (fondo = referencia).|Hakanson (1980); Tomlinson (1980)|
AG-016|B|Diagnóstico redox-proxy|Basico|MO%, textura y salinidad como proxies de condición redox del sedimento.|Berner (1981)|
AG-017|B|Clasificación textural (Folk)|Basico|Diagrama arena-finos y clase textural por estación.|Folk (1954)|
AG-018|B|Asociación MP-ambiente en sedimento|Intermedio|Matriz de correlaciones MP×variables del sedimento con corrección por empates.|Elaboración OMNIPLAST|
AG-019|C|Abundancia espacial de MP en agua|Doctoral|Boxplot por zona, KW, tamaño de efecto y ranking de sitios.|Hartmann et al. (2019)|
AG-020|C|Distribución de tamaños|Doctoral|Banda de tamaño con modelo log-decay y fracción <300 µm (riesgo de ingestión).|Kooi et al. (2016)|
AG-021|C|Composición de polímeros|Intermedio|Participación por polímero y matriz (torta + ranking), señal de fuentes.|Erni-Cassola et al. (2019)|
AG-022|C|Morfologías y dominancia de fibras|Intermedio|Fibras vs fragmentos por zona; índice de fibra urbana.|Cole et al. (2011)|
AG-023|C|Colores y foto-degradación|Basico|Espectro de colores; proporción de colores cálidos como proxy de exposición solar.|Enders et al. (2015)|
AG-024|C|Perfil vertical de MP|Intermedio|MP vs profundidad con GLM NB y sedimentación por densidad de polímero.|Kowalski et al. (2016)|
AG-025|C|Comparación estacional|Intermedio|Secas vs lluvias por matriz con pruebas pareadas por sitio.|GESAMP (2019)|
AG-026|C|Hotspots espaciales (Gi*)|Doctoral|Estadístico Getis-Ord Gi* por sitio con niveles de confianza.|Getis & Ord (1992)|
AG-027|C|Comparación entre matrices|Intermedio|Co-variación agua-sedimento-biota por sitio (Spearman).|Hartmann et al. (2019)|
AG-028|C|Ingestión en biota|Intermedio|Frecuencia de ingestión y carga por individuo y especie.|Lusher et al. (2013)|
AG-029|C|Stock lagunar total|Doctoral|Estimación de stock (partículas y masa) con incertidumbre Monte Carlo de volumen y concentración.|Elaboración OMNIPLAST|
AG-030|C|Auditoría QA/QC de laboratorio|Doctoral|Blancos, duplicados, tasas de identificación y alertas de contaminación de proceso.|Koelmans et al. (2019)|
AG-031|C|Tamaño efectivo por polímero|Intermedio|Mediana de tamaño por polímero con comparación no paramétrica.|Elaboración OMNIPLAST|
AG-032|C|Ratio fibra/fragmento vs carga urbana|Intermedio|Fibras como trazadora de vertimientos; correlación con DBO.|Browne et al. (2011)|
AG-033|D|Diagnóstico de distribución|Basico|Normalidad (Shapiro/AD), sesgo, curtosis y diagrama de decisión paramétrico/no paramétrico.|Zar (2010)|
AG-034|D|Transformaciones óptimas|Intermedio|Comparación log/sqrt/Box-Cox por p de Shapiro y recomendación aplicada.|Box & Cox (1964)|
AG-035|D|t de Welch con efecto|Basico|Comparación de dos grupos con d de Cohen e IC de la diferencia.|Welch (1947)|
AG-036|D|ANOVA + Tukey HSD|Basico|Diferencias entre zonas con post-hoc y eta².|Fisher (1925); Tukey (1949)|
AG-037|D|Kruskal-Wallis + Dunn|Intermedio|Alternativa robusta con comparaciones post-hoc corregidas por Holm.|Kruskal & Wallis (1952); Dunn (1964)|
AG-038|D|Games-Howell post-hoc|Intermedio|Post-hoc válido bajo heterocedasticidad y n desiguales.|Games & Howell (1976)|
AG-039|D|Comparaciones pareadas|Intermedio|Wilcoxon/Friedman por sitio entre campañas (diseño bloqueado).|Zar (2010)|
AG-040|D|Correlaciones multivariadas|Intermedio|Matriz MP×ambiente con heatmap de r y señalización de p.|Zar (2010)|
AG-041|D|Correlación parcial|Intermedio|Efecto de MP controlando covariables (evitar correlaciones espurias).|Zar (2010)|
AG-042|D|Regresión múltiple OLS|Intermedio|Modelo multivariado con SE robustas HC3, VIF e inferencia completa.|MacKinnon & White (1985)|
AG-043|D|Regresión cuantílica|Doctoral|Efecto heterogéneo del predictor en los extremos de la distribución de MP.|Koenker & Bassett (1978)|
AG-044|D|GLM Poisson para conteos|Intermedio|Modelado de conteos de MP con razones de cambio por unidad de covariable.|McCullagh & Nelder (1989)|
AG-045|D|GLM Binomial Negativa|Doctoral|Concurrencia de sobredispersión (α) y modelado correcto de conteos ambientales.|Cameron & Trivedi (2013)|
AG-046|D|Modelo mixto por sitio|Doctoral|Intercepto aleatorio por estación y ICC (pseudorreplicación evitada).|Zuur et al. (2009)|
AG-047|D|GEE por conglomerados|Doctoral|Inferencia poblacional robusta con estructura de correlación por sitio.|Liang & Zeger (1986)|
AG-048|D|PCA ambiental|Intermedio|Ordenación de sitios por variables fisicoquímicas, biplot y Kaiser.|Pearson (1901)|
AG-049|D|Clustering jerárquico de estaciones|Intermedio|Dendrograma de similitud fisicoquímica con cophenetic.|Sokal & Rohlf (1962)|
AG-050|D|PERMANOVA + ANOSIM|Doctoral|Diferencias multivariadas de microbioma/ambiente entre zonas (999 permutaciones).|Anderson (2001)|
AG-051|D|SIMPER de contribuciones|Doctoral|Descomposición de la disimilitud: taxones/variables que separan los grupos.|Clarke (1993)|
AG-052|D|NMDS con stress de Kruskal|Doctoral|Ordenación no métrica del microbioma con calidad de ajuste.|Kruskal (1964)|
AG-053|E|Random Forest: clasificación de impacto|Doctoral|Clasificación de zonas con CV estratificada, importancia por permutación y matriz de confusión.|Breiman (2001)|
AG-054|E|Random Forest: predecir MP|Doctoral|Modelo predictivo de carga de MP desde el ambiente con R²/RMSE.|Breiman (2001)|
AG-055|E|Importancia de variables|Intermedio|Ranking por permutación del modelo de bosques (interpretabilidad).|Breiman (2001)|
AG-056|E|SVM con kernel RBF|Intermedio|Clasificación de márgenes con estandarización previa.|Cortes & Vapnik (1995)|
AG-057|E|k-means + silueta|Intermedio|Segmentación de estaciones y k óptimo por silueta.|Rousseeuw (1987)|
AG-058|E|DBSCAN de atípicos espaciales|Intermedio|Detección de estaciones anómalas sin número previo de grupos.|Ester et al. (1996)|
AG-059|E|Anomalías por reconstrucción PCA|Intermedio|Error de reconstrucción para detectar muestras multivariadas atípicas.|Chandola et al. (2009)|
AG-060|E|Gradient Boosting de riesgo|Doctoral|Modelo aditivo para clasificar riesgo alto/bajo con importancia de variables.|Friedman (2001)|
AG-061|E|Bootstrap de métricas|Intermedio|Intervalos de confianza bootstrap de la carga media por sitio.|Efron (1979)|
AG-062|E|Curvas ROC/PR y AUC|Intermedio|Capacidad discriminativa del modelo con OR por predictor.|Hanley & McNeil (1982)|
AG-063|F|Expresión qPCR (ΔΔCt)|Doctoral|Tabla completa de ΔCt/ΔΔCt/FC por muestra con prueba de grupos (Yuan et al.).|Bustin et al. (2009); Yuan et al. (2006)|
AG-064|F|Curva estándar y eficiencia|Intermedio|Pendiente, R², eficiencia MIQE y veredicto del ensayo.|Bustin et al. (2009)|
AG-065|F|Curvas de fusión (Tm)|Intermedio|Detección de picos de Tm y alerta de inespecificidad (plantilla si faltan datos).|Bustin et al. (2009)|
AG-066|F|Genotipado por agrupamiento|Intermedio|Llamadas alélicas 2D, call rate y verificación HWE.|Livak et al. (1995)|
AG-067|F|Asociación gen-ambiente|Doctoral|Copias log10 de genes de degradación vs variables ambientales (tabla de r).|Elaboración OMNIPLAST|
AG-068|F|Diversidad alfa del microbioma|Intermedio|S, Chao1, Shannon, Simpson y Pielou por muestra/zona con comparación inferencial.|Chao (1984); Hill (1973)|
AG-069|F|Diversidad beta + PCoA|Intermedio|Distancias Bray-Curtis y ordenación con varianza explicada.|Legendre & Legendre (2012)|
AG-070|F|PERMANOVA del microbioma|Doctoral|Efecto de la zona sobre la estructura de comunidad (con betadisper).|Anderson (2001)|
AG-071|F|Abundancia diferencial 16S|Doctoral|DESeq2-lite: size factors + MWU + FDR; enriquecidos por grupo.|Love et al. (2014)|
AG-072|F|Filogenética de aislados|Intermedio|Alineamiento por pares + árbol NJ desde FASTA del usuario (o guía de pipeline).|Saitou & Nei (1987)|
AG-073|F|Redes de co-ocurrencia|Doctoral|Grafo con hubs, modularidad y comunidades (detección de consortios degradadores).|Faust & Raes (2012)|
AG-074|F|Anotación funcional (KO)|Intermedio|Panel de genes y su mapeo conceptual a rutas de degradación (referencias primarias).|Yoshida et al. (2016); Danso et al. (2019)|
AG-075|G|Enriquecimiento de taxa en MP|Doctoral|Taxones enriquecidos en sitios con alta carga (z-score por taxón).|Zettler et al. (2013)|
AG-076|G|Core microbioma de la plastísfera|Intermedio|Prevalencia y abundancia; núcleo del sistema de estudio.|Shade & Handelsman (2012)|
AG-077|G|Curvas de biodegradación (Gompertz)|Doctoral|Ajuste del modelo de Gompertz modificado a series de crecimiento sobre plástico.|Zwietering et al. (1990)|
AG-078|G|Biodegradabilidad por polímero|Intermedio|Ranking de biodegradabilidad relativa según genes detectados y literatura.|Danso et al. (2019)|
AG-079|G|Panel enzimático correlacionado|Intermedio|Correlaciones entre genes del panel (petase-mhetase-laccasa-alkB-phaC).|Yoshida et al. (2016)|
AG-080|G|Biomarcadores LEfSe-lite|Doctoral|Taxones discriminantes por zona con efecto estandarizado y FDR.|Segata et al. (2011)|
AG-081|G|Co-marcador de resistencia (intI1)|Intermedio|Asociación conceptual MP×AMR con guía de medición de intI1 (plantilla).|Wright et al. (2020)|
AG-082|G|Score de formadores de biofilm|Basico|Abundancia ponderada de taxa formadores de biofilm conocidos sobre MP.|Oberbeckmann et al. (2015)|
AG-083|H|Ingestión por especie|Intermedio|Carga de MP en ostras vs peces con comparación no paramétrica.|Lusher et al. (2013)|
AG-084|H|Intensidad de contaminación en tracto|Intermedio|Partículas por individuo, distribución y % con >5 partículas.|Lusher et al. (2013)|
AG-085|H|Fibras vs fragmentos en biota|Intermedio|Prueba χ² de morfología ingerida por especie (trampa selectiva de fibras).|Browne et al. (2008)|
AG-086|H|Biomarcadores de estrés oxidativo|Doctoral|SOD/CAT/GST/LPO por zona con KW y dirección del daño.|Regoli & Giuliani (2014)|
AG-087|H|IBR integrado (Beliaeff-Burgeot)|Doctoral|Índice integrado de respuesta biológica con perfil estandarizado por sitio.|Beliaeff & Burgeot (2002)|
AG-088|H|Índices de condición (CF/HSI/GSI)|Intermedio|Condición energética de la biota por zona de impacto.|Bagenal (1978)|
AG-089|H|Histopatología escalonada|Intermedio|Índice de alteraciones (Bernet) desde prevalencia×severidad (plantilla incluida).|Bernet et al. (1999)|
AG-090|H|Índices bióticos y diversidad|Intermedio|Shannon de morfotipos de MP y guía de aplicación AMBI/BMWP-Col con datos bentónicos.|Roldán (2003); Borja et al. (2000)|
AG-091|I|Índices acuáticos satelitales|Intermedio|NDWI/MNDWI/NDTI/NDCI sobre escena (raster propio o sintética de demostración).|McFeeters (1996); Xu (2006)|
AG-092|I|Detección de plásticos flotantes (FDI)|Doctoral|FDI de Biermann + umbral P98.5 y máscara de candidatos sobre agua.|Biermann et al. (2020)|
AG-093|I|Interpolación IDW/Kriging de MP|Doctoral|Superficie continua de concentración con variograma ajustado.|Isaaks & Srivastava (1989)|
AG-094|I|Autocorrelación espacial (Moran/LISA)|Doctoral|I global por permutación y clases HH/LL por estación.|Moran (1950); Anselin (1995)|
AG-095|I|Dispersión lagrangiana de MP|Doctoral|Monte Carlo de trayectorias desde la fuente principal (viento+corriente+Stokes).|Kowalski et al. (2016)|
AG-096|I|Búsqueda Sentinel-2 (STAC)|Intermedio|Escenas reales L2A de Planetary Computer sobre el área de estudio (o fallback).|Microsoft Planetary Computer|
AG-097|J|Cocientes de riesgo HQ|Doctoral|PEC/PNEC por polímero con clases y tabla editable de PNEC.|EU TGD-R.10 (2003); Frias & Nash (2019)|
AG-098|J|Exposición humana (One Health)|Doctoral|Monte Carlo de ingestión por consumo de pesca con sensibilidad de entradas.|WHO (2019); EFSA (2016)|
AG-099|J|Percepción social|Intermedio|Encuestas Likert: alfa de Cronbach, descriptivos y comparación por grupo.|Likert (1932)|
AG-100|J|Motor de recomendaciones de gestión|Doctoral|Síntesis priorizada de medidas según umbrales incumplidos (reglas trazables).|Elaboración OMNIPLAST (trazable a umbral)|
"""


def _parse():
    agents = {}
    for line in _SPEC.strip().splitlines():
        parts = line.split("|")
        agents[parts[0]] = {
            "id": parts[0], "categoria": parts[1],
            "nombre": parts[2], "nivel": parts[3], "desc": parts[4],
            "refs": [r for r in parts[5].split(";") if r.strip()],
        }
    return agents


REGISTRY = _parse()


def by_category(cat):
    return [a for a in REGISTRY.values() if a["categoria"] == cat]


def search(query):
    q = query.lower()
    return [a for a in REGISTRY.values()
            if q in a["nombre"].lower() or q in a["desc"].lower()
            or q in CATEGORIES[a["categoria"]][0].lower()]
