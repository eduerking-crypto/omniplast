# OMNIPLAST — BRIEF MAESTRO (v1.0)
### Planta de requisitos expandida: lo pedido + lo que faltaba

> Este documento es el "prompt maestro" expandido de la plataforma: la visión original del
> investigador aumentada con las capas que un comité doctoral exigiría. Sirve como contrato
> de alcance y hoja de ruta de versión 2.

---

## 1. Visión

Construir el sistema de referencia en español para investigar **microplásticos en cualquier
ecosistema** (lagunas costeras, manglares, ríos, suelos, acuicultura, océano abierto, urbano):
desde la toma de muestra hasta el informe publicable, con **trazabilidad total** (cada número
puede rastrearse a un método citado), **cero datos inventados** y capacidad de **comparar
estudios** bajo unidad y protocolo comunes.

## 2. Ecosistema de datos (lo pedido, formalizado)

- **Multiparámetro**: fisicoquímica (17 variables), bioquímica (biomarcadores), biología
  (biota, ingestión), molecular (qPCR, 16S), social (encuestas), remoto (satélite).
- **Profundidad** (0/1.5/3 m), **geolocalización** (lat/lon por estación), **temporada**
  (secas/lluvias/transición), **tipo de ecosistema** (14 clases configurables),
  **exposición** (zona impacto/transición/referencia), **QA/QC por muestra** (blancos,
  duplicados, método de identificación).
- **Casos guía**: calidad de la Ciénaga de La Virgen por carga de MP; bacterias
  plastísferas y sus genes de degradación asociados al ambiente.

## 3. Capacidades analíticas (matriz completa)

| Capa | Implementado en OMNIPLAST 1.0 |
|---|---|
| Descriptiva + supuestos | Shapiro-Wilk/AD, Levene/Bartlett, transformaciones Box-Cox, atípicos IQR |
| Comparaciones | t/Welch (d de Cohen), ANOVA+Tukey (η²/ω²), Kruskal+Dunn (Holm), Games-Howell, Wilcoxon/Friedman, χ² (V de Cramer) |
| Asociaciones | Pearson/Spearman/Kendall con IC, correlación parcial, matrices con p |
| Regresiones | OLS (HC3, VIF), cuantílica, GLM Poisson/NB (conteos de MP), logística con OR |
| Jerárquicos | Mixtos (intercepto aleatorio, ICC), GEE intercambiable |
| Multivariante | PCA (+Varimax), FA, CA, RDA con pseudo-F por permutación, PCoA, NMDS-SMACOF, PERMANOVA+betadisper, ANOSIM, SIMPER, Procrustes |
| Series temporales | Mann-Kendall (+estacional Hirsch), pendiente de Sen, STL, ACF/PACF, Holt amortiguado |
| Bayes | Conjugado normal, beta-binomial, MCMC Metropolis de 2 medias, BF por BIC |
| ML | RF (clasif/reg + permutación), SVM, Gradient Boosting, k-means/silueta, DBSCAN, anomalías por PCA, ROC/PR+AUC, bootstrap de métricas |
| Genética | ΔΔCt Livak, curva estándar MIQE, Tm calling, genotipado+HWE, Fst, SNP-asociación, panel gen-ambiente |
| Microbioma | Alfa (S, Chao1, ACE, Shannon, Simpson, Pielou), beta Bray/Jaccard/Aitchison, DESeq2-lite, LEfSe-lite, core, redes (hubs, modularidad), UPGMA |
| Eco-tox | Ingestión, biomarcadores, IBR/2 (Beliaeff), histopatología Bernet, CF/HSI/GSI, AMBI/BMWP-Col |
| Riesgo | HQ por polímero, SSD lognormal + HC5 + PAF, msPAF, exposición pesca/agua Monte Carlo + tornado |
| Geo | IDW, kriging con variograma, Moran I + LISA, Gi*, drift lagrangiano con Stokes |
| Remoto | NDWI/MNDWI/NDTI/NDCI/NDVI/NDMI, FDI flotantes, STAC Planetary Computer, escena sintética demo |
| Síntesis | IICAM compuesto + Monte Carlo, apportionment de fuentes, motor de recomendaciones trazable |

## 4. Lo añadido que el brief original no pedía (valor doctoral)

1. **Cadena de custodia digital**: cada agente bitacora estado/duración/resumen en BD.
2. **Plantillas Excel oficiales** con validaciones desplegables y catálogos (LEEME incluido)
   + ejemplo lleno: garantiza datos comparables entre estudios (el mayor déficit de la
   literatura MP según Koelmans 2019).
3. **Asesor de decisión estadística** (árbol replicable) y **diseño muestral a priori**
   con potencia — evita el error más común: ANOVA sobre conteos sobredispersos.
4. **Scripts R canónicos** (vegan/DESeq2/lavaan/glmmTMB) para reproducibilidad en RStudio.
5. **Regla anti-invención**: sin DOI/cebador/valor normativo fabricado; lo que falta se
   declara con plantilla y protocolo exactos.
6. **API REST + Swagger** para integrarse con laboratorios (LIMS) y tesis colaborativas.
7. **QA automatizado de UI** (Playwright): 22 páginas sin excepciones; tests pytest del
   núcleo: 17/17.
8. **Modo degradado elegante**: sin internet, los conectores explican qué falta; sin
   Ollama, PLASTI responde con reglas+KB.

## 5. Versión 2.0 (roadmap)

- Importador real de plantillas Excel → SQLite (módulo `importador.py` + validador de
  claves foráneas) — la estructura ya está lista.
- Autenticación multiusuario + roles (investigador/técnico/revisor) sobre FastAPI.
- Metabarcoding real: carga de FASTQ + DADA2 vía subprocess opcional.
- CCA completa en Python y capscale (hoy vía script R generado).
- Serie temporal Sentinel-2 automatizada (stackstac) y detección de cambios de manglar
  conectada a GMW.
- Panel de revisión por pares internos (aprobación de resultados antes de informe).
- Publicación de datasets FAIR (DOI con Zenodo desde la página Reportes).
- Ollama con dominio específico fino-tuneado con papers del proyecto (RAG sobre
  biblioteca local).

## 6. Criterios de aceptación doctorales (checklist)

- [x] Todo número con cita o supuestos explícitos
- [x] Supuestos de cada prueba verificados antes de recomendarla
- [x] Incertidumbre reportada (IC bootstrap, Monte Carlo, bandas)
- [x] Datos demo etiquetados como sintéticos en cada exporte
- [x] Reproducibilidad: semilla fija, tests, scripts R, API
- [ ] Importador bidireccional completo (v2)
- [ ] Validación de IICAM con panel de expertos (documentar en METODOLOGIA.md)
