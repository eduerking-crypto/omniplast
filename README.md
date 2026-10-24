# 🌊 OMNIPLAST — Plataforma Doctoral de Análisis Integral de Microplásticos

> Detectar · Monitorear · Asociar · Comparar · Comunicar — microplásticos en **cualquier ecosistema**.

OMNIPLAST es un sistema de investigación de nivel doctoral que integra estadística inferencial y bayesiana, ecología molecular (qPCR/16S), teledetección satelital, geoestadística, evaluación de riesgo y One Health en una sola plataforma trazable, con **100 agentes especialistas** y un asistente científico (PLASTI) que cita fuentes primarias.

```
Python 3.11+ · Streamlit · FastAPI · SQLite · Plotly · Folium · scikit-learn · statsmodels · Biopython
```

---

## Inicio rápido

```bat
:: Windows (doble clic también sirve)
run.bat          → http://localhost:8787  (app + base demo sembrada automáticamente)
run_api.bat      → http://localhost:8010/docs  (API REST con Swagger)
```

```bash
# Linux/macOS
python -m streamlit run app/Home.py
python -m uvicorn api.main:app --port 8010
```

Primer arranque: la base SQLite se crea y se siembra sola con el **proyecto demo "Ciénaga de La Virgen" (Cartagena de Indias)** — 10 estaciones, 4 campañas, 3 matrices (agua, sedimento, biota), 1.248 registros de MP, panel qPCR de genes de degradación y metataxonomía 16S. Los datos son **sintéticos verosímiles**, etiquetados como DEMO en cada exporte.

## Los entregables clave

| Módulo | Qué hace | Dónde |
|---|---|---|
| **100 agentes** | Cada uno ejecuta un flujo analítico completo con interpretación y referencias | Página *Agentes* |
| **Motor estadístico** | t/Welch, ANOVA+Tukey, Kruskal+Dunn, Games-Howell, GLM Poisson/NB, mixtos, GEE, cuantílica, mediación bootstrap, PERMANOVA+betadisper, ANOSIM, SIMPER, PCA/RDA, NMDS (SMACOF), Mann-Kendall/Sen, MCMC, potencia | Páginas 06-10 |
| **Genética & molecular** | ΔΔCt (Livak), eficiencia MIQE, Tm calling, genotipado+HWE, asociación gen-ambiente, diversidad alfa/beta, DESeq2-lite, LEfSe, redes de co-ocurrencia | Páginas 11-12 |
| **Teledetección & GIS** | NDWI/MNDWI/NDTI/NDCI, detección de plásticos flotantes (FDI, Biermann 2020), kriging con variograma, Moran/LISA, hotspots Gi*, modelo lagrangiano de dispersión | Páginas 13-14 |
| **Riesgo & One Health** | HQ (PEC/PNEC), SSD con HC5, exposición humana Monte Carlo, percepción social | Página 16 |
| **Índices de calidad** | CF (Hakanson), PLI (Tomlinson), PHI (Lithner), WQI-NSF, TSI (Carlson) e **IICAM** — índice compuesto propuesto con incertidumbre Monte Carlo | Página 17 |
| **Bases públicas** | OBIS, GBIF, iNaturalist, Crossref, PubMed, NASA POWER, OpenAQ (fallan explícitamente sin red: cero registros inventados) | Páginas 15 y 21 |
| **PLASTI** | Asistente con base de conocimiento curada + Ollama local + asesor de decisión estadística y diseño muestral | Página 19 |
| **Reportes** | Informe DOCX académico, scripts R (vegan/DESeq2/lavaan/glmmTMB) para RStudio | Página 20 |

### Plantillas Excel oficiales
Descárguelas desde *Proyectos → Cargar datos propios* (o pág. 20):
- **PLANTILLA_OMNIPLAST.xlsx** — estructura vacía con validaciones desplegables, catálogo de parámetros con unidades y guía LEEME.
- **EJEMPLO_OMNIPLAST.xlsx** — la misma estructura **llena con el dataset demo** para reproducir cualquier análisis de inmediato.

## Honestidad científica por diseño

- Umbrales e índices **con cita primaria y editables** (`omniplast/config.py`); donde no existe norma oficial se dice explícitamente (p. ej., no hay valor regulatorio ambiental vigente para MP; la UE reguló *metodología* — Reglamento (UE) 2024/1448).
- **No se inventan DOIs ni cebadores de genes funcionales**: las secuencias de cebadores son solo las universales clásicas publicadas (515F-Y/806R-A, 27F/1492R, ITS1F/ITS4); el resto se diseña desde secuencias de referencia con el módulo de diseño.
- QA/QC de laboratorio embebido (Koelmans et al. 2019): blancos, duplicados, % UNKNOWN, método de identificación por registro.

## Arquitectura

```
MICROPLASTICO/
├── omniplast/               # núcleo científico (sin UI)
│   ├── config.py            # catálogos, umbrales citados, paleta, tokens
│   ├── db/                  # modelos SQLAlchemy + DAO + sembrado demo determinista
│   ├── core/
│   │   ├── statistics/      # descriptiva, inferencial, multivariante, series, ml, mixtos, bayes
│   │   ├── genetics/        # qPCR, comunidad 16S, cebadores, genotipado
│   │   ├── eco/             # biomarcadores/IBR, One Health, toxicología (HQ/SSD)
│   │   ├── geo/             # kriging/Moran/Gi*, índices satelitales/FDI, drift lagrangiano
│   │   ├── connectors/      # 8 bases públicas gratuitas con cache y fallback honesto
│   │   ├── indices.py       # PLI, PHI, WQI, TSI, AMBI, BMWP-Col, IICAM
│   │   ├── sources.py       # apportionment polímero×morfología
│   │   └── uncertainty.py   # Monte Carlo genérico + tornado
│   ├── agents/              # registro de 100 agentes + handlers + engine + bitácora
│   ├── assistant/           # PLASTI: KB curada + asesor estadístico + Ollama
│   └── reporting/           # figuras con estilo, DOCX, scripts R, plantillas Excel
├── app/                     # Streamlit: Home + 21 páginas + sistema de diseño
├── api/                     # FastAPI (Swagger en :8010/docs)
├── scripts/                 # smoke tests, QA de navegador, generación de plantillas
├── tests/                   # 17 pruebas pytest del núcleo
└── docs/                    # brief maestro, metodología con fórmulas y citas
```

## Verificación incluida

```bash
python -m pytest tests -q        # 17 pruebas del núcleo científico
python scripts/smoke_agents.py   # ejecuta los 100 agentes contra la base demo
python scripts/qa_navegador.py   # Playwright: 22 páginas, cero excepciones (con app viva)
```

Estado verificado en este repositorio: **17/17 tests · 100/100 agentes (97 ok, 3 requieren insumos externos por diseño) · 22/22 páginas sin excepciones · 0 errores JS**.

## Créditos normativos y científicos principales

Hartmann et al. (2019) · GESAMP (2019) · Koelmans et al. (2019) · Hakanson (1980) · Tomlinson et al. (1980) · Lithner et al. (2011) · Brown et al. (1970) · Carlson (1977) · Anderson (2001) · Beliaeff & Burgeot (2002) · Yoshida et al. (2016) · Zettler et al. (2013) · Bustin et al. (2009, MIQE) · Borja et al. (2000) · Roldán (2003, BMWP-Col) · Biermann et al. (2020) · OMS (2019) · Reglamento (UE) 2024/1448.

La lista completa con fórmulas está en `docs/METODOLOGIA.md`.
