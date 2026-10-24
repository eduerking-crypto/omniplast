import sys
from pathlib import Path

_here = Path(__file__).resolve()
for _p in _here.parents:
    if (_p / "app" / "ui" / "bootstrap.py").exists():
        for _s in (str(_p), str(_p / "app")):
            if _s not in sys.path:
                sys.path.insert(0, _s)
        break

import streamlit as st
from ui.bootstrap import footer, kpi_row, page_header, section_title

from omniplast.agents import REGISTRY
from omniplast.db import dao
from omniplast.utils import fmt_num

st.markdown("""
<div class="omp-hero">
  <h1>OMNIPLAST <span class="accent">·</span> Observatorio Omni-Ecológico de Plásticos</h1>
  <p>Plataforma doctoral para <b>detectar, monitorear, asociar y comparar</b> microplásticos
  en cualquier ecosistema: estadística inferencial y bayesiana, multivariante comunitaria,
  qPCR y genes de degradación, teledetección satelital, geoestadística, riesgo ecológico y
  One Health — integrado en un solo sistema con trazabilidad científica.</p>
  <span class="omp-chip">100 agentes especialistas</span>
  <span class="omp-chip">3 matrices (agua · sedimento · biota)</span>
  <span class="omp-chip">qPCR / 16S / plastísfera</span>
  <span class="omp-chip">PERMANOVA · NMDS · GLM NB · mixtos · bayes</span>
  <span class="omp-chip">Sentinel-2 · kriging · trayectorias</span>
  <span class="omp-chip">API REST + scripts RStudio</span>
</div>
""", unsafe_allow_html=True)

try:
    summary = dao.project_summary()
except Exception:
    summary = {}

kpi_row([
    (fmt_num(summary.get("sitios", 0), 0), "Estaciones", "georreferenciadas", "#00A896"),
    (fmt_num(summary.get("muestras", 0), 0), "Muestras", "4 campañas estacionales", "#0A2540"),
    (fmt_num(summary.get("mp_registros", 0), 0), "Registros de MP", "polímero × morfología × tamaño", "#05BFDB"),
    (fmt_num(summary.get("geneticos", 0), 0), "Datos moleculares", "qPCR + metataxonomía 16S", "#E9B44C"),
    (fmt_num(len(REGISTRY), 0), "Agentes", "listos para ejecutar", "#9B5DE5"),
])

page_header("Un laboratorio doctoral completo", "Todo lo que el panel izquierdo ofrece está respaldado por funciones reales ejecutándose sobre su base de datos.")

left, right = st.columns(2, gap="large")
with left:
    st.markdown('<div class="section-title-container">', unsafe_allow_html=True)
    section_title("Investigación por área")
    cards = [
        ("Química & fisicoquímica", "WQI-NSF, estado trófico Carlson, DBO/DQO, nutrientes, salinidad, oxígeno con fórmulas de saturación publicadas."),
        ("Genética & molecular", "ΔΔCt con guías MIQE, curvas estándar y fusión, genotipado con HWE, panel petase/alkB/laccasa/phaC, filogenética."),
        ("Microbioma & plastísfera", "Diversidad alfa/beta, PERMANOVA + betadisper, abundancia diferencial tipo DESeq2, redes de co-ocurrencia, LEfSe."),
        ("Biota & ecotoxicología", "Ingestión por especie, biomarcadores SOD/CAT/GST/LPO, IBR integrado (Beliaeff-Burgeot), índices de condición."),
    ]
    for titulo, texto in cards:
        st.markdown(f'<div class="omp-card tall"><h4>{titulo}</h4><p>{texto}</p></div>',
                    unsafe_allow_html=True)
with right:
    section_title("Ecosistema & contexto")
    cards = [
        ("Teledetección", "Índices NDWI/MNDWI/NDTI/NDCI y detección de plásticos flotantes por FDI (Biermann 2020) con escenas reales o sintéticas."),
        ("Geoestadística & riesgo", "Kriging con variograma, Moran/LISA, hotspots Gi*, dispersión lagrangiana, HQ PEC/PNEC, SSD con HC5."),
        ("One Health & social", "Exposición humana por consumo de pesca (Monte Carlo), percepción social Likert con alfa de Cronbach."),
        ("Estadística doctoral", "Desde t/Welch hasta GLMM, GEE, cuantílica, mediación bootstrap, MCMC bayesiano y potencia a priori."),
    ]
    for titulo, texto in cards:
        st.markdown(f'<div class="omp-card tall"><h4>{titulo}</h4><p>{texto}</p></div>',
                    unsafe_allow_html=True)

section_title("Cómo empezar (3 pasos)")
st.markdown("""
<div class="omp-card tall">
<p><b>1 · Explore el Dashboard</b> — el proyecto demo (Ciénaga de La Virgen, <b>datos sintéticos verosímiles</b>)
ya está sembrado con 10 estaciones, 4 campañas y 3 matrices.</p>
<p><b>2 · Use sus propios datos</b> — cargue CSVs con la misma estructura de
<code>data/demo/*.csv</code> (plantillas exportadas) o conéctese a las bases públicas
(OBIS, GBIF, iNaturalist, Crossref, PubMed…).</p>
<p><b>3 · Delegue en los agentes</b> — la página <i>Agentes</i> ejecuta cualquiera de los
100 análisis con un clic; el <i>Asistente PLASTI</i> responde preguntas metodológicas
citando fuentes primarias.</p>
</div>
""", unsafe_allow_html=True)

section_title("Honestidad científica por diseño")
st.markdown("""
<div class="omp-card tall">
<p>· Los <b>umbrales</b> (PLI, PHI, WQI, HQ, TSI) provienen de literatura citada y son <b>editables</b>;
donde no existe norma oficial (microplásticos ambientales), OMNIPLAST lo dice explícitamente.</p>
<p>· No se inventan <b>DOIs, cebadores de genes funcionales ni valores normativos</b>: cuando un insumo es
taxon-específico, el sistema entrega la guía y la plantilla exactas.</p>
<p>· Los datos demo están etiquetados como <b>DATOS SINTÉTICOS</b> en cada informe exportado.</p>
</div>
""", unsafe_allow_html=True)

footer()
