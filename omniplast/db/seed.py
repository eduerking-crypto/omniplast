"""Genera el proyecto demo 'Cienaga de La Virgen' con datos simulados verosimiles.

Dataset con fines academicos: reproduce patrones reportados en literatura para lagunas
costeras tropicales hiperurbanizadas (ver config.py). Todos los valores son sinteticos
y deterministas (rng seed 42). La UI lo etiqueta siempre como 'DATOS DEMO'.
"""
from datetime import datetime, timedelta

import numpy as np
import pandas as pd

from ..config import DEFAULT_STATIONS, DEMO_REGION
from .dao import get_engine, get_session
from .models import (BioParam, Campaign, EnvParam, GeneticResult, GeneticRun,
                     MPRecord, Project, Sample, Site)

RNG = np.random.default_rng(42)

CAMPAIGNS = [
    ("CAMP-24S", "Campaña Secas 2024", "secas", datetime(2024, 2, 12), datetime(2024, 2, 16)),
    ("CAMP-24L", "Campaña Lluvias 2024", "lluvias", datetime(2024, 10, 21), datetime(2024, 10, 25)),
    ("CAMP-25S", "Campaña Secas 2025", "secas", datetime(2025, 2, 10), datetime(2025, 2, 14)),
    ("CAMP-25L", "Campaña Lluvias 2025", "lluvias", datetime(2025, 10, 13), datetime(2025, 10, 17)),
]

WATER_PARAMS = {
    "impacto": {"temperatura_c": (29.8, .8), "ph": (7.2, .3), "od_mg_l": (2.8, 1.1),
                "dbo_mg_l": (68, 26), "dqo_mg_l": (160, 55), "turbidez_ntu": (95, 38),
                "salinidad_psu": (9, 4), "conductividad_us": (14000, 3500), "tss_mg_l": (85, 32),
                "nitrato_mg_l": (1.8, .7), "nitrito_mg_l": (.18, .08), "fosfato_mg_l": (.85, .35),
                "amonio_mg_l": (1.6, .7), "clorofila_a_ug_l": (28, 11), "toc_mg_l": (9.5, 2.8),
                "coliformes_totales_nmp": (46000, 15000), "e_coli_nmp": (12000, 5000)},
    "transicion": {"temperatura_c": (29.5, .7), "ph": (7.5, .25), "od_mg_l": (4.2, .9),
                   "dbo_mg_l": (32, 12), "dqo_mg_l": (90, 30), "turbidez_ntu": (55, 20),
                   "salinidad_psu": (14, 5), "conductividad_us": (21000, 4000), "tss_mg_l": (55, 18),
                   "nitrato_mg_l": (1.0, .4), "nitrito_mg_l": (.08, .03), "fosfato_mg_l": (.35, .15),
                   "amonio_mg_l": (.7, .3), "clorofila_a_ug_l": (15, 6), "toc_mg_l": (6.0, 1.8),
                   "coliformes_totales_nmp": (9000, 3500), "e_coli_nmp": (1800, 900)},
    "referencia_interna": {"temperatura_c": (29.3, .6), "ph": (7.8, .2), "od_mg_l": (5.5, .8),
                           "dbo_mg_l": (14, 5), "dqo_mg_l": (45, 15), "turbidez_ntu": (32, 12),
                           "salinidad_psu": (20, 6), "conductividad_us": (28000, 4500), "tss_mg_l": (35, 12),
                           "nitrato_mg_l": (.5, .2), "nitrito_mg_l": (.04, .02), "fosfato_mg_l": (.12, .05),
                           "amonio_mg_l": (.3, .15), "clorofila_a_ug_l": (7, 3), "toc_mg_l": (3.5, 1.1),
                           "coliformes_totales_nmp": (900, 400), "e_coli_nmp": (120, 60)},
    "referencia": {"temperatura_c": (28.9, .5), "ph": (8.1, .12), "od_mg_l": (6.3, .5),
                   "dbo_mg_l": (6, 2), "dqo_mg_l": (22, 7), "turbidez_ntu": (12, 4),
                   "salinidad_psu": (33, 2), "conductividad_us": (48000, 3000), "tss_mg_l": (14, 5),
                   "nitrato_mg_l": (.15, .08), "nitrito_mg_l": (.01, .005), "fosfato_mg_l": (.03, .015),
                   "amonio_mg_l": (.08, .04), "clorofila_a_ug_l": (2.2, 1), "toc_mg_l": (1.6, .5),
                   "coliformes_totales_nmp": (60, 30), "e_coli_nmp": (8, 5)},
}

SED_PARAMS = {
    "impacto": {"toc_pct": (7.8, 1.8), "arena_pct": (38, 10), "limo_arcilla_pct": (62, 10)},
    "transicion": {"toc_pct": (5.2, 1.4), "arena_pct": (52, 9), "limo_arcilla_pct": (48, 9)},
    "referencia_interna": {"toc_pct": (3.4, 1.0), "arena_pct": (66, 8), "limo_arcilla_pct": (34, 8)},
    "referencia": {"toc_pct": (1.2, .4), "arena_pct": (84, 5), "limo_arcilla_pct": (16, 5)},
}

MP_TOTALS = {
    "agua": {"impacto": (620, 1500), "transicion": (250, 700),
             "referencia_interna": (120, 350), "referencia": (40, 120)},
    "sedimento": {"impacto": (1800, 3400), "transicion": (900, 1600),
                  "referencia_interna": (400, 900), "referencia": (150, 400)},
    "biota": {"impacto": (8, 22), "transicion": (4, 9),
              "referencia_interna": (2, 5), "referencia": (1, 3)},
}

POLYMER_W = {
    "agua": {"PE": .28, "PP": .20, "PET": .18, "PS": .07, "PA": .07, "PVC": .04,
             "EPS": .06, "PMMA": .02, "UNKNOWN": .08},
    "sedimento": {"PE": .22, "PP": .16, "PET": .16, "PS": .07, "PA": .06, "PVC": .09,
                  "EPS": .04, "PMMA": .03, "PLA": .01, "UNKNOWN": .16},
    "biota": {"PE": .20, "PP": .14, "PET": .30, "PS": .06, "PA": .12, "PVC": .05,
              "EPS": .03, "PMMA": .02, "UNKNOWN": .08},
}

MORPH_BY_POLYMER = {
    "PET": [(.75, "fibra"), (.20, "fragmento"), (.05, "film")],
    "PA": [(.80, "fibra"), (.15, "fragmento"), (.05, "film")],
    "PE": [(.42, "fragmento"), (.38, "film"), (.15, "fibra"), (.05, "espuma")],
    "PP": [(.45, "fragmento"), (.35, "film"), (.15, "fibra"), (.05, "espuma")],
    "PS": [(.35, "fragmento"), (.30, "espuma"), (.35, "esfera")],
    "EPS": [(.85, "espuma"), (.15, "fragmento")],
    "PVC": [(.80, "fragmento"), (.20, "film")],
    "PMMA": [(.70, "esfera"), (.30, "fragmento")],
    "PLA": [(.60, "fragmento"), (.40, "film")],
    "UNKNOWN": [(.40, "fragmento"), (.35, "fibra"), (.15, "film"), (.10, "espuma")],
}
COLORS_BY_MORPH = {
    "fibra": ["azul", "transparente", "negro", "rojo", "azul", "verde"],
    "fragmento": ["blanco", "azul", "transparente", "verde", "rojo", "negro"],
    "film": ["transparente", "blanco", "azul", "negro"],
    "espuma": ["blanco", "amarillo", "naranja"],
    "gránulo": ["blanco", "azul", "negro"],
    "esfera": ["transparente", "blanco", "azul"],
}
SIZE_BANDS = [(1, 100), (100, 300), (300, 500), (500, 1000), (1000, 5000)]
SIZE_W_WATER = [.06, .22, .26, .24, .22]
SIZE_W_SED = [.03, .18, .27, .28, .24]
SIZE_W_BIOTA = [.08, .34, .30, .20, .08]

GENES = ["petase", "mhetase", "laccasa", "alkb", "phac", "16s"]
GENE_ZONE_CENTERS = {
    "petase": {"impacto": 29.5, "transicion": 31.5, "referencia_interna": 33.5, "referencia": 35.0},
    "mhetase": {"impacto": 31.5, "transicion": 33.5, "referencia_interna": 35.5, "referencia": 37.0},
    "laccasa": {"impacto": 25.5, "transicion": 27.5, "referencia_interna": 29.5, "referencia": 31.0},
    "alkb": {"impacto": 28.5, "transicion": 30.5, "referencia_interna": 32.5, "referencia": 34.0},
    "phac": {"impacto": 26.5, "transicion": 28.5, "referencia_interna": 30.5, "referencia": 32.0},
    "16s": {"impacto": 15.8, "transicion": 15.8, "referencia_interna": 15.6, "referencia": 16.2},
}
GENE_DETECT_P = {"petase": .60, "mhetase": .35, "laccasa": .78, "alkb": .50,
                 "phac": .70, "16s": 1.0}

GENERA = [
    ("Pseudomonas", 900, 2.4), ("Vibrio", 850, 2.2), ("Bacillus", 800, 1.6),
    ("Rhodococcus", 420, 2.0), ("Stenotrophomonas", 380, 2.1), ("Acinetobacter", 500, 1.8),
    ("Micrococcus", 300, 1.4), ("Alcanivorax", 140, 2.6), ("Marinobacter", 260, 1.7),
    ("Shewanella", 240, 1.9), ("Flavobacterium", 520, 1.2), ("Cytophaga", 300, 1.1),
    ("Tenacibaculum", 160, 1.6), ("Oleispira", 90, 2.2), ("Thalassospira", 150, 1.8),
    ("Photobacterium", 320, 1.5), ("Aeromonas", 260, 1.9), ("Enterococcus", 300, 2.0),
    ("Escherichia/Shigella", 350, 2.1), ("Bacteroides", 260, 1.4), ("Nitrosomonas", 180, 1.0),
    ("Nitrospira", 200, 1.0), ("Planctomyces", 150, 1.1), ("Roseobacter-clade", 280, 1.2),
]
MP_ENRICHED = {"Pseudomonas", "Vibrio", "Bacillus", "Rhodococcus", "Stenotrophomonas",
               "Acinetobacter", "Micrococcus", "Alcanivorax", "Marinobacter", "Shewanella",
               "Oleispira", "Thalassospira", "Aeromonas"}

BIOTA_SPECIES = {
    "ostra": ("Crassostrea rhizophorae", "Ostra de manglar"),
    "pez": ("Mugil incilis", "Chula / lebranche (mugílido)"),
}


def _zone_factor(zona):
    return {"impacto": 1.0, "transicion": .66, "referencia_interna": .45, "referencia": .22}[zona]


def _draw_morph(polymer):
    opts = MORPH_BY_POLYMER.get(polymer, MORPH_BY_POLYMER["UNKNOWN"])
    u = RNG.random()
    acc = 0.0
    for p, morph in opts:
        acc += p
        if u <= acc:
            return morph
    return opts[-1][1]


def _mp_records_for(sample_id, total, unidad, pol_w, size_w, method):
    n_combos = int(RNG.integers(7, 13))
    polymers = list(pol_w.keys())
    probs = np.array([pol_w[p] for p in polymers])
    probs = probs / probs.sum()
    rows = []
    if total < n_combos:
        n_combos = max(3, int(total))
    for _ in range(n_combos):
        polymer = RNG.choice(polymers, p=probs)
        morph = _draw_morph(polymer)
        color = RNG.choice(COLORS_BY_MORPH[morph])
        lo, hi = SIZE_BANDS[int(RNG.choice(len(SIZE_BANDS), p=size_w))]
        rows.append({
            "sample_id": sample_id, "polimero": polymer, "morfologia": morph, "color": color,
            "tamano_min_um": float(lo), "tamano_max_um": float(hi),
            "abundancia": 0.0, "unidad": unidad, "metodo_id": method, "peso_mg": 0.0,
        })
    fracs = RNG.dirichlet(np.ones(n_combos) * 3.2)
    counts = np.maximum(1, np.round(total * fracs)).astype(float)
    for row, c in zip(rows, counts):
        row["abundancia"] = float(c)
        row["peso_mg"] = float(c * RNG.uniform(.0006, .0021))
    return [MPRecord(**r) for r in rows]


def seed_demo_if_empty():
    engine = get_engine()
    from .models import create_all
    create_all(engine)
    with get_session() as s:
        if s.query(Sample).count() > 0:
            return False
    _seed()
    return True


def _seed():
    with get_session() as s:
        proj = Project(codigo="CLV-DEMO", nombre="Microplásticos en la Ciénaga de La Virgen",
                       descripcion=("Proyecto doctoral demo: carga, fuentes, plastísfera y efectos "
                                    "ecotoxicológicos de microplásticos en un humedal costero "
                                    "hiperurbanizado del Caribe colombiano. DATOS SIMULADOS."),
                       ecosistema="Laguna costera", region=DEMO_REGION,
                       responsable="Equipo OMNIPLAST", fecha_inicio=datetime(2024, 2, 1),
                       fecha_fin=datetime(2025, 11, 30))
        s.add(proj)
        s.flush()
        sites = {}
        for code, name, lat, lon, zone in DEFAULT_STATIONS:
            st = Site(project_id=proj.id, codigo=code, nombre=name, lat=lat, lon=lon,
                      zona=zone, descripcion=f"Estación {code} — zona {zone}")
            s.add(st)
            sites[code] = st
        camps = {}
        for code, name, season, d0, d1 in CAMPAIGNS:
            c = Campaign(project_id=proj.id, codigo=code, nombre=name, temporada=season,
                         fecha_inicio=d0, fecha_fin=d1,
                         objetivo="Muestreo multimatriz con perfiles de profundidad y QA/QC completo")
            s.add(c)
            camps[code] = c
        s.flush()

        for code, c in camps.items():
            season = c.temporada
            year = c.fecha_inicio.year
            for st in sites.values():
                zone = st.zona
                depth_set = [0.0]
                if st.codigo in ("MED-01", "MED-02", "NOR-01", "JPC-01"):
                    depth_set = [0.0, 1.5, 3.0]
                for depth in depth_set:
                    smp = Sample(
                        campaign_id=c.id, site_id=st.id,
                        codigo=f"{st.codigo}-{code}-W{int(depth*10):02d}",
                        fecha=c.fecha_inicio + timedelta(days=int(RNG.integers(0, 5))),
                        matriz="agua", subtipo="superficie" if depth == 0 else "columna",
                        profundidad_m=depth, replicas=3,
                        metodo="Red de neustón 330 um" if depth == 0 else "Niskin",
                        colector="Equipo de campo CLV",
                        qa_qc={"blanco_proceso": True, "duplicado": True},
                        notas="")
                    s.add(smp)
                    s.flush()
                    season_f = 1.4 if season == "lluvias" else 1.0
                    year_f = 1.12 if year == 2025 else 1.0
                    depth_f = {0.0: 1.0, 1.5: .82, 3.0: .66}[depth]
                    base = MP_TOTALS["agua"][zone]
                    total = float(RNG.uniform(*base) * season_f * year_f * depth_f)
                    s.add_all(_mp_records_for(smp.id, total, "p/m3", POLYMER_W["agua"],
                                             SIZE_W_WATER, "FTIR microscopio (FPA)"))
                    for param, (mu, sd) in WATER_PARAMS[zone].items():
                        val = RNG.normal(mu, sd)
                        if season == "lluvias":
                            val *= {"temperatura_c": .97, "turbidez_ntu": 1.25, "tss_mg_l": 1.3,
                                    "nitrato_mg_l": 1.4, "salinidad_psu": .6, "dbo_mg_l": 1.15,
                                    "conductividad_us": .75}.get(param, 1.0)
                        if depth > 0:
                            val *= {"od_mg_l": 1 - depth * .11, "clorofila_a_ug_l": .55,
                                    "turbidez_ntu": 1.1}.get(param, 1.0)
                        if param in ("coliformes_totales_nmp", "e_coli_nmp"):
                            val = max(2, val)
                        s.add(EnvParam(sample_id=smp.id, parametro=param, valor=round(float(val), 4),
                                       unidad="", detectado=True))

                smp = Sample(campaign_id=c.id, site_id=st.id,
                             codigo=f"{st.codigo}-{code}-SED",
                             fecha=c.fecha_inicio + timedelta(days=int(RNG.integers(0, 5))),
                             matriz="sedimento", subtipo="sedimento superficial (0-5 cm)",
                             profundidad_m=0.0, replicas=3, metodo="Nucleo de sedimento",
                             colector="Equipo de campo CLV", qa_qc={"blanco_proceso": True}, notas="")
                s.add(smp)
                s.flush()
                season_f = 1.18 if season == "lluvias" else 1.0
                year_f = 1.10 if year == 2025 else 1.0
                base = MP_TOTALS["sedimento"][zone]
                total = float(RNG.uniform(*base) * season_f * year_f)
                s.add_all(_mp_records_for(smp.id, total, "p/kg dw", POLYMER_W["sedimento"],
                                         SIZE_W_SED, "FTIR-ATR"))
                for param, (mu, sd) in SED_PARAMS[zone].items():
                    val = RNG.normal(mu, sd)
                    if param == "arena_pct":
                        val = np.clip(val, 5, 98)
                    if param == "limo_arcilla_pct":
                        val = 100 - np.clip(val, 2, 95)
                    s.add(EnvParam(sample_id=smp.id, parametro=param, valor=round(float(val), 3),
                                   unidad="", detectado=True))

                if st.codigo in ("TUR-01", "TUR-02", "VIR-01", "MED-01"):
                    for sp_key, (sci, _) in BIOTA_SPECIES.items():
                        bsmp = Sample(campaign_id=c.id, site_id=st.id,
                                      codigo=f"{st.codigo}-{code}-BIO-{sp_key.upper()}",
                                      fecha=c.fecha_inicio + timedelta(days=int(RNG.integers(0, 5))),
                                      matriz="biota", subtipo=sci,
                                      profundidad_m=0.0, replicas=10, metodo="Dissección biota",
                                      colector="Equipo de campo CLV",
                                      qa_qc={"blanco_proceso": True, "n_individuos": 10}, notas="")
                        s.add(bsmp)
                        s.flush()
                        base = MP_TOTALS["biota"][zone]
                        total = float(RNG.uniform(*base) * (1.15 if sp_key == "ostra" else 1.0))
                        s.add_all(_mp_records_for(bsmp.id, total, "p/ind", POLYMER_W["biota"],
                                                  SIZE_W_BIOTA, "UV-visual + tinción (Nile Red)"))
                        zone_f = _zone_factor(zone)
                        biomarkers = {
                            "sod_u_mg_prot": RNG.uniform(45, 95) * (1.9 - zone_f),
                            "cat_u_mg_prot": RNG.uniform(30, 70) * (1.8 - zone_f),
                            "gst_u_mg_prot": RNG.uniform(20, 48) * (1.7 - zone_f),
                            "lpo_nmol_tba_mg": RNG.uniform(2.0, 4.5) * (1.8 - zone_f),
                        }
                        for param, val in biomarkers.items():
                            s.add(BioParam(sample_id=bsmp.id, taxon=sci, tipo="biomarcador",
                                           parametro=param, valor=round(float(val), 3),
                                           unidad="", tejido="hepatopáncreas" if sp_key == "ostra" else "hígado"))
                        if sp_key == "pez":
                            s.add(BioParam(sample_id=bsmp.id, taxon=sci, tipo="indice",
                                           parametro="factor_condicion", valor=round(float(RNG.uniform(0.72, 1.28) * (0.6 + .5 * zone_f)), 3), unidad="", tejido=""))
                            s.add(BioParam(sample_id=bsmp.id, taxon=sci, tipo="indice",
                                           parametro="hsi_pct", valor=round(float(RNG.uniform(0.9, 2.6)), 3), unidad="", tejido=""))
                            s.add(BioParam(sample_id=bsmp.id, taxon=sci, tipo="indice",
                                           parametro="gsi_pct", valor=round(float(RNG.uniform(1.5, 8.0)), 3), unidad="", tejido=""))

                for target_sample in [x for x in (s.query(Sample).filter_by(site_id=st.id)
                                                  .filter(Sample.campaign_id == c.id).all())]:
                    if target_sample.matriz == "sedimento":
                        run = GeneticRun(sample_id=target_sample.id, tipo="qpcr",
                                         objetivo="Panel de genes de degradación (plastísfera)",
                                         plataforma="qPCR SYBR Green", protocolo="MIQE-compliant")
                        s.add(run)
                        s.flush()
                        for gene in GENES:
                            center = GENE_ZONE_CENTERS[gene][zone]
                            center -= (0.8 if season == "lluvias" else 0.0)
                            detected = RNG.random() < GENE_DETECT_P[gene]
                            if detected:
                                ct = float(RNG.normal(center, .55))
                                copies = float(10 ** (11 - 0.30103 * ct))
                                s.add(GeneticResult(run_id=run.id, sample_id=target_sample.id,
                                                    target=gene, metrica="copies_g",
                                                    valor=round(copies, 1), ct=round(ct, 2), unidad="copias/g",
                                                    detectado=True))
                            else:
                                s.add(GeneticResult(run_id=run.id, sample_id=target_sample.id,
                                                    target=gene, metrica="copies_g", valor=0.0,
                                                    ct=None, unidad="copias/g", detectado=False))
                        run16 = GeneticRun(sample_id=target_sample.id, tipo="16S",
                                           objetivo="Amplicón V4 515F-Y/806R-A (metabarcoding)",
                                           plataforma="Illumina MiSeq 2x250", protocolo="Earth Microbiome Project")
                        s.add(run16)
                        s.flush()
                        mp_load = _zone_factor(zone)
                        for genus, base_ab, enrich in GENERA:
                            ab = RNG.normal(base_ab, base_ab * .35)
                            if genus in MP_ENRICHED:
                                ab *= (1 + enrich * mp_load * RNG.uniform(.6, 1.1))
                            ab = max(0, ab)
                            s.add(GeneticResult(run_id=run16.id, sample_id=target_sample.id,
                                                target=genus, metrica="reads", valor=int(ab),
                                                unidad="lecturas", detectado=True))
        s.commit()
    _export_csvs()


def _export_csvs():
    from .dao import bio_df, env_df, mp_df, qpcr_df, samples_df, sites_df
    from ..config import DEMO_DIR
    exports = {
        "sitios.csv": sites_df(), "muestras.csv": samples_df(),
        "microplasticos.csv": mp_df(), "parametros_ambientales.csv": env_df(),
        "biota.csv": bio_df(), "qpcr.csv": qpcr_df(),
    }
    for name, df in exports.items():
        df.to_csv(DEMO_DIR / name, index=False, encoding="utf-8-sig")
