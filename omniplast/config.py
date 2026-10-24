"""Configuracion central: rutas, catalagos cientificos, umbrales y tokens de disenio.

Referencias clave de los umbrales (no inventadas; verificar vigencia antes de publicar):
- Hakanson (1980) Limnol. Surv. — clases del factor de contaminacion CF.
- Tomlinson et al. (1980) Helgoland Mar. Res. — Pollution Load Index (PLI).
- Lithner et al. (2011) Environ. Sci. Technol. 45:4727-4735 — puntajes de riesgo quimico de polimeros.
- Frias & Nash (2019) Mar. Pollut. Bull. 138:316-322 — marco de riesgo relativo de MPs.
- Carlson (1977) — Trophic State Index.
- Borja et al. (2000) Mar. Pollut. Bull. 40:1100-1108 — clases AMBI.
- Roldan (2003) — BMWP adaptado para Colombia (BMWP-Col).
- Hartmann et al. (2019) Environ. Sci. Technol. 53:1039-1047 — definicion de tamanos (1 um - 5 mm).
- GESAMP (2019) Guidelines for monitoring plastic litter — protocolos de muestreo.
- Bustin et al. (2009) Clin. Chem. 55:611-622 — guias MIQE para qPCR.
- Reglamento (UE) 2024/1448 — metodologia de medida de microplasticos en agua de consumo.
"""
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data"
EXPORT_DIR = DATA_DIR / "exports"
CACHE_DIR = DATA_DIR / "cache"
DEMO_DIR = DATA_DIR / "demo"
DB_PATH = DATA_DIR / "omniplast.db"
for _d in (DATA_DIR, EXPORT_DIR, CACHE_DIR, DEMO_DIR):
    _d.mkdir(parents=True, exist_ok=True)

APP_NAME = "OMNIPLAST"
TAGLINE = "Plataforma Doctoral de Análisis Integral de Microplásticos"
VERSION = "1.0.0"
API_PORT = 8010
STREAMLIT_PORT = 8787

PALETTE = {
    "ink": "#0E2433",
    "navy": "#0A2540",
    "navy2": "#123A5C",
    "teal": "#00A896",
    "aqua": "#05BFDB",
    "coral": "#FF6B6B",
    "amber": "#E9B44C",
    "green": "#3E9B4F",
    "sand": "#F5F8FB",
    "surface": "#FFFFFF",
    "line": "#E1E8F0",
    "muted": "#5B7186",
}

CHART_SEQUENCE = ["#0A2540", "#00A896", "#05BFDB", "#E9B44C", "#FF6B6B",
                  "#3E9B4F", "#9B5DE5", "#4CC9F0", "#F07CAD", "#7C9EB2"]

MORPHOLOGY_COLORS = {
    "fibra": "#FF6B6B", "fragmento": "#05BFDB", "film": "#00A896",
    "espuma": "#E9B44C", "gránulo": "#9B5DE5", "esfera": "#4CC9F0",
}

QUALITY_CLASSES = [
    ("Excelente", 80, 100), ("Buena", 60, 79), ("Aceptable", 45, 59),
    ("Deteriorada", 30, 44), ("Crítica", 0, 29),
]
QUALITY_COLORS = {"Excelente": "#3E9B4F", "Buena": "#00A896", "Aceptable": "#E9B44C",
                  "Deteriorada": "#F08A24", "Crítica": "#FF6B6B"}

POLYMER_HAZARD = {
    "PE": 11, "PP": 11, "PS": 30, "EPS": 30, "PET": 51, "PC": 68,
    "PVC": 5000, "PA": 2056, "PU": 1156, "PUR": 1156, "PVDC": 5000,
    "ABS": 8552, "PMMA": 2500, "PLA": 4, "PVA": 1, "POM": 1000, "EVA": 11,
}
POLYMER_HAZARD_FLAGGED = {"PMMA", "PLA", "PVA", "POM", "EVA", "PVDC", "EPS"}
POLYMER_NAMES = {
    "PE": "Polietileno", "PP": "Polipropileno", "PS": "Poliestireno",
    "EPS": "Poliestireno expandido", "PET": "Polietilentereftalato (polyester)",
    "PC": "Policarbonato", "PVC": "Policloruro de vinilo", "PA": "Poliamida (nylon)",
    "PU": "Poliuretano", "PUR": "Poliuretano", "PVDC": "Policloruro de vinilideno",
    "ABS": "Acrilonitrilo butadieno estireno", "PMMA": "Polimetilmetacrilato",
    "PLA": "Acido polilactico", "PVA": "Alcohol polivinilico", "POM": "Poliacetal",
    "EVA": "Etilen vinil acetato", "UNKNOWN": "No identificado",
}
POLYMER_DENSITY = {
    "PE": 0.94, "PP": 0.91, "PS": 1.05, "EPS": 0.03, "PET": 1.38,
    "PC": 1.20, "PVC": 1.40, "PA": 1.14, "PU": 1.20, "PUR": 1.20,
    "PMMA": 1.18, "PLA": 1.25, "PVA": 1.19, "EVA": 0.95, "ABS": 1.06,
}

PHI_CLASSES = [
    ("Bajo", 0, 1), ("Moderado", 1, 10), ("Alto", 10, 100),
    ("Peligroso", 100, 1000), ("Extremadamente peligroso", 1000, float("inf")),
]

CF_CLASSES = [
    ("Baja", 0, 1), ("Moderada", 1, 3), ("Considerable", 3, 6), ("Muy alta", 6, float("inf")),
]

HQ_CLASSES = [
    ("Riesgo despreciable", 0, 0.1), ("Riesgo bajo", 0.1, 1),
    ("Riesgo medio", 1, 10), ("Riesgo alto", 10, float("inf")),
]

MP_TIERS = {
    "agua": [("Baja", 0, 100), ("Moderada", 100, 1000), ("Alta", 1000, 5000),
             ("Extremadamente alta", 5000, float("inf"))],
    "sedimento": [("Baja", 0, 500), ("Moderada", 500, 2000), ("Alta", 2000, 5000),
                  ("Muy alta", 5000, float("inf"))],
    "biota": [("Baja", 0, 1), ("Moderada", 1, 5), ("Alta", 5, 15),
              ("Muy alta", 15, float("inf"))],
}
MP_TIERS_NOTE = ("Clasificacion referencial construida a partir de rangos publicados en "
                 "estudios costeros/estuarinos tropicales (revisar Hartmann et al. 2019; "
                 "GESAMP 2019). No existe valor normativo internacional vigente para MPs "
                 "ambientales; el unico marco metodologico normativo de la UE es el "
                 "Reglamento (UE) 2024/1448 sobre medicion en agua de consumo humano.")

WQI_WEIGHTS = {
    "od_mg_l": 0.17, "coliformes_totales_nmp": 0.16, "ph": 0.11, "dbo_mg_l": 0.11,
    "temperatura_c": 0.10, "fosfato_mg_l": 0.10, "nitrito_mg_l": 0.10,
    "turbidez_ntu": 0.08, "tss_mg_l": 0.07,
}
WQI_CLASSES = [
    ("Excelente", 91, 100), ("Buena", 71, 90), ("Media", 51, 70),
    ("Mala", 26, 50), ("Muy mala", 0, 25),
]

TSI_CLASSES = [
    ("Oligotrófico", 0, 30), ("Oligo-mesotrófico", 30, 40), ("Mesotrófico", 40, 50),
    ("Meso-eutrófico", 50, 60), ("Eutrófico", 60, 70), ("Hipereutrófico", 70, 200),
]

AMBI_CLASSES = [
    ("No perturbado", 0.0, 1.2), ("Ligeramente perturbado", 1.2, 3.3),
    ("Moderadamente perturbado", 3.3, 4.3), ("Fuertemente perturbado", 4.3, 5.3),
    ("Muy fuertemente perturbado", 5.3, 6.0), ("Extremadamente perturbado", 6.0, 7.0),
]
AMBI_FORMULA_NOTE = "AMBI = (0*%GI + 1.5*%GII + 3*%GIII + 4.5*%GIV + 6*%GV) / 10 (Borja et al. 2000)."

BMWP_COL_CLASSES = [
    ("Muy buena calidad", 121, 500), ("Buena calidad", 101, 120),
    ("Calidad aceptable", 61, 100), ("Calidad dudosa", 36, 60),
    ("Calidad crítica", 16, 35), ("Muy crítica", 0, 15),
]

SIZE_CLASSES_UM = [
    ("mega", 20000, None), ("macro", 5000, 20000), ("micro", 1000, 5000),
    ("micro-pequeno", 500, 1000), ("micro-min", 100, 500), ("nano-micro", 1, 100),
]
MP_DEFINITION = "Microplastico: particulas de 1 um a 5 mm (Hartmann et al. 2019)."

MORPHOLOGIES = ["fibra", "fragmento", "film", "espuma", "gránulo", "esfera"]
MP_COLORS = ["azul", "transparente", "blanco", "negro", "rojo", "verde", "amarillo", "naranja", "rosa"]
MATRICES = ["agua", "sedimento", "biota", "aire"]
DEPTHS_M = [0.0, 1.5, 3.0]
ECOSYSTEM_TYPES = [
    "Laguna costera", "Ciénaga estuarina", "Manglar", "Arrecife coralino",
    "Océano abierto", "Río", "Lago", "Pantano freshwater", "Playa arenosa",
    "Suelo agrícola", "Urbano", "Acuícola", "Estuario", "Aguas subterráneas",
]
SEASONS = {"secas": "Secas (dic-abr)", "lluvias": "Lluvias (sep-nov)", "transicion": "Transición"}
SAMPLING_METHODS = ["Red de neustón 330 um", "Red de plancton 100 um", "Muestreo superficial (cuchara)",
                    "Niskin", "Nucleo de sedimento", "Cuadrante de sedimento", "Dissección biota",
                    "Trampa pasiva de depositación", "Bombas de alto volumen"]
ID_METHODS = ["FTIR-ATR", "FTIR microscopio (FPA)", "Raman", "UV-visual + tinción (Nile Red)",
              "Piro-GC/MS", "SEM-EDX", "Identificación visual (hot needle test)"]

GENE_CATALOG = {
    "petase": {"gen": "petase (ISF6_4831)", "organismo": "Ideonella sakaiensis 201-F6",
               "funcion": "Hidrolisis de PET a MHET",
               "ref": "Yoshida et al. (2016) Science 351:1196-1203"},
    "mhetase": {"gen": "mhetase (ISF6_0224)", "organismo": "Ideonella sakaiensis 201-F6",
                "funcion": "Hidrolisis de MHET a ácido tereftálico",
                "ref": "Yoshida et al. (2016) Science 351:1196-1203"},
    "laccasa": {"gen": "lcc (laccasa bacteriana)", "organismo": "Themonas/Bacillus/Pseudomonas spp.",
                "funcion": "Oxidacion de polimeros aromaticos y aditivos",
                "ref": "Danso, Chow & Streit (2019) Appl. Environ. Microbiol. 85:e01095-19"},
    "alkb": {"gen": "alkB (monooxigenasa de alcanos)", "organismo": "Alcanivorax, Rhodococcus spp.",
             "funcion": "Oxidacion de cadenas alifaticas (PE-like)",
             "ref": "van Beilen & Funhoff (2007) Appl. Microbiol. Biotechnol. 74:13-21"},
    "phac": {"gen": "phaC (sintasa de PHA)", "organismo": "Amplia (bacterias acumuladoras de PHA)",
             "funcion": "Sintesis de bioplasticos; indicador de plasticidad metabolica",
             "ref": "Tokiwa et al. (2009) Appl. Microbiol. Biotechnol. 80:153-174"},
    "est_lip": {"gen": "estA/lipA (esterasas/lipasas)", "organismo": "Pseudomonas, Bacillus spp.",
                "funcion": "Hidrolisis de enlaces ester (PET, PUR)",
                "ref": "Tokiwa et al. (2009) Appl. Microbiol. Biotechnol. 80:153-174"},
    "mnp": {"gen": "mnp (peroxidasa de manganeso)", "organismo": "Hongos ligninoliticos",
            "funcion": "Degradacion oxidativa de polimeros recalcitrantes",
            "ref": "Kirk & Farrell (1987) Annu. Rev. Microbiol. 41:465-505"},
    "16s": {"gen": "16S rRNA (V4)", "organismo": "Universal bacteriano",
            "funcion": "Gen de referencia / normalizador",
            "ref": "Apprill et al. (2015) AEM 81:e01441-15 (515F-Y/806R-A)"},
}

UNIVERSAL_PRIMERS = [
    {"nombre": "515F-Y", "secuencia": "GTGYCAGCMGCCGCGGTAA", "uso": "16S rRNA V4, forward",
     "ref": "Apprill et al. (2015) AEM 81:e01441-15"},
    {"nombre": "806R-A", "secuencia": "GGACTACNVGGGTWTCTAAT", "uso": "16S rRNA V4, reverse",
     "ref": "Parada et al. (2015) AEM 81:e01518-15"},
    {"nombre": "27F", "secuencia": "AGAGTTTGATCMTGGCTCAG", "uso": "16S rRNA casi completo, forward",
     "ref": "Lane (1991); Turner et al. (1999) J. Appl. Microbiol. 47:297-308"},
    {"nombre": "1492R", "secuencia": "TACGGYTACCTTGTTACGACTT", "uso": "16S rRNA, reverse",
     "ref": "Turner et al. (1999) J. Appl. Microbiol. 47:297-308"},
    {"nombre": "ITS1F", "secuencia": "TCCGTAGGTGAACCTGCGG", "uso": "ITS fungal, forward",
     "ref": "White et al. (1990) PCR Protocols, Academic Press"},
    {"nombre": "ITS4", "secuencia": "TCCTCCGCTTATTGATATGC", "uso": "ITS fungal, reverse",
     "ref": "White et al. (1990) PCR Protocols, Academic Press"},
]
PRIMER_DISCLAIMER = ("Secuencias universales clásicas publicadas; verificar contra la fuente "
                     "original antes de uso experimental. Para genes funcionales (petase, alkB, "
                     "lcc) los cebadores reportados son taxon-especificos y degenerados: diseñe "
                     "los suyos con el módulo de diseño (Tm/Wallace, GC, dímeros) a partir de "
                     "secuencias de referencia del NCBI.")

CONNECTORS = {
    "obis": "https://api.obis.org/v3/occurrence",
    "gbif": "https://api.gbif.org/v1/occurrence/search",
    "gbif_match": "https://api.gbif.org/v1/species/match",
    "worms": "https://www.marinespecies.org/rest/AphiaRecordsByMatchName",
    "inaturalist": "https://api.inaturalist.org/v1/observations",
    "crossref": "https://api.crossref.org/works",
    "pubmed": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
    "openaq": "https://api.openaq.org/v2/latest",
    "nasa_power": "https://power.larc.nasa.gov/api/temporal/daily/point",
    "stac_pc": "https://planetarycomputer.microsoft.com/api/stac/v1/search",
}

OLLAMA_URL = "http://localhost:11434"
OLLAMA_MODEL = "qwen3.5:4b"
OLLAMA_VISION_MODEL = "qwen3-vl:8b"

REFERENCE_VALUES_WATER = {
    "agua_potable_vigilancia": {"valor": 2.4, "unidad": "partículas/L",
                                 "nota": "Valor de vigilancia UE para agua de consumo humano "
                                         "(verificar estado vigente del Reglamento (UE) 2024/1448 "
                                         "y la Directiva (UE) 2020/2184 antes de citar)."},
}

DEFAULT_STATIONS = [
    ("JPC-01", "Canal Juan Polo (descarga urbana)", 10.4238, -75.5091, "impacto"),
    ("JPC-02", "Canal Juan Polo medio", 10.4285, -75.5060, "impacto"),
    ("VIR-01", "Sector La Esperanza", 10.4182, -75.5010, "impacto"),
    ("TUR-01", "Borde de manglar Turbana", 10.4098, -75.5012, "transicion"),
    ("TUR-02", "Turbana sur", 10.4052, -75.5077, "transicion"),
    ("MED-01", "Centro lagunar", 10.4308, -75.4953, "referencia_interna"),
    ("MED-02", "Norte centro", 10.4402, -75.4918, "referencia_interna"),
    ("NOR-01", "Sector norte", 10.4521, -75.4945, "referencia_interna"),
    ("NOR-02", "Punta norte", 10.4577, -75.4902, "transicion"),
    ("REF-01", "Boca marina (referencia externa)", 10.4650, -75.4760, "referencia"),
]
DEMO_REGION = "Ciénaga de La Virgen, Cartagena de Indias, Colombia"
DEMO_HULL = {"lat_min": 10.400, "lat_max": 10.470, "lon_min": -75.515, "lon_max": -75.470}

R_AVAILABLE_HINT = ("Instale R 4.x + paquetes vegan, DESeq2, lavaan, glmmTMB para ejecutar "
                    "los scripts generados en RStudio (https://posit.co).")
