"""Plantillas Excel de OMNIPLAST: estructura del proyecto, catalogos y ejemplo lleno.

Genera dos libros:
- PLANTILLA_OMNIPLAST.xlsx : estructura vacia con validaciones desplegables y guia LEEME.
- EJEMPLO_OMNIPLAST.xlsx   : la misma estructura llena con el dataset demo (Cienaga de
  La Virgen, DATOS SINTETICOS) para reproducir cualquier analisis de inmediato.
"""
import io
from datetime import datetime

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.datavalidation import DataValidation

from ..config import POLYMER_HAZARD, POLYMER_NAMES

NAVY = "0A2540"
TEAL = "00A896"
SAND = "F5F8FB"
GRAY = "808B96"

HEAD_FILL = PatternFill("solid", fgColor=NAVY)
HEAD_FONT = Font(bold=True, color="FFFFFF", size=10, name="Calibri")
SUB_FILL = PatternFill("solid", fgColor=TEAL)
NOTE_FONT = Font(italic=True, color=GRAY, size=9)
TITLE_FONT = Font(bold=True, size=14, color=NAVY, name="Calibri")
THIN = Border(*[Side(style="thin", color="D8E2EC")] * 4)

DICT_POLIMEROS = ",".join(list(POLYMER_NAMES.keys()))
DICT_MORF = "fibra,fragmento,film,espuma,gránulo,esfera"
DICT_COLOR = "azul,transparente,blanco,negro,rojo,verde,amarillo,naranja,rosa"
DICT_MATRIZ = "agua,sedimento,biota,aire"
DICT_ZONA = "impacto,transicion,referencia_interna,referencia"
DICT_TEMP = "secas,lluvias,transicion"
DICT_SI_NO = "SI,NO"
DICT_UNIDAD = "p/m3,p/kg dw,p/ind,p/L,particulas"
DICT_METODO_ID = "FTIR-ATR,FTIR microscopio (FPA),Raman,UV-visual + tinción (Nile Red),Piro-GC/MS,SEM-EDX,Identificación visual (hot needle test)"
DICT_GENES = "petase,mhetase,laccasa,alkb,phac,16s,intI1"
DICT_TIPO_BIO = "biomarcador,indice,morfometria,abundancia"
PARAM_CATALOG = [
    ("temperatura_c", "°C", "Temperatura del agua"),
    ("ph", "—", "pH"),
    ("od_mg_l", "mg/L", "Oxígeno disuelto"),
    ("dbo_mg_l", "mg/L", "Demanda bioquímica de O₂ (5 días)"),
    ("dqo_mg_l", "mg/L", "Demanda química de O₂"),
    ("turbidez_ntu", "NTU", "Turbidez"),
    ("salinidad_psu", "PSU", "Salinidad"),
    ("conductividad_us", "µS/cm", "Conductividad"),
    ("tss_mg_l", "mg/L", "Sólidos suspendidos totales"),
    ("nitrato_mg_l", "mg/L", "Nitrato (N-NO₃)"),
    ("nitrito_mg_l", "mg/L", "Nitrito (N-NO₂)"),
    ("amonio_mg_l", "mg/L", "Amonio (N-NH₄)"),
    ("fosfato_mg_l", "mg/L", "Fosfato (P-PO₄)"),
    ("clorofila_a_ug_l", "µg/L", "Clorofila a"),
    ("toc_mg_l", "mg/L", "Carbono orgánico total (agua)"),
    ("toc_pct", "%", "Materia orgánica (sedimento)"),
    ("arena_pct", "%", "Fracción arena"),
    ("limo_arcilla_pct", "%", "Fracción finos"),
    ("coliformes_totales_nmp", "NMP/100 mL", "Coliformes totales"),
    ("e_coli_nmp", "NMP/100 mL", "Escherichia coli"),
    ("sod_u_mg_prot", "U/mg prot", "Superóxido dismutasa"),
    ("cat_u_mg_prot", "U/mg prot", "Catalasa"),
    ("gst_u_mg_prot", "U/mg prot", "Glutatión-S-transferasa"),
    ("lpo_nmol_tba_mg", "nmol/mg", "Peroxidación lipídica (TBARS)"),
    ("factor_condicion", "—", "(peso/largo³)×100"),
    ("hsi_pct", "%", "Índice hepatosomático"),
    ("gsi_pct", "%", "Índice gonadosomático"),
]

SHEETS = {
    "Sitios": ["codigo", "nombre", "lat", "lon", "zona", "descripcion"],
    "Campanas": ["codigo", "nombre", "temporada", "fecha_inicio", "fecha_fin", "objetivo"],
    "Muestras": ["codigo", "sitio_codigo", "campana_codigo", "fecha", "matriz", "subtipo",
                 "profundidad_m", "replicas", "metodo", "colector", "blanco_proceso",
                 "duplicado", "notas"],
    "Microplasticos": ["muestra_codigo", "polimero", "morfologia", "color",
                       "tamano_min_um", "tamano_max_um", "abundancia", "unidad",
                       "peso_mg", "metodo_id", "notas"],
    "Parametros_Ambientales": ["muestra_codigo", "parametro", "valor", "unidad", "detectado"],
    "Biota": ["muestra_codigo", "especie", "tipo", "parametro", "valor", "unidad", "tejido"],
    "qPCR": ["muestra_codigo", "gen", "ct", "copias_g", "detectado"],
}

VALIDACIONES = {
    "Sitios": {"zona": DICT_ZONA},
    "Campanas": {"temporada": DICT_TEMP},
    "Muestras": {"matriz": DICT_MATRIZ, "blanco_proceso": DICT_SI_NO,
                 "duplicado": DICT_SI_NO},
    "Microplasticos": {"polimero": DICT_POLIMEROS, "morfologia": DICT_MORF,
                       "color": DICT_COLOR, "unidad": DICT_UNIDAD,
                       "metodo_id": DICT_METODO_ID},
    "Parametros_Ambientales": {"detectado": DICT_SI_NO},
    "Biota": {"tipo": DICT_TIPO_BIO},
    "qPCR": {"gen": DICT_GENES, "detectado": DICT_SI_NO},
}

ANCHOS = {
    "codigo": 18, "nombre": 42, "lat": 10, "lon": 10, "zona": 20, "descripcion": 34,
    "temporada": 14, "fecha_inicio": 12, "fecha_fin": 12, "objetivo": 34, "fecha": 12,
    "matriz": 12, "subtipo": 30, "profundidad_m": 13, "replicas": 10, "metodo": 30,
    "colector": 18, "notas": 24, "muestra_codigo": 26, "polimero": 12,
    "morfologia": 12, "color": 13, "tamano_min_um": 13, "tamano_max_um": 13,
    "abundancia": 12, "unidad": 12, "peso_mg": 11, "metodo_id": 26, "parametro": 22,
    "valor": 12, "unidad.1": 14, "detectado": 11, "especie": 30, "tipo": 14,
    "tejido": 16, "gen": 12, "ct": 8, "copias_g": 13,
}

LEEME = [
    ("OMNIPLAST — Plantilla de base de datos para análisis de microplásticos", "title"),
    (f"Generada: {datetime.now():%Y-%m-%d} · Estructura alineada con el modelo de datos "
     "de la plataforma (SQLite) y con protocolos GESAMP (2019) / Koelmans et al. (2019).",
     "sub"),
    ("", None),
    ("CÓMO USAR ESTA PLANTILLA (5 pasos)", "h2"),
    ("1. Complete la hoja 'Sitios' con sus estaciones georreferenciadas (use zona de "
     "impacto/transición/referencia según su diseño experimental).", None),
    ("2. Registre las 'Campanas' (temporadas: secas/lluvias) y luego cada 'Muestra' "
     "con su código único (ej.: JPC-01-C1-W00 = sitio-campaña-matriz-profundidad).", None),
    ("3. En 'Microplasticos' registre UNA FILA por combinación muestra×polímero×morfología"
     "×color×banda de tamaño. La columna abundancia es el N.º de partículas.", None),
    ("4. 'Parametros_Ambientales', 'Biota', 'qPCR' y 'Microbioma_16S' se enlazan por "
     "muestra_codigo. En qPCR use Ct=42 (o en blanco + detectado=NO) para no-detecciones.", None),
    ("5. Importe el archivo desde OMNIPLAST (página Proyectos → Cargar datos propios) o "
     "reemplace la base demo ejecutando el importador incluido.", None),
    ("", None),
    ("CONTROL DE CALIDAD OBLIGATORIO (Koelmans et al. 2019)", "h2"),
    ("· Blanco de proceso por lote (SI/NO en Muestras) y blancos de aire en laboratorio.", None),
    ("· Duplicados de campo ≥10% de las muestras.", None),
    ("· Registre SIEMPRE volumen filtrado, límite inferior de tamaño y método de "
     "identificación (FTIR/Raman/visual).", None),
    ("· Confirmación espectral ≥10% de partículas; reporte % de UNKNOWN.", None),
    ("", None),
    ("¿QUÉ ANÁLISIS QUIERE HACER? → DÓNDE EN OMNIPLAST", "h2"),
    ("Calidad del agua (WQI, TSI, nutrientes) → Dashboard · página 04 · agentes AG-007/008", None),
    ("Carga espacial de MP y hotspots → páginas 01 y 03 · agentes AG-019/026", None),
    ("Comparar zonas/temporadas estadísticamente → página 07 · agentes AG-035..039", None),
    ("Ordenaciones (PCA/NMDS/PERMANOVA) → página 08 · agentes AG-048..052", None),
    ("Genes de degradación (qPCR ΔΔCt) → página 11 · agentes AG-063/067", None),
    ("Microbioma y plastísfera → página 12 · agentes AG-068..082", None),
    ("Mapas, kriging, Moran, trayectorias → página 14 · agentes AG-093..096", None),
    ("Riesgo ecológico y One Health → página 16 · agentes AG-097/098", None),
    ("Índices integrados (PLI/PHI/IICAM) e informe DOCX → página 17", None),
    ("", None),
    ("UNIDADES Y FRANQUEZA", "h2"),
    ("Agua: partículas/m³ · Sedimento: partículas/kg peso seco · Biota: partículas/individuo.", None),
    ("Los umbrales e índices de OMNIPLAST son referenciales con cita primaria editable "
     "(config.py); no existe valor normativo ambiental internacional vigente para MP.", "note"),
    ("La hoja 'Ejemplo_Analisis' del libro EJEMPLO_OMNIPLAST.xlsx muestra un análisis "
     "guiado completo con estos mismos datos.", "note"),
]

GENERA_EJEMPLO = ["Pseudomonas", "Vibrio", "Bacillus", "Rhodococcus", "Stenotrophomonas",
                  "Acinetobacter", "Micrococcus", "Alcanivorax", "Marinobacter",
                  "Shewanella", "Flavobacterium", "Cytophaga", "Tenacibaculum",
                  "Oleispira", "Thalassospira", "Photobacterium", "Aeromonas",
                  "Enterococcus", "Escherichia/Shigella", "Bacteroides", "Nitrosomonas",
                  "Nitrospira", "Planctomyces", "Roseobacter-clade"]


def _estilar(ws, ncols, nrows=1):
    for j in range(1, ncols + 1):
        c = ws.cell(row=1, column=j)
        c.fill = HEAD_FILL
        c.font = HEAD_FONT
        c.border = THIN
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    for j in range(1, ncols + 1):
        ancho = ANCHOS.get(ws.cell(row=1, column=j).value, 16)
        ws.column_dimensions[get_column_letter(j)].width = max(ancho, 10)
    ws.freeze_panes = "A2"
    if nrows > 50:
        ws.auto_filter.ref = f"A1:{get_column_letter(ncols)}{nrows}"


def _validar(ws, ncols, nrows):
    headers = {ws.cell(row=1, column=j).value: j for j in range(1, ncols + 1)}
    for col_name, lista in VALIDACIONES.get(ws.title, {}).items():
        j = headers.get(col_name)
        if not j or len(lista) > 250:
            continue
        dv = DataValidation(type="list", formula1=f'"{lista}"', allow_blank=True,
                            showErrorMessage=True,
                            error="Use un valor del catálogo (ver hoja Catalogos).")
        ws.add_data_validation(dv)
        letra = get_column_letter(j)
        dv.add(f"{letra}2:{letra}{max(nrows, 200)}")


def _hoja_leeme(wb):
    ws = wb.create_sheet("LEEME")
    ws.sheet_properties.tabColor = TEAL
    ws.column_dimensions["A"].width = 118
    for i, (texto, estilo) in enumerate(LEEME, start=1):
        c = ws.cell(row=i, column=1, value=texto)
        if estilo == "title":
            c.font = TITLE_FONT
        elif estilo == "h2":
            c.font = Font(bold=True, size=11, color="FFFFFF")
            c.fill = SUB_FILL
        elif estilo == "sub":
            c.font = Font(italic=True, size=9, color=GRAY)
        elif estilo == "note":
            c.font = NOTE_FONT
        else:
            c.font = Font(size=10)
        c.alignment = Alignment(wrap_text=True, vertical="top")
    return ws


def _hoja_catalogos(wb):
    ws = wb.create_sheet("Catalogos")
    ws.sheet_properties.tabColor = NAVY
    fila = 1
    def tabla(titulo, headers, rows):
        nonlocal fila
        c = ws.cell(row=fila, column=1, value=titulo)
        c.font = Font(bold=True, size=11, color="FFFFFF")
        c.fill = SUB_FILL
        fila += 1
        for j, h in enumerate(headers, start=1):
            hc = ws.cell(row=fila, column=j, value=h)
            hc.fill = HEAD_FILL
            hc.font = HEAD_FONT
        fila += 1
        for r in rows:
            for j, v in enumerate(r, start=1):
                ws.cell(row=fila, column=j, value=v)
            fila += 1
        fila += 2
    tabla("POLÍMEROS (con puntaje de peligro — Lithner et al. 2011)",
          ["código", "nombre", "peligro"],
          [(k, POLYMER_NAMES[k], POLYMER_HAZARD.get(k, "sin puntaje publicado"))
           for k in POLYMER_NAMES])
    tabla("MORFOLOGÍAS", ["morfologia", "interpretación"],
          [("fibra", "textiles / aguas residuales"), ("fragmento", "degradación de objetos 3D"),
           ("film", "bolsas / empaques / mulch agrícola"), ("espuma", "EPS construcción y foodware"),
           ("gránulo", "materia prima industrial"), ("esfera", "microperlas cosméticas / pintura")])
    tabla("PARÁMETROS AMBIENTALES (código, unidad, descripción)",
          ["parametro", "unidad", "descripción"], PARAM_CATALOG)
    tabla("GENES DEL PANEL qPCR", ["gen", "uso"],
          [("petase", "Hidrólisis de PET (Yoshida 2016)"), ("mhetase", "MHET → TPA"),
           ("laccasa", "oxidación de aromáticos"), ("alkb", "ω-oxidación de alcanos"),
           ("phac", "síntesis PHA"), ("16s", "gen de referencia/normalizador"),
           ("intI1", "co-marcador de resistencia antimicrobiana (opcional)")])
    for col, ancho in zip("ABC", (24, 34, 52)):
        ws.column_dimensions[col].width = ancho
    return ws


def _hoja_ejemplo_analisis(wb, sites, samples, mp, env):
    """Hoja con un análisis guiado YA CALCULADO por el motor sobre este mismo dataset."""
    from ..core import indices
    ws = wb.create_sheet("Ejemplo_Analisis")
    ws.sheet_properties.tabColor = "E9B44C"
    fila = 1

    def titulo(texto, color=TEAL, size=11):
        nonlocal fila
        c = ws.cell(row=fila, column=1, value=texto)
        c.font = Font(bold=True, size=size, color="FFFFFF")
        c.fill = PatternFill("solid", fgColor=color)
        fila += 1

    def linea(texto, note=False):
        nonlocal fila
        c = ws.cell(row=fila, column=1, value=texto)
        c.font = NOTE_FONT if note else Font(size=10)
        c.alignment = Alignment(wrap_text=True, vertical="top")
        fila += 1

    def tabla_desde(df, titulo_txt):
        nonlocal fila
        titulo(titulo_txt)
        for j, col in enumerate(df.columns, start=1):
            hc = ws.cell(row=fila, column=j, value=str(col))
            hc.fill = HEAD_FILL
            hc.font = HEAD_FONT
            hc.border = THIN
        fila += 1
        for _, r in df.iterrows():
            for j, v in enumerate(r, start=1):
                if isinstance(v, float):
                    v = round(v, 2)
                ws.cell(row=fila, column=j, value=v).border = THIN
            fila += 1
        fila += 2

    titulo("EJEMPLO DE ANÁLISIS COMPLETO (calculado por el motor OMNIPLAST sobre ESTE dataset)",
           NAVY, 12)
    linea("Propósito: mostrar un flujo doctoral de extremo a extremo que usted puede "
          "reproducir en la app: Dashboard → página 17 (Índices) → agente AG-100. "
          "Todos los números de abajo salen de las demás hojas de este libro.", note=True)
    fila += 1

    agua = mp[mp["matriz"] == "agua"]
    sed = mp[mp["matriz"] == "sedimento"]
    por_agua = agua.groupby("sitio")["abundancia"].median()
    por_sed = sed.groupby("sitio")["abundancia"].median()
    fondo = por_sed.min()
    env_agua = env[env["matriz"] == "agua"]
    filas = []
    for sitio in sites["codigo"]:
        cf = indices.contamination_factor(por_sed.get(sitio, np.nan), fondo)
        wqi = indices.wqi_nsf(env_agua[env_agua["sitio"] == sitio]
                              .select_dtypes("number").mean(numeric_only=True).to_dict())
        phi = indices.polymer_hazard_index(sed[sed["sitio"] == sitio])
        filas.append({
            "sitio": sitio, "MP_agua_mediana": por_agua.get(sitio, np.nan),
            "clase_agua": indices.mp_tier(por_agua.get(sitio, 0), "agua"),
            "MP_sedimento_mediana": por_sed.get(sitio, np.nan),
            "CF_sedimento": cf,
            "clase_CF": indices.cf_classify(cf) if cf == cf else "—",
            "PHI": phi["phi"], "clase_PHI": phi["clase"],
            "WQI": wqi["wqi"], "clase_WQI": wqi["clase"],
        })
    resumen = pd.DataFrame(filas)
    tabla_desde(resumen, "1 · RADIOGRAFÍA POR ESTACIÓN (MP, contaminación, peligro y calidad del agua)")
    linea("Fórmulas: CF = C/C_fondo (Hakanson 1980) · PHI = Σ(puntaje×%abundancia)/100 "
          "(Lithner 2011) · WQI-NSF (Brown 1970). Clases editables en omniplast/config.py.",
          note=True)

    ranking = (resumen.dropna(subset=["MP_agua_mediana"])
               .sort_values("MP_agua_mediana", ascending=False)[["sitio", "MP_agua_mediana"]]
               .head(5).reset_index(drop=True))
    ranking.insert(1, "lectura", ["máxima carga", "segunda carga", "tercera carga",
                                  "cuarta carga", "quinta carga"][:len(ranking)])
    tabla_desde(ranking, "2 · RANKING DE CARGA EN AGUA (top 5) — prioriza hotspots Gi* (agente AG-026)")

    linea("3 · LECTURA DOCTORAL SUGERIDA (motor de recomendaciones, agente AG-100):", note=False)
    fila += 0
    if not resumen.empty:
        pli_max = indices.pli_site([resumen["CF_sedimento"].max()])
        wqi_medio = resumen["WQI"].mean()
        fibra_pct = (agua["morfologia"] == "fibra").mean() * 100
        pasos = [
            f"PLI del sistema = {pli_max:.2f} ({'deterioro' if pli_max > 1 else 'sin deterioro'}) → "
            "si >1: programa de reducción en fuentes terrestres focalizado en hotspots.",
            f"WQI medio = {wqi_medio:.0f}/100 → "
            + ("<50 exige saneamiento del sistema hidráulico." if wqi_medio < 50
               else "50-70: mantener vigilancia fisicoquímica trimestral."),
            f"Fibras = {fibra_pct:.0f}% del registro → "
            + (">40%: firma de vertimientos domésticos; filtros de microfibra y rejas en PTAR."
               if fibra_pct > 40 else "≤40%: mezcla de fuentes; completar con apportionment (página 03)."),
            "Siguiente paso en la app: página 17 → 'IICAM compuesto' para integrar todo con "
            "incertidumbre Monte Carlo, y página 20 → informe DOCX con descargo metodológico.",
        ]
        for p in pasos:
            linea("· " + p)
    for col, ancho in zip("ABCDEFGHIJ", (16, 16, 18, 16, 14, 18, 10, 18, 10, 18)):
        ws.column_dimensions[col].width = ancho
    ws.column_dimensions["A"].width = 22
    return ws


def _libro_base(con_datos):
    wb = Workbook()
    wb.remove(wb.active)
    _hoja_leeme(wb)
    return wb


def generar_plantilla_bytes():
    """Libro vacío con estructura, validaciones y catálogos (+3 filas ejemplo grises)."""
    wb = _libro_base(con_datos=False)
    ejemplos = {
        "Sitios": [["JPC-01", "Canal Juan Polo (descarga urbana)", 10.4238, -75.5091,
                    "impacto", "EJEMPLO — borre esta fila"],
                   ["REF-01", "Boca marina (referencia externa)", 10.4650, -75.4760,
                    "referencia", "EJEMPLO — borre esta fila"]],
        "Campanas": [["CAMP-24S", "Campaña Secas 2024", "secas", "2024-02-12",
                      "2024-02-16", "EJEMPLO — borre esta fila"]],
        "Muestras": [["JPC-01-C1-W00", "JPC-01", "CAMP-24S", "2024-02-12", "agua",
                      "superficie", 0.0, 3, "Red de neustón 330 um", "Su nombre",
                      "SI", "SI", "EJEMPLO"]],
        "Microplasticos": [["JPC-01-C1-W00", "PE", "fragmento", "azul", 100, 300, 14.0,
                            "p/m3", 0.01, "FTIR-ATR", "EJEMPLO"]],
        "Parametros_Ambientales": [["JPC-01-C1-W00", "od_mg_l", 2.8, "mg/L", "SI"]],
        "Biota": [["JPC-01-C1-BIO", "Crassostrea rhizophorae", "biomarcador",
                   "lpo_nmol_tba_mg", 9.4, "nmol/mg", "hepatopáncreas"]],
        "qPCR": [["JPC-01-C1-SED", "petase", 29.1, 400.5, "SI"]],
    }
    for nombre, columnas in SHEETS.items():
        ws = wb.create_sheet(nombre)
        ws.sheet_properties.tabColor = GRAY if nombre != "Sitios" else TEAL
        ws.append(columnas)
        for fila in ejemplos.get(nombre, []):
            ws.append(fila)
            for j in range(1, len(columnas) + 1):
                ws.cell(row=ws.max_row, column=j).font = NOTE_FONT
        _estilar(ws, len(columnas))
        _validar(ws, len(columnas), 3)
    ws = wb.create_sheet("Microbioma_16S")
    ws.append(["muestra_codigo"] + GENERA_EJEMPLO)
    ws.append(["JPC-01-C1-SED", 120, 340, 210, 60, 85, 90, 12, 40, 55, 45, 120, 80,
               30, 8, 25, 95, 70, 60, 110, 45, 35, 25, 70])
    _estilar(ws, len(GENERA_EJEMPLO) + 1)
    _hoja_catalogos(wb)
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


def generar_ejemplo_bytes():
    """Libro lleno con el dataset demo (Ciénaga de La Virgen, DATOS SINTÉTICOS)."""
    from ..db import dao
    wb = _libro_base(con_datos=True)
    sites = dao.sites_df()
    camps = dao.campaigns_df()
    samples = dao.samples_df()
    mp = dao.mp_df()
    env = dao.env_df()
    bio = dao.bio_df()
    qpcr = dao.qpcr_df()
    mb = dao.microbiome_wide()

    def _dump(nombre, df):
        ws = wb.create_sheet(nombre)
        ws.append(list(df.columns))
        for row in df.itertuples(index=False):
            ws.append(list(row))
        _estilar(ws, len(df.columns), len(df) + 1)
        _validar(ws, len(df.columns), len(df) + 1)

    _dump("Sitios", sites)
    camps2 = camps.copy()
    camps2["fecha_inicio"] = camps2["fecha_inicio"].dt.strftime("%Y-%m-%d")
    camps2["fecha_fin"] = camps2["fecha_fin"].dt.strftime("%Y-%m-%d")
    _dump("Campanas", camps2[["codigo", "nombre", "temporada", "fecha_inicio",
                              "fecha_fin", "objetivo"]])
    smp = samples.copy()
    smp["fecha"] = pd.to_datetime(smp["fecha"]).dt.strftime("%Y-%m-%d")
    qa = smp["qa_qc"].apply(lambda d: d if isinstance(d, dict) else {})
    smp = smp.assign(blanco_proceso=qa.apply(lambda d: "SI" if d.get("blanco_proceso") else "NO"),
                     duplicado=qa.apply(lambda d: "SI" if d.get("duplicado") else "NO"))
    _dump("Muestras", smp[["codigo", "sitio", "campana", "fecha", "matriz", "subtipo",
                           "profundidad_m", "replicas", "metodo", "colector",
                           "blanco_proceso", "duplicado"]].rename(
                               columns={"sitio": "sitio_codigo", "campana": "campana_codigo"}))
    _dump("Microplasticos", mp[["codigo", "polimero", "morfologia", "color",
                                "tamano_min_um", "tamano_max_um", "abundancia", "unidad",
                                "peso_mg", "metodo_id"]].rename(columns={"codigo": "muestra_codigo"}))
    _dump("Parametros_Ambientales", env[["codigo", "parametro", "valor", "unidad"]].rename(
        columns={"codigo": "muestra_codigo"}))
    if not bio.empty:
        _dump("Biota", bio[["codigo", "taxon", "tipo", "parametro", "valor", "unidad",
                            "tejido"]].rename(columns={"codigo": "muestra_codigo",
                                                       "taxon": "especie"}))
    if not qpcr.empty:
        _dump("qPCR", qpcr[["codigo", "target", "ct", "valor", "detectado"]].rename(
            columns={"codigo": "muestra_codigo", "target": "gen", "valor": "copias_g"}))
    if not mb.empty:
        taxa = [c for c in mb.columns if c not in
                ("sample_id", "sitio", "zona", "campana", "temporada", "matriz", "anio")]
        _dump("Microbioma_16S", mb[["codigo"] + taxa].rename(columns={"codigo": "muestra_codigo"}))
    _hoja_ejemplo_analisis(wb, sites, samples, mp, dao.env_wide())
    _hoja_catalogos(wb)
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


def guardar_en_exports():
    from ..config import EXPORT_DIR
    EXPORT_DIR.mkdir(parents=True, exist_ok=True)
    p1 = EXPORT_DIR / "PLANTILLA_OMNIPLAST.xlsx"
    p1.write_bytes(generar_plantilla_bytes())
    p2 = EXPORT_DIR / "EJEMPLO_OMNIPLAST.xlsx"
    p2.write_bytes(generar_ejemplo_bytes())
    return p1, p2
