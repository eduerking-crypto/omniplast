"""API REST de OMNIPLAST (FastAPI + SQLite).

Ejecución: python -m uvicorn api.main:app --port 8010
Documentación interactiva: http://127.0.0.1:8010/docs
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pandas as pd  # noqa: E402
from fastapi import FastAPI, HTTPException, Query  # noqa: E402
from pydantic import BaseModel  # noqa: E402

from omniplast.agents import REGISTRY, run_agent  # noqa: E402
from omniplast.assistant import ask, recommend_test  # noqa: E402
from omniplast.core import indices  # noqa: E402
from omniplast.core.connectors import connectors  # noqa: E402
from omniplast.db import dao  # noqa: E402
from omniplast.db.seed import seed_demo_if_empty  # noqa: E402
from omniplast.reporting import r_scripts  # noqa: E402
from omniplast.utils import to_jsonable  # noqa: E402

app = FastAPI(
    title="OMNIPLAST API",
    description="API de la plataforma doctoral de análisis integral de microplásticos. "
                "Base de datos: SQLite con proyecto demo (Ciénaga de La Virgen, datos "
                "sintéticos etiquetados).",
    version="1.0.0",
)


@app.on_event("startup")
def startup():
    seed_demo_if_empty()


@app.get("/api/health")
def health():
    return {"ok": True, "app": "OMNIPLAST", "version": "1.0.0"}


@app.get("/api/resumen")
def resumen():
    return to_jsonable(dao.project_summary())


@app.get("/api/sitios")
def sitios():
    return to_jsonable(dao.sites_df().to_dict("records"))


@app.get("/api/muestras")
def muestras(matriz: str = Query(default=None, description="agua|sedimento|biota")):
    df = dao.samples_df()
    if matriz:
        df = df[df["matriz"] == matriz]
    return to_jsonable(df.head(500).to_dict("records"))


@app.get("/api/microplasticos")
def microplasticos(matriz: str = None, sitio: str = None, limit: int = 2000):
    df = dao.mp_df()
    if matriz:
        df = df[df["matriz"] == matriz]
    if sitio:
        df = df[df["sitio"] == sitio]
    return to_jsonable(df.head(limit).to_dict("records"))


@app.get("/api/ambiente")
def ambiente(matriz: str = None, limit: int = 2000):
    df = dao.env_wide()
    if matriz:
        df = df[df["matriz"] == matriz]
    return to_jsonable(df.head(limit).to_dict("records"))


@app.get("/api/agentes")
def listar_agentes(categoria: str = None):
    if categoria:
        return {k: v for k, v in REGISTRY.items() if v["categoria"] == categoria}
    return REGISTRY


class AgentRequest(BaseModel):
    agente_id: str
    params: dict = {}


@app.post("/api/agentes/ejecutar")
def ejecutar_agente(req: AgentRequest):
    if req.agente_id not in REGISTRY:
        raise HTTPException(404, f"Agente {req.agente_id} no existe")
    r = run_agent(req.agente_id, req.params)
    r["tablas"] = {k: (v.to_dict("records") if isinstance(v, pd.DataFrame) else v)
                   for k, v in r.get("tablas", {}).items()}
    r["figuras"] = {k: fig.to_json() for k, fig in r.get("figuras", {}).items()}
    return to_jsonable(r)


class AsistenteRequest(BaseModel):
    pregunta: str
    usar_llm: bool = True


@app.post("/api/asistente")
def asistente(req: AsistenteRequest):
    return ask(req.pregunta, use_llm=req.usar_llm)


class TestRequest(BaseModel):
    n_grupos: int = 2
    pareado: bool = False
    normalidad: bool | None = None
    heterocedasticidad: bool | None = None
    conteos: bool = False
    multivariado: bool = False
    repetido: bool = False


@app.post("/api/asesor/estadistico")
def asesor(req: TestRequest):
    return recommend_test(n_grupos=req.n_grupos, pareado=req.pareado,
                           normalidad=req.normalidad,
                           heterocedasticidad=req.heterocedasticidad,
                           conteos=req.conteos, multivariado=req.multivariado,
                           repetido=req.repetido)


@app.get("/api/indices/qli")
def indices_sitio():
    mp = dao.mp_df()
    if mp.empty:
        raise HTTPException(404, "Sin datos de MP")
    out = {}
    for matriz in ("agua", "sedimento"):
        per = mp[mp["matriz"] == matriz].groupby("sitio")["abundancia"].sum()
        if per.empty:
            continue
        bg = per.min()
        out[matriz] = {s: {"total": float(t),
                           "cf": indices.contamination_factor(t, bg),
                           "pli": indices.pli_site([indices.contamination_factor(t, bg)]),
                           "clase_mp": indices.mp_tier(t, matriz)}
                       for s, t in per.items()}
    return to_jsonable(out)


@app.get("/api/conectores/{nombre}")
def conector(nombre: str, q: str = None, lat: float = None, lon: float = None,
             limit: int = 100):
    if nombre == "obis":
        r = connectors.search_obis(scientific_name=q, limit=limit)
    elif nombre == "gbif":
        r = connectors.search_gbif(scientific_name=q, limit=limit)
    elif nombre == "inaturalist":
        r = connectors.search_inaturalist(taxon_name=q, lat=lat, lon=lon, limit=limit)
    elif nombre == "crossref":
        r = connectors.search_crossref(q or "microplastics", rows=limit)
    elif nombre == "pubmed":
        r = connectors.search_pubmed(q or "microplastics", rows=limit)
    elif nombre == "openaq":
        r = connectors.fetch_openaq(lat=lat, lon=lon)
    elif nombre == "nasa_power":
        r = connectors.fetch_nasa_power(lat or 10.43, lon or -75.49)
    else:
        raise HTTPException(404, f"Conector {nombre} no existe")
    if not r.get("ok"):
        raise HTTPException(502, r.get("error", "fallo del conector"))
    if "df" in r:
        r["df"] = r["df"].to_dict("records")
    return to_jsonable(r)


@app.get("/api/reportes/scripts_r")
def scripts_r():
    r = r_scripts.generate_all()
    return {"ok": True, "archivos": {k: str(v) for k, v in r["scripts"].items()},
            "nota": r["nota"]}


class WQIRequest(BaseModel):
    od_mg_l: float | None = None
    coliformes_totales_nmp: float | None = None
    ph: float | None = None
    dbo_mg_l: float | None = None
    temperatura_c: float | None = None
    fosfato_mg_l: float | None = None
    nitrito_mg_l: float | None = None
    turbidez_ntu: float | None = None
    tss_mg_l: float | None = None


@app.post("/api/indices/wqi")
def calcular_wqi(req: WQIRequest):
    return indices.wqi_nsf(req.model_dump())
