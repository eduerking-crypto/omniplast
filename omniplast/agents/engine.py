"""Motor de ejecucion de agentes: contexto de datos, despacho y bitacora."""
import time
import traceback

from ..db import dao
from .handlers_a import HANDLERS_A
from .handlers_b import HANDLERS_B
from .handlers_c1 import HANDLERS_C1
from .handlers_c2 import HANDLERS_C2
from .handlers_c3 import HANDLERS_C3
from .registry import CATEGORIES, REGISTRY

HANDLERS = {**HANDLERS_A, **HANDLERS_B, **HANDLERS_C1, **HANDLERS_C2, **HANDLERS_C3}


def build_context(params=None):
    return {
        "mp": dao.mp_df(),
        "env": dao.env_wide(),
        "bio": dao.bio_df(),
        "qpcr": dao.qpcr_df(),
        "mb": dao.microbiome_wide(),
        "samples": dao.samples_df(),
        "sites": dao.sites_df(),
        "params": params or {},
    }


def run_agent(agent_id, params=None):
    spec = REGISTRY.get(agent_id, {"id": agent_id, "nombre": agent_id})
    handler = HANDLERS.get(agent_id)
    if handler is None:
        return {"estado": "error", "resumen": f"Agente {agent_id} no registrado.",
                "tablas": {}, "figuras": {}, "refs": [], "pasos": [], "agente": spec}
    t0 = time.time()
    try:
        ctx = build_context(params)
        result = handler(ctx)
    except Exception as e:
        result = {"estado": "error", "resumen": f"Error interno: {e}", "tablas": {},
                  "figuras": {}, "refs": spec.get("refs", []),
                  "pasos": [traceback.format_exc(limit=3)]}
    dur = round(time.time() - t0, 2)
    result["duracion_s"] = dur
    result["agente"] = spec
    try:
        dao.log_agent_run(agent_id, spec.get("nombre", agent_id), result.get("estado", "?"),
                          str(result.get("resumen", ""))[:400], {}, dur)
    except Exception:
        pass
    return result


def run_category(categoria, params=None):
    return {aid: run_agent(aid, params)
            for aid, meta in REGISTRY.items() if meta["categoria"] == categoria}
