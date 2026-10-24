"""Cliente Ollama local para PLASTI con fallback degradado y elegante.

GOTCHAS de Ollama aplicados (memoria de sesiones previas):
- qwen3-vl enruta texto al campo 'thinking'; probar message.content y luego thinking.
- num_predict bajo trunca: elevarlo a 700+.
"""
import requests

from ..config import OLLAMA_MODEL, OLLAMA_URL
from .knowledge_base import context_for_prompt, search_kb
from . import advisor

SYSTEM_PROMPT = (
    "Eres PLASTI, asistente experto doctoral de la plataforma OMNIPLAST (análisis integral "
    "de microplásticos en ecosistemas). Respondes en español, con rigor científico, "
    "citas de las fuentes provistas, sin inventar valores normativos y aclarando cuando "
    "algo requiere verificación. Estructura: respuesta directa, luego 'Cómo hacerlo en "
    "OMNIPLAST' (agente/página recomendada), y 'Referencias'. Sé conciso (máx. 250 palabras)."
)


def ollama_available():
    try:
        r = requests.get(f"{OLLAMA_URL}/api/tags", timeout=3)
        return r.status_code == 200
    except Exception:
        return False


def chat_ollama(prompt, model=None, timeout=180):
    try:
        r = requests.post(
            f"{OLLAMA_URL}/api/chat",
            json={"model": model or OLLAMA_MODEL,
                  "messages": [{"role": "system", "content": SYSTEM_PROMPT},
                               {"role": "user", "content": prompt}],
                  "stream": False,
                  "options": {"num_predict": 700, "temperature": 0.4}},
            timeout=timeout)
        r.raise_for_status()
        data = r.json()
        msg = data.get("message", {})
        content = (msg.get("content") or "").strip()
        if not content:
            content = (msg.get("thinking") or "").strip()
        return {"ok": bool(content), "respuesta": content,
                "modelo": data.get("model", model or OLLAMA_MODEL), "fuente": "ollama"}
    except Exception as e:
        return {"ok": False, "respuesta": "", "error": str(e), "fuente": "ollama"}


def offline_answer(prompt):
    """Respuesta sin LLM: recuperacion de KB + asesor estadistico + indices clave."""
    kb_hits = search_kb(prompt, top=3)
    partes = []
    if kb_hits:
        partes.append("**PLASTI (modo sin LLM — base de conocimiento local):**")
        for h in kb_hits:
            partes.append(f"\n**{h['tema']}**\n{h['resumen']}\n{h['detalle']}\n"
                          f"Fuentes: {'; '.join(h['fuentes'])}")
    q = prompt.lower()
    if any(w in q for w in ("prueba", "test", "anova", "kruskal", "comparar", "estadíst")):
        rec = advisor.recommend_test(multivariado="comunidad" in q or "multivari" in q,
                                     n_grupos=3 if "zonas" in q else 2,
                                     conteos="conteo" in q or "abundancia" in q)
        partes.append(f"\n**Recomendación del asesor estadístico:** {rec['prueba']}\n"
                      f"Justificación: {rec['justificacion']}\nReferencia: {rec['ref']}")
    if any(w in q for w in ("cuántas", "cuantas", "muestras", "potencia", "n mínimo")):
        dis = advisor.sampling_design()
        partes.append(f"\n**Diseño muestral sugerido:** {dis['n_por_grupo']} por grupo, "
                      f"{dis['campanas_recomendadas']} campañas, {dis['total_muestras_aprox']} "
                      f"muestras totales aproximadas.")
    if not partes:
        partes.append(
            "PLASTI (modo sin LLM): no encontré una entrada exacta. Pruebe con palabras "
            "clave como 'protocolo agua', 'QA/QC', 'PERMANOVA', 'PHI', 'PETasa', "
            "'diseño muestral', o active Ollama para respuestas generativas "
            "(ollama pull qwen3.5:4b).")
    return {"ok": True, "respuesta": "\n".join(partes), "fuente": "kb_local",
            "modelo": "reglas+KB"}


def ask(prompt, use_llm=True):
    if use_llm and ollama_available():
        contexto = context_for_prompt(prompt)
        full = f"{prompt}\n\n{contexto}" if contexto else prompt
        res = chat_ollama(full)
        if res["ok"]:
            return res
    return offline_answer(prompt)
