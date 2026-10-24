import sys
sys.path.insert(0, ".")
import collections
from omniplast.agents import REGISTRY, run_agent

estados = collections.Counter()
fallos = []
for aid in sorted(REGISTRY):
    try:
        r = run_agent(aid)
        estados[r["estado"]] += 1
        if r["estado"] == "error":
            fallos.append((aid, r["resumen"], r.get("pasos", [])))
    except Exception as e:
        estados["excepcion"] += 1
        fallos.append((aid, str(e), []))
print("ESTADOS:", dict(estados), "total:", len(REGISTRY))
for aid, msg, pasos in fallos:
    print(f"\n[{aid}] {msg}")
    if pasos:
        print(pasos[0][:600])
