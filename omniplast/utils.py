"""Utilidades genericas: formato, correccion multiple, geo, serializacion."""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


def safe_div(a, b, default=np.nan):
    try:
        if b is None or (isinstance(b, float) and math.isnan(b)) or b == 0:
            return default
        return a / b
    except Exception:
        return default


def holm(pvals):
    """Correccion de Holm (1979) para comparaciones multiples."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adj = np.empty(n)
    running = 0.0
    for rank, idx in enumerate(order):
        val = (n - rank) * p[idx]
        running = max(running, val)
        adj[idx] = min(1.0, running)
    return adj


def benjamini_hochberg(pvals):
    """Correccion de falsa tasa de descubrimiento (Benjamini-Hochberg 1995)."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.minimum(1.0, ranked)
    return adj


def interpret_p(p, alpha=0.05):
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return "no computable"
    return f"p = {p:.4f} → {'significativo' if p < alpha else 'no significativo'} (α = {alpha})"


def fmt_num(x, dec=2):
    try:
        if x is None or (isinstance(x, float) and math.isnan(x)):
            return "—"
        if abs(x) >= 1e5 or (0 < abs(x) < 1e-3):
            return f"{x:.2e}"
        return f"{x:,.{dec}f}"
    except Exception:
        return str(x)


def classify_value(x, classes, labels=None):
    """Clasifica un valor contra una lista [(etiqueta, min, max), ...]."""
    for i, (label, lo, hi) in enumerate(classes):
        if lo <= x <= hi:
            return label
    return classes[-1][0] if classes else None


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0088
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = p2 - p1
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def robust_z(x):
    x = np.asarray(x, dtype=float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med)) * 1.4826
    if mad == 0 or np.isnan(mad):
        sd = np.nanstd(x)
        return (x - med) / (sd if sd > 0 else 1.0)
    return (x - med) / mad


def to_jsonable(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj) if np.isfinite(obj) else None
    if isinstance(obj, (np.ndarray,)):
        return obj.tolist()
    if isinstance(obj, (pd.Timestamp,)):
        return obj.isoformat()
    if isinstance(obj, pd.DataFrame):
        return {"dataframe": True, "rows": len(obj)}
    if isinstance(obj, dict):
        return {k: to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    return obj


def save_payload(payload, name, subdir="exports"):
    """Guarda un payload JSON en data/<subdir> con nombre seguro."""
    out_dir = Path(__file__).resolve().parents[1] / "data" / subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in str(name))
    path = out_dir / f"{safe}.json"
    path.write_text(json.dumps(to_jsonable(payload), ensure_ascii=False, indent=2),
                    encoding="utf-8")
    return path


def percentile_score(x, lo, hi):
    """Score 0-100 por posicion de x en un rango de referencia (clamp)."""
    if hi == lo:
        return 50.0
    return float(np.clip((x - lo) / (hi - lo) * 100, 0, 100))


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)
    return Path(path)


def col_ci95(series):
    x = pd.to_numeric(series, errors="coerce").dropna().values
    if len(x) < 2:
        return (np.nan, np.nan)
    m = np.mean(x)
    se = np.std(x, ddof=1) / np.sqrt(len(x))
    from scipy import stats
    t = stats.t.ppf(0.975, len(x) - 1)
    return (m - t * se, m + t * se)
