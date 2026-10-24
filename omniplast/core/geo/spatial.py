"""Geoestadistica: IDW, kriging ordinario, Moran I global/local, hotspots Gi*."""
import numpy as np
import pandas as pd

from ...utils import haversine_km


def build_grid(bbox, n=40):
    lats = np.linspace(bbox["lat_min"], bbox["lat_max"], n)
    lons = np.linspace(bbox["lon_min"], bbox["lon_max"], n)
    return np.meshgrid(lons, lats)


def idw_interpolate(lat, lon, values, grid_lon, grid_lat, power=2):
    """Inverse Distance Weighting (Shepard 1968)."""
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    v = np.asarray(values, float)
    ok = np.isfinite(v)
    lat, lon, v = lat[ok], lon[ok], v[ok]
    if len(v) == 0:
        return np.full(grid_lon.shape, np.nan)
    dist = np.sqrt(((grid_lon[:, :, None] - lon) * 111.32 * np.cos(np.radians(lat))) ** 2 +
                   ((grid_lat[:, :, None] - lat) * 110.57) ** 2)
    dist = np.maximum(dist, 1e-6)
    w = 1 / dist ** power
    return (w @ v) / w.sum(axis=2)


def variogram_empirical(lat, lon, values, n_bins=10, max_dist_km=None):
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    v = np.asarray(values, float)
    n = len(v)
    if n < 4:
        return pd.DataFrame(), None
    d = []
    g = []
    for i in range(n):
        for j in range(i + 1, n):
            d.append(haversine_km(lat[i], lon[i], lat[j], lon[j]))
            g.append((v[i] - v[j]) ** 2 / 2)
    d = np.asarray(d)
    g = np.asarray(g)
    max_d = max_dist_km or np.percentile(d, 95)
    edges = np.linspace(0, max_d, n_bins + 1)
    rows = []
    for b in range(n_bins):
        m = (d > edges[b]) & (d <= edges[b + 1])
        if m.sum() >= 1:
            rows.append({"dist_km": float(d[m].mean()), "semivarianza": float(g[m].mean()),
                         "n_pares": int(m.sum())})
    return pd.DataFrame(rows), float(max_d)


def _spherical(h, c0, c1, a):
    return np.where(h == 0, 0, c0 + c1 * (1.5 * h / a - 0.5 * (h / a) ** 3) * (h <= a)
                    + c1 * (h > a))


def _exponential(h, c0, c1, a):
    return c0 + c1 * (1 - np.exp(-3 * h / a))


def _gaussian(h, c0, c1, a):
    return c0 + c1 * (1 - np.exp(-3 * (h / a) ** 2))


def fit_variogram(emp_variogram, model="spherical"):
    from scipy.optimize import curve_fit
    if emp_variogram is None or emp_variogram.empty:
        return None
    h = emp_variogram["dist_km"].values
    g = emp_variogram["semivarianza"].values
    fn = {"spherical": _spherical, "exponential": _exponential,
          "gaussian": _gaussian}[model]
    c0_0 = float(np.min(g) / 2)
    c1_0 = float(np.max(g) - c0_0)
    a_0 = float(np.max(h) / 2) or 1.0
    try:
        popt, _ = curve_fit(fn, h, g, p0=[c0_0, c1_0, a_0], maxfev=8000)
        ss_res = float(((g - fn(h, *popt)) ** 2).sum())
        ss_tot = float(((g - g.mean()) ** 2).sum()) or 1e-12
        return {"modelo": model, "nugget": float(popt[0]), "meseta": float(popt[0] + popt[1]),
                "alcance_km": float(popt[2]), "r2": 1 - ss_res / ss_tot, "fn": lambda x: fn(x, *popt)}
    except Exception:
        return None


def ordinary_kriging(lat, lon, values, grid_lon, grid_lat, model="spherical"):
    """Kriging ordinario con variograma ajustado (solucion matricial por punto)."""
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    v = np.asarray(values, float)
    ok = np.isfinite(v)
    lat, lon, v = lat[ok], lon[ok], v[ok]
    n = len(v)
    if n < 5:
        return np.full(grid_lon.shape, np.nan), None
    emp, _ = variogram_empirical(lat, lon, v)
    fit = fit_variogram(emp, model)
    fn = (fit["fn"] if fit else (lambda x: np.maximum(x * 0, np.nanvar(v))))
    coords = np.column_stack([lat * 110.57, lon * 111.32 * np.cos(np.radians(lat.mean()))])
    A = np.zeros((n + 1, n + 1))
    for i in range(n):
        for j in range(n):
            A[i, j] = fn(np.hypot(*(coords[i] - coords[j]))) if i != j else 0
    A[n, :n] = 1
    A[:n, n] = 1
    A[n, n] = 0
    try:
        inv_A = np.linalg.pinv(A)
    except np.linalg.LinAlgError:
        return np.full(grid_lon.shape, np.nan), fit
    est = np.full(grid_lon.shape, np.nan)
    flat_lon = grid_lon.ravel()
    flat_lat = grid_lat.ravel()
    for k in range(len(flat_lon)):
        p = np.array([flat_lat[k] * 110.57, flat_lon[k] * 111.32 * np.cos(np.radians(lat.mean()))])
        d = np.array([fn(np.hypot(*(coords[i] - p))) for i in range(n)])
        b = np.append(d, 1)
        w = inv_A @ b
        est.flat[k] = w[:n] @ v
    return est, fit


def morans_i(lat, lon, values, k=8, permutations=999, seed=42):
    """I de Moran global con pesos k-NN y test por permutacion."""
    rng = np.random.default_rng(seed)
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    v = np.asarray(values, float)
    ok = np.isfinite(v)
    lat, lon, v = lat[ok], lon[ok], v[ok]
    n = len(v)
    if n < 5:
        return {"test": "Moran I", "error": "n≥5 puntos"}
    W = np.zeros((n, n))
    for i in range(n):
        d = np.array([haversine_km(lat[i], lon[i], lat[j], lon[j]) for j in range(n)])
        d[i] = np.inf
        nn = np.argsort(d)[:k]
        W[i, nn] = 1
    W = W / W.sum(axis=1, keepdims=True)
    z = v - v.mean()
    num = (W * np.outer(z, z)).sum()
    den = (z ** 2).sum()
    I = (n / W.sum()) * num / den if den > 0 else np.nan
    count = 0
    for _ in range(permutations):
        vp = rng.permutation(v)
        zp = vp - vp.mean()
        Ip = (n / W.sum()) * (W * np.outer(zp, zp)).sum() / (zp ** 2).sum() if (zp ** 2).sum() > 0 else 0
        if Ip >= I:
            count += 1
    p = (count + 1) / (permutations + 1)
    return {"test": "I de Moran global (k-NN)", "I": float(I), "p": float(p),
            "permutaciones": permutations, "k_vecinos": k,
            "interpretacion": f"I={I:.3f} ({'agrupamiento espacial significativo' if p < .05 and I > 0 else 'patrón aleatorio o sin autocorrelación significativa'}), p={p:.4f}"}


def local_moran(lat, lon, values, k=8):
    """LISA: clasificacion HH/LL/HL/LH por punto."""
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    v = np.asarray(values, float)
    ok = np.isfinite(v)
    lat, lon, v = lat[ok], lon[ok], v[ok]
    n = len(v)
    if n < 5:
        return pd.DataFrame()
    W = np.zeros((n, n))
    for i in range(n):
        d = np.array([haversine_km(lat[i], lon[i], lat[j], lon[j]) for j in range(n)])
        d[i] = np.inf
        nn = np.argsort(d)[:k]
        W[i, nn] = 1
    W = W / W.sum(axis=1, keepdims=True)
    z = v - v.mean()
    zl = W @ z
    rows = []
    for i in range(n):
        if z[i] > 0 and zl[i] > 0:
            cls = "HH (hotspot)"
        elif z[i] < 0 and zl[i] < 0:
            cls = "LL (coldspot)"
        elif z[i] > 0 and zl[i] < 0:
            cls = "HL (atípico alto)"
        else:
            cls = "LH (atípico bajo)"
        rows.append({"lat": lat[i], "lon": lon[i], "valor": v[i], "lisa_z": z[i],
                     "lisa_lag": zl[i], "clase": cls})
    return pd.DataFrame(rows)


def getis_ord_gi_star(lat, lon, values, distance_km=None, k=8):
    """Hotspot analysis Gi* con vecindario por distancia o k-NN."""
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    v = np.asarray(values, float)
    ok = np.isfinite(v)
    lat, lon, v = lat[ok], lon[ok], v[ok]
    n = len(v)
    if n < 5:
        return pd.DataFrame()
    W = np.zeros((n, n))
    for i in range(n):
        if distance_km:
            for j in range(n):
                W[i, j] = 1 if haversine_km(lat[i], lon[i], lat[j], lon[j]) <= distance_km else 0
        else:
            d = np.array([haversine_km(lat[i], lon[i], lat[j], lon[j]) for j in range(n)])
            d[i] = np.inf
            W[i, np.argsort(d)[:k]] = 1
    m = v.mean()
    s = v.std(ddof=1)
    rows = []
    for i in range(n):
        wi = W[i]
        Wi = wi.sum()
        if Wi == 0 or s == 0:
            z = 0.0
        else:
            z = (wi @ v - Wi * m) / (s * np.sqrt((n * (wi ** 2).sum() - Wi ** 2) / (n - 1)))
        conf = ("99%" if z > 2.58 else "95%" if z > 1.96 else
                "90%" if z > 1.645 else "—")
        rows.append({"lat": lat[i], "lon": lon[i], "valor": v[i], "gi_z": float(z),
                     "tipo": ("hotspot" if z > 1.645 else "coldspot" if z < -1.645 else "no significativo"),
                     "confianza": conf})
    return pd.DataFrame(rows)
