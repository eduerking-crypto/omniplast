"""Teledeteccion: indices espectrales acuaticos, deteccion de plasticos, STAC.

Indices con longitudes de onda Sentinel-2 MSI (B2 490, B3 560, B4 665, B5 705,
B6 740, B7 783, B8 842, B8A 865, B11 1610, B12 2190 nm).
FDI (Floating Debris Index): Biermann et al. (2020) Sci. Rep. 10:19836,
variante del FAI de Hu (2009).
"""
import numpy as np

S2_WAVELENGTHS = {"B2": 490, "B3": 560, "B4": 665, "B5": 705, "B6": 740,
                  "B7": 783, "B8": 842, "B8A": 865, "B11": 1610, "B12": 2190}


def _norm_diff(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        return (a - b) / (a + b)


def compute_index(bands, name):
    """bands: dict Banda -> array 2D. Devuelve array del índice."""
    if name == "NDWI":
        return _norm_diff(bands["B3"], bands["B8"])
    if name == "MNDWI":
        return _norm_diff(bands["B3"], bands["B11"])
    if name == "NDTI":
        return _norm_diff(bands["B4"], bands["B3"])
    if name == "NDCI":
        return _norm_diff(bands["B5"], bands["B4"])
    if name == "NDVI":
        return _norm_diff(bands["B8"], bands["B4"])
    if name == "NDMI":
        return _norm_diff(bands["B8A"], bands["B11"])
    if name == "FDI":
        red, nir, swir = bands["B4"], bands["B8"], bands["B11"]
        lam = (S2_WAVELENGTHS["B8"] - S2_WAVELENGTHS["B4"]) / \
              (S2_WAVELENGTHS["B11"] - S2_WAVELENGTHS["B4"])
        return nir - (red + (swir - red) * lam)
    if name == "PI":
        with np.errstate(divide="ignore"):
            return bands["B8"] / np.maximum(bands["B4"], 1e-9)
    if name == "2BDA_chla":
        with np.errstate(divide="ignore"):
            return bands["B5"] / np.maximum(bands["B4"], 1e-9)
    if name == "SABI":
        with np.errstate(divide="ignore"):
            return (bands["B8"] - bands["B4"]) / (bands["B8"] + bands["B4"] + 1e-9)
    raise ValueError(f"Índice {name} no implementado")


INDEX_INFO = {
    "NDWI": ("NDWI (McFeeters 1996)", "Cuerpos de agua; >0.2 muy probable agua"),
    "MNDWI": ("MNDWI (Xu 2006)", "Agua con SWIR; discrimina vegetación mejor que NDWI"),
    "NDTI": ("NDTI (Lacaux et al. 2007)", "Turbidez; mayor valor = agua más turbia"),
    "NDCI": ("NDCI (Mishra & Mishra 2012)", "Clorofila-a en aguas productivas"),
    "FDI": ("FDI (Biermann et al. 2020)", "Basura/plásticos flotantes; valores altos = candidate debris"),
    "PI": ("PI (ratio NIR/RED)", "Proxy simple de materiales flotantes flotantes/vegetación"),
    "NDVI": ("NDVI (Rouse 1974)", "Vegetación; en humedales delimita manglar"),
    "NDMI": ("NDMI (Gao 1996)", "Humedad de vegetación; estrés hídrico en manglar"),
    "2BDA_chla": ("2BDA (chla)", "Clorofila-a empírico en aguas turbias costeras"),
}


def classify_floating_plastics(fdi, threshold=None, water_mask=None):
    """Mascara booleana de candidatos a escombros flotantes sobre el FDI."""
    if threshold is None:
        if water_mask is not None:
            vals = fdi[water_mask]
            threshold = float(np.nanpercentile(vals, 98.5)) if np.isfinite(vals).any() else 0
        else:
            threshold = float(np.nanpercentile(fdi, 98.5))
    mask = fdi > threshold
    if water_mask is not None:
        mask = mask & water_mask
    return mask, float(threshold)


def synthetic_lagoon_bands(n=240, seed=7):
    """Escena sintetica tipo S2 de una laguna costera con pluma urbana y debris flotantes.

    SOLO para demostración cuando no hay conectividad: valores 0-1 tipo reflectancia L2A.
    """
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:n, 0:n]
    cx, cy = n * 0.78, n * 0.68
    lagoon = ((xx - n * 0.45) ** 2 / (n * 0.42) ** 2 + (yy - n * 0.5) ** 2 / (n * 0.34) ** 2) < 1
    plume = np.exp(-(((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * (n * 0.09) ** 2))) * lagoon
    base = rng.normal(0.02, 0.004, (n, n))
    bands = {}
    bands["B3"] = base + lagoon * (0.09 + plume * 0.16) + (~lagoon) * 0.08
    bands["B4"] = base + lagoon * (0.07 + plume * 0.22) + (~lagoon) * 0.09
    bands["B5"] = base + lagoon * (0.085 + plume * 0.12) + (~lagoon) * 0.15
    bands["B8"] = base + lagoon * 0.028 + (~lagoon) * 0.24
    bands["B8A"] = base + lagoon * 0.03 + (~lagoon) * 0.22
    bands["B11"] = base + lagoon * (0.008 + plume * 0.02) + (~lagoon) * 0.13
    debris_rows = rng.integers(int(cy - 8), int(cy + 8), 40)
    debris_cols = rng.integers(int(cx - 14), int(cx + 14), 40)
    for ch, mult in [("B4", .35), ("B8", .30), ("B11", .22), ("B5", .30), ("B8A", .3)]:
        for r, c in zip(debris_rows, debris_cols):
            bands[ch][r:r + 2, c:c + 3] += mult
    return {"bands": bands, "lagoon_mask": lagoon, "plume_intensity": plume,
            "nota": "Escena SINTÉTICA tipo Sentinel-2 para modo demostración sin conexión"}


def fetch_sentinel2_items(bbox, date_range="2024-01-01/2025-12-31",
                          cloud_max=20, limit=10, api_key=None):
    """Busca items Sentinel-2 L2A en Microsoft Planetary Computer (STAC, gratis)."""
    import requests
    url = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
    body = {"collections": ["sentinel-2-l2a"],
            "bbox": [bbox["lon_min"], bbox["lat_min"], bbox["lon_max"], bbox["lat_max"]],
            "datetime": date_range,
            "query": {"eo:cloud_cover": {"lt": cloud_max}},
            "limit": limit}
    headers = {}
    if api_key:
        headers["Ocp-Apim-Subscription-Key"] = api_key
    try:
        r = requests.post(url, json=body, headers=headers, timeout=25)
        r.raise_for_status()
        items = r.json().get("features", [])
        out = []
        for it in items:
            props = it.get("properties", {})
            out.append({"id": it.get("id"), "fecha": props.get("datetime"),
                        "nubes_pct": round(props.get("eo:cloud_cover", -1), 1),
                        "escena": props.get("s2:mgrs_tile", ""),
                        "thumb": (it.get("assets", {}).get("rendered_preview", {}) or {}).get("href")})
        return {"ok": True, "items": out}
    except Exception as e:
        return {"ok": False, "error": str(e),
                "hint": "Sin conectividad a Planetary Computer; use la escena sintética de demostración"}


def raster_info(path):
    """Metadatos de un raster (requiere rasterio opcional)."""
    try:
        import rasterio
    except ImportError:
        return {"ok": False, "error": "rasterio no instalado (pip install rasterio)",
                "hint": "Disponible el modo sintético sin rasterio"}
    try:
        with rasterio.open(path) as src:
            return {"ok": True, "ancho": src.width, "alto": src.height,
                    "n_bandas": src.count, "crs": str(src.crs), "bounds": list(src.bounds),
                    "res": float(src.res[0])}
    except Exception as e:
        return {"ok": False, "error": str(e)}


def read_raster_bands(path, bands=(1, 2, 3, 4)):
    try:
        import rasterio
    except ImportError:
        return None
    try:
        with rasterio.open(path) as src:
            return {f"B{i}": src.read(i).astype(float) for i in bands}
    except Exception:
        return None
