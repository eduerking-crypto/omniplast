"""Conectores a bases de datos publicas (todas gratuitas).

Cada conector: manejo de errores y cache local (1 h) en data/cache. Nunca se inventan
registros: si falla la red se devuelve ok=False con orientacion. Claves de API de pago
(OpenAQ v3, Copernicus) se leen de variables de entorno para el futuro.
"""
import json
import os
import time
from pathlib import Path

import pandas as pd
import requests

from ...config import CACHE_DIR, CONNECTORS

CACHE_TTL = 3600
UA = {"User-Agent": "OMNIPLAST/1.0 (plataforma académica de microplásticos)"}


def _cache_path(key):
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in key)
    return CACHE_DIR / f"{safe[:180]}.json"


def fetch_json(url, params=None, headers=None, timeout=20, cache_key=None, method="get",
               body=None):
    key = cache_key or f"{url}|{json.dumps(params or body, sort_keys=True)}"
    cpath = _cache_path(key)
    if cpath.exists() and time.time() - cpath.stat().st_mtime < CACHE_TTL:
        return {"ok": True, "data": json.loads(cpath.read_text(encoding="utf-8")), "cache": True}
    try:
        h = {**UA, **(headers or {})}
        if method == "get":
            r = requests.get(url, params=params, headers=h, timeout=timeout)
        else:
            r = requests.post(url, json=body, headers=h, timeout=timeout)
        r.raise_for_status()
        data = r.json()
        cpath.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
        return {"ok": True, "data": data, "cache": False}
    except Exception as e:
        return {"ok": False, "error": str(e)}


def search_obis(scientific_name=None, bbox=None, year=None, limit=200):
    params = {"limit": min(limit, 1000)}
    if scientific_name:
        params["scientificname"] = scientific_name
    if year:
        params["year"] = year
    if bbox:
        params["geometry"] = (f"POLYGON(({bbox['lon_min']} {bbox['lat_min']},"
                              f"{bbox['lon_max']} {bbox['lat_min']},"
                              f"{bbox['lon_max']} {bbox['lat_max']},"
                              f"{bbox['lon_min']} {bbox['lat_max']},"
                              f"{bbox['lon_min']} {bbox['lat_min']}))")
    res = fetch_json(CONNECTORS["obis"], params=params)
    if not res["ok"]:
        return {"ok": False, "error": res["error"],
                "hint": "OBIS (api.obis.org) requiere conexión a internet"}
    records = res["data"].get("results", [])
    df = pd.DataFrame(records)
    if not df.empty:
        keep = [c for c in ["scientificName", "decimalLongitude", "decimalLatitude",
                            "eventDate", "basisOfRecord", "occurrenceStatus",
                            "marineRegion", "obisID", "individualCount"] if c in df.columns]
        df = df[keep]
    return {"ok": True, "df": df, "n": len(df),
            "meta": {"fuente": "OBIS (Ocean Biodiversity Information System)",
                     "licencia": "CC-BY 4.0", "url": "https://obis.org"}}


def search_gbif(scientific_name=None, bbox=None, limit=200, country=None):
    params = {"limit": min(limit, 300), "hasCoordinate": "true"}
    if scientific_name:
        params["scientificName"] = scientific_name
    if country:
        params["country"] = country
    if bbox:
        params["decimalLatitude"] = f"{bbox['lat_min']},{bbox['lat_max']}"
        params["decimalLongitude"] = f"{bbox['lon_min']},{bbox['lon_max']}"
    res = fetch_json(CONNECTORS["gbif"], params=params)
    if not res["ok"]:
        return {"ok": False, "error": res["error"],
                "hint": "GBIF (api.gbif.org) requiere conexión a internet"}
    records = res["data"].get("results", [])
    df = pd.DataFrame(records)
    if not df.empty:
        keep = [c for c in ["scientificName", "decimalLongitude", "decimalLatitude",
                            "eventDate", "basisOfRecord", "country",
                            "datasetName", "occurrenceID"] if c in df.columns]
        df = df[keep]
    return {"ok": True, "df": df, "n": len(df),
            "meta": {"fuente": "GBIF", "licencia": "CC-BY 4.0", "url": "https://gbif.org"}}


def worms_lookup(scientific_name):
    res = fetch_json(f"{CONNECTORS['worms']}?scientificnames[]={requests.utils.quote(scientific_name)}&marine_only=false")
    if not res["ok"]:
        return {"ok": False, "error": res["error"],
                "hint": "WoRMS (marinespecies.org) requiere conexión"}
    data = res["data"]
    if isinstance(data, list) and data:
        rec = data[0][0] if isinstance(data[0], list) else data[0]
        return {"ok": True,
                "record": {"AphiaID": rec.get("AphiaID"), "nombre_valido": rec.get("valid_name"),
                           "autoridad": rec.get("valid_authority"), "estado": rec.get("status"),
                           "clase": rec.get("class"), "orden": rec.get("order"),
                           "familia": rec.get("family"), "ISIN": rec.get("isMarine"),
                           "url": f"https://www.marinespecies.org/aphia.php?p=taxdetails&id={rec.get('AphiaID')}"}}
    return {"ok": False, "error": "Sin coincidencias en WoRMS (¿especie marina?)"}


def search_inaturalist(taxon_name=None, lat=None, lon=None, radius_km=5, limit=100):
    params = {"per_page": min(limit, 200), "photos": "true", "verifiable": "true"}
    if taxon_name:
        params["taxon_name"] = taxon_name
    if lat is not None and lon is not None:
        params["lat"], params["lng"] = lat, lon
        params["radius"] = radius_km
    res = fetch_json(CONNECTORS["inaturalist"], params=params)
    if not res["ok"]:
        return {"ok": False, "error": res["error"],
                "hint": "iNaturalist (api.inaturalist.org) requiere conexión"}
    results = res["data"].get("results", [])
    rows = []
    for obs in results:
        rows.append({"especie": obs.get("species_guess"),
                     "lat": obs.get("latitude"), "lon": obs.get("longitude"),
                     "fecha": obs.get("observed_on"), "lugar": obs.get("place_guess"),
                     "n_obs": obs.get("num_identification_agreements", 0),
                     "imagen": ((obs.get("photos") or [{}])[0].get("url", "")).replace("square", "medium")})
    return {"ok": True, "df": pd.DataFrame(rows), "n": len(rows),
            "meta": {"fuente": "iNaturalist", "licencia": "CC-BY-NC-4.0 (por observación)",
                     "url": "https://inaturalist.org"}}


def search_crossref(query, rows=20, from_year=None):
    params = {"query.bibliographic": query, "rows": min(rows, 100),
              "select": "DOI,title,author,issued,container-title,is-referenced-by-count,type,URL"}
    if from_year:
        params["from-pub-date"] = f"{from_year}-01-01"
    res = fetch_json(CONNECTORS["crossref"], params=params)
    if not res["ok"]:
        return {"ok": False, "error": res["error"],
                "hint": "Crossref (api.crossref.org) requiere conexión"}
    items = res["data"].get("message", {}).get("items", [])
    rows = []
    for it in items:
        title = (it.get("title") or [""])[0]
        authors = "; ".join(f"{a.get('given', '')} {a.get('family', '')}".strip()
                           for a in (it.get("author") or [])[:4])
        rows.append({"titulo": title, "autores": authors,
                     "anio": (it.get("issued", {}).get("date-parts") or [[None]])[0][0],
                     "revista": (it.get("container-title") or [""])[0],
                     "doi": it.get("DOI", ""), "citas": it.get("is-referenced-by-count", 0),
                     "tipo": it.get("type", ""), "url": it.get("URL", "")})
    df = pd.DataFrame(rows).sort_values("citas", ascending=False)
    return {"ok": True, "df": df, "n": len(df),
            "meta": {"fuente": "Crossref", "licencia": "CC0 metadata", "url": "https://crossref.org"}}


def search_pubmed(term, rows=20):
    res = fetch_json(f"{CONNECTORS['pubmed']}esearch.fcgi",
                     params={"db": "pubmed", "term": term, "retmode": "json", "retmax": rows})
    if not res["ok"]:
        return {"ok": False, "error": res["error"],
                "hint": "NCBI E-utilities requiere conexión"}
    ids = res["data"].get("esearchresult", {}).get("idlist", [])
    if not ids:
        return {"ok": True, "df": pd.DataFrame(), "n": 0, "meta": {"fuente": "PubMed"}}
    res2 = fetch_json(f"{CONNECTORS['pubmed']}esummary.fcgi",
                      params={"db": "pubmed", "id": ",".join(ids), "retmode": "json"})
    if not res2["ok"]:
        return {"ok": False, "error": res2["error"]}
    rows = []
    for pid in ids:
        rec = res2["data"].get("result", {}).get(pid, {})
        rows.append({"titulo": rec.get("title", ""), "revista": rec.get("fulljournalname", ""),
                     "anio": rec.get("pubdate", "")[:4], "pmid": pid,
                     "url": f"https://pubmed.ncbi.nlm.nih.gov/{pid}/"})
    return {"ok": True, "df": pd.DataFrame(rows), "n": len(rows),
            "meta": {"fuente": "PubMed/NCBI", "licencia": "dominio público (EE. UU.)"}}


def fetch_openaq(lat=None, lon=None, city=None):
    api_key = os.environ.get("OPENAQ_API_KEY")
    params = {"limit": 100}
    if lat is not None and lon is not None:
        params["coordinates"] = f"{lat},{lon}"
        params["radius"] = 25000
    if city:
        params["city"] = city
    headers = {"X-API-Key": api_key} if api_key else {}
    res = fetch_json(CONNECTORS["openaq"], params=params, headers=headers)
    if not res["ok"]:
        return {"ok": False, "error": res["error"],
                "hint": "OpenAQ v3 requiere clave gratuita: registre una en openaq.org y "
                        "exportela como OPENAQ_API_KEY (los datos de calidad del aire "
                        "contextualizan fuentes urbanas de MP)"}
    results = res["data"].get("results", [])
    rows = []
    for m in results:
        rows.append({"ubicacion": m.get("location", ""),
                     "ciudad": (m.get("city") or ""),
                     "parametro": (m.get("parameter") or m.get("measurements", [{}])[0].get("parameter", "")),
                     "valor": ((m.get("measurements") or [{}])[0].get("lastValue")),
                     "fecha": (m.get("measurements") or [{}])[0].get("lastUpdated")})
    return {"ok": True, "df": pd.DataFrame(rows), "n": len(rows),
            "meta": {"fuente": "OpenAQ", "licencia": "CC-BY 4.0", "url": "https://openaq.org"}}


def fetch_nasa_power(lat, lon, start="20240101", end="20251231"):
    params = {"parameters": "PRECTOTCORR,TS,T2M,RH2M,WS2M",
              "community": "RE", "latitude": lat, "longitude": lon,
              "start": start, "end": end, "format": "JSON"}
    res = fetch_json(CONNECTORS["nasa_power"], params=params, timeout=40,
                    cache_key=f"power_{lat}_{lon}_{start}_{end}")
    if not res["ok"]:
        return {"ok": False, "error": res["error"],
                "hint": "NASA POWER (power.larc.nasa.gov) requiere conexión"}
    props = res["data"].get("properties", {}).get("parameter", {})
    df = pd.DataFrame(props)
    if not df.empty:
        df.index = pd.to_datetime(df.index, format="%Y%m%d", errors="coerce")
        df = df.dropna(how="all")
    return {"ok": True, "df": df, "n": len(df),
            "meta": {"fuente": "NASA POWER",
                     "variables": {"PRECTOTCORR": "precipitación (mm/día)",
                                   "TS": "temperatura de suelo (K)", "T2M": "T aire 2m (K)",
                                   "RH2M": "humedad relativa (%)", "WS2M": "viento 2m (m/s)"},
                     "licencia": "dominio público NASA",
                     "nota": "Use estos datos para modelar resuspensión y descarga pluvial de MP"}}
