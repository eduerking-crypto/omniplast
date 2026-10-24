import sys
from pathlib import Path

_here = Path(__file__).resolve()
for _p in _here.parents:
    if (_p / "app" / "ui" / "bootstrap.py").exists():
        for _s in (str(_p), str(_p / "app")):
            if _s not in sys.path:
                sys.path.insert(0, _s)
        break

import numpy as np
import pandas as pd
import streamlit as st
from ui.bootstrap import chart, footer, page_header, resumen_box, section_title, tabla

from omniplast.core.geo import remote_sensing
from omniplast.db import dao
from omniplast.reporting import figures

page_header("Teledetección satelital",
            "Sentinel-2: índices acuáticos, turbidez, clorofila y detección de plásticos flotantes (FDI)")

t_idx, t_fdi, t_stac = st.tabs(["Índices espectrales", "Plásticos flotantes (FDI)",
                                "Escenas Sentinel-2 reales (STAC)"])

upload = st.file_uploader("Raster propio (GeoTIFF multibanda B3,B4,B5,B8,B11 — opcional)",
                          type=["tif", "tiff"])
bands, mask, nota = None, None, "ESCENA SINTÉTICA DE DEMOSTRACIÓN"
if upload is not None:
    rb = remote_sensing.read_raster_bands(upload, bands=(2, 3, 4, 8, 12))
    if rb:
        rename = {list(rb.keys())[0]: "B3", list(rb.keys())[1]: "B4", list(rb.keys())[2]: "B5",
                 list(rb.keys())[3]: "B8", list(rb.keys())[4]: "B11"}
        bands = {rename[k]: v for k, v in rb.items()}
        mask = None
        nota = "RASTER DEL USUARIO"

with t_idx:
    if bands is None:
        escena = remote_sensing.synthetic_lagoon_bands()
        bands, mask = escena["bands"], escena["lagoon_mask"]
    sel = st.multiselect("Índices", list(remote_sensing.INDEX_INFO.keys()),
                         default=["NDWI", "NDTI", "NDCI"])
    vals = []
    for name in sel:
        try:
            arr = remote_sensing.compute_index(bands, name)
        except KeyError as e:
            resumen_box(f"El índice {name} requiere la banda {e} (no disponible).", "warn")
            continue
        media = float(np.nanmean(arr[mask])) if mask is not None else float(np.nanmean(arr))
        label, desc = remote_sensing.INDEX_INFO[name]
        vals.append({"indice": name, "nombre": label, "media": round(media, 4),
                     "interpretación": desc})
        chart(figures.heatmap(arr, title=f"{label} — {nota}", colorscale="Teal"))
    tabla(pd.DataFrame(vals))

with t_fdi:
    if bands is None:
        escena = remote_sensing.synthetic_lagoon_bands()
        bands, mask = escena["bands"], escena["lagoon_mask"]
    try:
        fdi = remote_sensing.compute_index(bands, "FDI")
        mask_fdi, thr = remote_sensing.classify_floating_plastics(fdi, water_mask=mask)
        agua_px = int(mask.sum()) if mask is not None else fdi.size
        pct = float(mask_fdi.sum()) / agua_px * 100 if agua_px else 0
        chart(figures.heatmap(fdi, title=f"FDI — candidatos a escombros flotantes ({nota})",
                             colorscale="Viridis"))
        resumen_box(f"Umbral P98.5 = {thr:.4f} · <b>{pct:.2f}%</b> del área de agua con "
                    "firma flotante. El FDI (Biermann et al. 2020) detecta 'floating debris' "
                    "GENÉRICO: valide en campo/dron antes de atribuirlo a plástico.",
                    "warn")
        tabla(pd.DataFrame([{"umbral": round(thr, 4), "candidatos_pct": round(pct, 3),
                             "pixeles_agua": agua_px}]))
    except KeyError as e:
        resumen_box(f"FDI requiere bandas B4/B8/B11 (falta {e}).", "error")

with t_stac:
    from omniplast.config import DEMO_HULL
    st.markdown("Búsqueda real de escenas **Sentinel-2 L2A** (Microsoft Planetary Computer, "
                "gratuito) sobre la Ciénaga de La Virgen. El área y el rango son "
                "editables desde la API (`/api/conectores`).")
    if st.button("Buscar escenas 2024-2025 (nubes < 20%)", type="primary"):
        r = remote_sensing.fetch_sentinel2_items(DEMO_HULL, limit=10)
        if r["ok"]:
            items = pd.DataFrame(r["items"])
            tabla(items)
            st.success(f"{len(items)} escenas encontradas. La descarga de píxeles "
                       "requiere token de Planetary Computer (registro gratuito).")
        else:
            resumen_box(f"Sin acceso al catálogo: {r['error']}. {r.get('hint', '')}", "warn")
    st.info("Para trabajar con píxeles reales: `pip install pystac-client planetary-computer` "
            "y use `stackstac` u `odc-stac`; los scripts generados en la página Reportes "
            "incluyen la estructura recomendada.")

footer()
