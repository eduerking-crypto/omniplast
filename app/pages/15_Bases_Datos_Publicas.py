import sys
from pathlib import Path

_here = Path(__file__).resolve()
for _p in _here.parents:
    if (_p / "app" / "ui" / "bootstrap.py").exists():
        for _s in (str(_p), str(_p / "app")):
            if _s not in sys.path:
                sys.path.insert(0, _s)
        break

import pandas as pd
import streamlit as st
from ui.bootstrap import chart, footer, page_header, resumen_box, section_title, tabla

from omniplast.core.connectors import connectors
from omniplast.db import dao

page_header("Bases de datos públicas",
            "Biodiversidad (OBIS · GBIF · iNaturalist) · literatura (Crossref · PubMed) · "
            "clima (NASA POWER) · aire (OpenAQ) — todas gratuitas")

st.markdown('<div class="omp-resumen">Los conectores fallan <b>explícitamente</b> sin '
            'conexión: nunca se muestran registros inventados. OBIS/GBIF/Crossref son CC-BY '
            'y gratuitos; OpenAQ v3 requiere clave gratuita (variable de entorno '
            '<code>OPENAQ_API_KEY</code>) como ejemplo de escalado a APIs de pago.</div>',
            unsafe_allow_html=True)

fuente = st.selectbox("Conector", ["OBIS (biodiversidad marina)",
                                   "GBIF (todas las especies)",
                                   "iNaturalist (cercano al sitio)",
                                   "Crossref (artículos)",
                                   "PubMed (biomedicina)",
                                   "NASA POWER (clima del sitio)",
                                   "OpenAQ (calidad del aire)"])
q = st.text_input("Consulta (taxón, palabras clave…)", "Crassostrea rhizophorae")

if st.button("Consultar", type="primary"):
    with st.spinner("Consultando el servicio…"):
        if fuente.startswith("OBIS"):
            r = connectors.search_obis(scientific_name=q, limit=200)
        elif fuente.startswith("GBIF"):
            r = connectors.search_gbif(scientific_name=q, limit=200)
        elif fuente.startswith("iNaturalist"):
            r = connectors.search_inaturalist(taxon_name=q, lat=10.43, lon=-75.49,
                                               radius_km=8)
        elif fuente.startswith("Crossref"):
            r = connectors.search_crossref(q, rows=25)
        elif fuente.startswith("PubMed"):
            r = connectors.search_pubmed(q, rows=25)
        elif fuente.startswith("NASA POWER"):
            r = connectors.fetch_nasa_power(10.43, -75.49)
        else:
            r = connectors.fetch_openaq(lat=10.39, lon=-75.51)
    if not r.get("ok"):
        resumen_box(f"<b>{r.get('error', 'fallo')}</b> — {r.get('hint', '')}", "error")
    else:
        df = r.get("df")
        st.success(f"{r['n']} registros · {r['meta']['fuente']} · "
                   f"licencia: {r['meta'].get('licencia', '—')}")
        if isinstance(df, pd.DataFrame) and not df.empty:
            tabla(df)
            if {"lat", "lon"}.issubset(df.columns) and df[["lat", "lon"]].notna().all().all():
                section_title("Mapa de ocurrencias")
                import folium
                from streamlit_folium import st_folium
                m = folium.Map(location=[df["lat"].mean(), df["lon"].mean()],
                               zoom_start=9, tiles="CartoDB positron")
                for _, row in df.head(200).iterrows():
                    folium.CircleMarker((row["lat"], row["lon"]), radius=4,
                                        color="#00A896", fill=True, fill_opacity=.6,
                                        popup=str(row.get("especie")
                                                  or row.get("scientificName")
                                                  or row.get("titulo", ""))).add_to(m)
                st_folium(m, height=380, returned_objects=[])
        else:
            st.info("Respuesta sin tabla (revise el JSON crudo en la API).")

section_title("Guardado en la base local (literatura)")
if st.button("Guardar resultados de Crossref como referencias del proyecto"):
    r = connectors.search_crossref("microplastics OR microplastic", rows=20, from_year=2019)
    if r.get("ok") and len(r.get("df", [])):
        refs = []
        for row in r["df"].itertuples():
            refs.append({"titulo": row.titulo, "autores": row.autores,
                         "anio": row.anio, "revista": row.revista,
                         "doi": row.doi, "fuente": "Crossref",
                         "resumen": "", "etiquetas": q})
        dao.save_literature(refs)
        st.success(f"{len(refs)} referencias guardadas (página Literatura).")
    else:
        resumen_box(r.get("error", "sin resultados"), "error")

footer()
