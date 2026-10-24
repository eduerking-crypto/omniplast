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
from ui.components import mapa_estaciones

from omniplast.core.geo import spatial, trajectory
from omniplast.db import dao
from omniplast.reporting import figures

page_header("GIS, geoestadística y dispersión",
            "Kriging/IDW · Moran/LISA · hotspots Gi* · modelo lagrangiano de partículas")

mp = dao.mp_df()
if mp.empty:
    st.stop()

mpa = mp[mp["matriz"] == "agua"].groupby(["sitio", "lat", "lon", "zona"],
                                         observed=True)["abundancia"].mean().reset_index()

t_mapa, t_interp, t_auto, t_drift = st.tabs(
    ["Mapa de estaciones", "Interpolación espacial", "Autocorrelación y hotspots",
     "Dispersión de partículas"])

with t_mapa:
    fig = mapa_estaciones(mpa.assign(abundancia=mpa["abundancia"]), matriz="agua")
    if fig is not None:
        from streamlit_folium import st_folium
        st_folium(fig, height=430, returned_objects=[])
    tabla(mpa.round(2))

with t_interp:
    if len(mpa) >= 5:
        from omniplast.config import DEMO_HULL
        grid_lon, grid_lat = spatial.build_grid(DEMO_HULL, n=36)
        metodo = st.radio("Método", ["IDW (rápido)", "Kriging ordinario"], horizontal=True)
        if metodo.startswith("IDW"):
            grid = spatial.idw_interpolate(mpa["lat"], mpa["lon"], mpa["abundancia"],
                                           grid_lon, grid_lat, power=2)
            resumen_box("IDW con potencia 2. Para publicación use kriging con "
                        "validación cruzada leave-one-out.")
        else:
            grid, fit = spatial.ordinary_kriging(mpa["lat"], mpa["lon"], mpa["abundancia"])
            if fit:
                resumen_box(f"Variograma {fit['modelo']}: nugget={fit['nugget']:.3g}, "
                            f"meseta={fit['meseta']:.3g}, alcance={fit['alcance_km']:.2f} km, "
                            f"R²={fit['r2']:.2f}.")
                emp, _ = spatial.variogram_empirical(mpa["lat"], mpa["lon"], mpa["abundancia"])
                chart(figures.scatter(emp, x="dist_km", y="semivarianza",
                                      title="Variograma empírico", xtitle="distancia (km)"))
            else:
                resumen_box("Variograma no ajustable con estos puntos; usando IDW.", "warn")
                grid = spatial.idw_interpolate(mpa["lat"], mpa["lon"], mpa["abundancia"],
                                              grid_lon, grid_lat)
        chart(figures.heatmap(grid, title="Superficie de concentración de MP (part/m³)",
                             colorscale="Teal"))

with t_auto:
    if len(mpa) >= 5:
        mi = spatial.morans_i(mpa["lat"], mpa["lon"], mpa["abundancia"], k=8)
        resumen_box(f"<b>Moran global:</b> {mi.get('interpretacion', mi.get('error', ''))}")
        lisa = spatial.local_moran(mpa["lat"], mpa["lon"], mpa["abundancia"], k=8)
        lisa.insert(0, "sitio", mpa["sitio"].values)
        tabla(lisa.round(3))
        gi = spatial.getis_ord_gi_star(mpa["lat"], mpa["lon"], mpa["abundancia"], k=6)
        gi.insert(0, "sitio", mpa["sitio"].values)
        section_title("Hotspots Getis-Ord Gi*")
        chart(figures.bars(gi.sort_values("gi_z", ascending=False), "sitio", "gi_z",
                          title="z de Gi* (hotspots >1.645)", ytitle="z"))
        tabla(gi.round(3))

with t_drift:
    st.markdown("Modelo lagrangiano: advección (viento + corriente) + difusión turbulenta "
                "+ hundimiento por Stokes (densidad del polímero). Origen: estación más "
                "contaminada o coordenadas manuales.")
    c1, c2, c3 = st.columns(3)
    lat_o = c1.number_input("Lat origen", 10.40, 10.47, 10.4238, step=0.0001, format="%.4f")
    lon_o = c2.number_input("Lon origen", -75.52, -75.46, -75.5091, step=0.0001, format="%.4f")
    horas = c3.selectbox("Horas de simulación", [6, 12, 24, 48], index=2)
    c4, c5, c6 = st.columns(3)
    polimero = c4.selectbox("Polímero", ["PE (0.91 flotante)", "PP (0.91 flotante)",
                                         "PET (1.38 sedimenta)", "PVC (1.40 sedimenta)"])
    viento = c5.slider("Viento (m/s)", 0.5, 12.0, 3.5, 0.5)
    corriente = c6.slider("Corriente (m/s)", 0.0, 0.5, 0.12, 0.01)
    densidad = {"PE": 960, "PP": 960, "PET": 1380, "PVC": 1400}[polimero.split()[0]]
    from omniplast.config import DEMO_HULL
    if st.button("Simular dispersión (1500 partículas)", type="primary"):
        r = trajectory.simulate_drift(lat_o, lon_o, hours=int(horas), n_particles=1500,
                                      bbox=DEMO_HULL, wind_speed_ms=float(viento),
                                      current_speed_ms=float(corriente),
                                      polymer_density=float(densidad))
        resumen_box(r["interpretacion"])
        pos = r["posiciones_finales"]
        chart(figures.scatter(pos, x="lon", y="lat", color="estado",
                              title=f"Posiciones tras {horas} h", xtitle="longitud",
                              ytitle="latitud"))
        st.caption("Nota: modelo didáctico de orden de magnitud; para estudios de "
                   "transporte formal use ROMS/OpenDrift con batimetría real.")

footer()
