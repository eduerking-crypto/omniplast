"""Componentes visuales compartidos: mapas folium y tarjetas de agente."""
import branca
import folium
import pandas as pd

from omniplast.config import PALETTE

ZONA_COLORS = {"impacto": "#FF6B6B", "transicion": "#E9B44C",
               "referencia_interna": "#05BFDB", "referencia": "#3E9B4F"}


def mapa_estaciones(per_sample, matriz="agua"):
    """Mapa folium con círculos proporcionales a la mediana de MP por estación."""
    sites = per_sample.groupby(["sitio", "zona"], as_index=False, observed=True).agg(
        lat=("lat", "first") if "lat" in per_sample.columns else ("abundancia", "first"),
        lon=("lon", "first") if "lon" in per_sample.columns else ("abundancia", "first"),
        abundancia=("abundancia", "median"))
    if "lat" not in per_sample.columns or sites.empty:
        return None
    center = [sites["lat"].mean(), sites["lon"].mean()]
    m = folium.Map(location=center, zoom_start=12, tiles="CartoDB positron")
    vmax = max(sites["abundancia"].max(), 1)
    for _, r in sites.iterrows():
        radius = 8 + 22 * (r["abundancia"] / vmax)
        color = ZONA_COLORS.get(r["zona"], PALETTE["teal"])
        folium.CircleMarker(
            location=(r["lat"], r["lon"]), radius=radius,
            popup=folium.Popup(
                f"<b>{r['sitio']}</b><br>zona: {r['zona']}<br>MP mediana: "
                f"{r['abundancia']:,.0f} ({matriz})", max_width=220),
            color=color, weight=2, fill=True, fill_opacity=.45, fill_color=color).add_to(m)
    legend_html = branca.element.Element("""
    <div style="position: fixed; bottom: 22px; left: 22px; z-index: 9999;
                background: white; padding: 9px 13px; border-radius: 9px;
                border: 1px solid #E1E8F0; font: 600 12px Inter, sans-serif;">
      <div style="margin-bottom:4px; color:#0A2540;">Zonas</div>
      <span style="color:#FF6B6B">●</span> impacto &nbsp;
      <span style="color:#E9B44C">●</span> transición<br>
      <span style="color:#05BFDB">●</span> ref. interna &nbsp;
      <span style="color:#3E9B4F">●</span> referencia
    </div>""")
    m.get_root().add_child(legend_html)
    return m


def agente_card(agente, key_idx, on_run=None):
    """Renderiza la tarjeta de un agente. Devuelve True si se pulsó ejecutar."""
    import streamlit as st
    nivel_badge = {"Basico": "info", "Intermedio": "ok", "Doctoral": "warn"}.get(
        agente["nivel"], "info")
    st.markdown(
        f"""<div class="omp-agent" style="border-left-color:{_color(agente)}">
        <div class="nombre">{agente['id']} · {agente['nombre']}
        &nbsp;<span class="omp-badge {nivel_badge}">{agente['nivel']}</span></div>
        <div class="desc">{agente['desc']}</div></div>""", unsafe_allow_html=True)
    col1, col2, _ = st.columns([1, 1, 3])
    ejecutar = col1.button("Ejecutar", key=f"run_{agente['id']}_{key_idx}", type="primary",
                           use_container_width=True)
    return ejecutar


def _color(agente):
    from omniplast.agents.registry import CATEGORIES
    return CATEGORIES.get(agente["categoria"], ("", "#00A896"))[1]
