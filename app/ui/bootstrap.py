"""Bootstrap de la UI: rutas, page config, tema, semilla y helpers compartidos."""
import sys
from pathlib import Path

_here = Path(__file__).resolve()
ROOT = _here.parents[2]
for _p in (str(ROOT), str(ROOT / "app")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import json

import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402
from streamlit.errors import StreamlitAPIException  # noqa: E402

from omniplast.config import APP_NAME, TAGLINE, VERSION  # noqa: E402
from omniplast.db.seed import seed_demo_if_empty  # noqa: E402
from ui.theme import CSS  # noqa: E402

st.set_page_config(page_title=f"{APP_NAME} · {TAGLINE}", page_icon="🌊",
                   layout="wide", initial_sidebar_state="expanded")


def _inyectar_tema():
    """Inyecta el CSS en el documento padre: Streamlit 1.59 vacía los <style> de
    st.markdown, por lo que se usa un componente que agrega la hoja a document.head."""
    try:
        import streamlit.components.v1 as components
        css_json = json.dumps(CSS)
        components.html(
            f"""
            <script>
            (function() {{
              var doc = window.parent.document;
              if (!doc.getElementById('omp-fonts')) {{
                var link = doc.createElement('link');
                link.id = 'omp-fonts';
                link.rel = 'stylesheet';
                link.href = 'https://fonts.googleapis.com/css2?family=Sora:wght@600;700;800&family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;600&display=swap';
                doc.head.appendChild(link);
              }}
              if (!doc.getElementById('omniplast-theme')) {{
                var style = doc.createElement('style');
                style.id = 'omniplast-theme';
                style.textContent = {css_json};
                doc.head.appendChild(style);
              }}
            }})();
            </script>
            """, height=0)
    except Exception:
        st.markdown(f"<style>{CSS}</style>", unsafe_allow_html=True)


_inyectar_tema()

pd.set_option("display.max_columns", 60)
pd.set_option("display.width", 220)

try:
    seed_demo_if_empty()
except Exception:
    pass


def chart(fig, **kwargs):
    try:
        st.plotly_chart(fig, width="stretch", **kwargs)
    except TypeError:
        try:
            st.plotly_chart(fig, use_container_width=True, **kwargs)
        except TypeError:
            st.plotly_chart(fig, **kwargs)


def tabla(df, height=None, **kwargs):
    if df is None:
        return
    if height is not None:
        kwargs.setdefault("height", height)
    try:
        st.dataframe(df, width="stretch", **kwargs)
    except (TypeError, StreamlitAPIException):
        try:
            st.dataframe(df, use_container_width=True, **kwargs)
        except (TypeError, StreamlitAPIException):
            st.dataframe(df, **kwargs)


def editor_seguro(df, **kwargs):
    try:
        return st.data_editor(df, width="stretch", **kwargs)
    except (TypeError, StreamlitAPIException):
        try:
            return st.data_editor(df, use_container_width=True, **kwargs)
        except (TypeError, StreamlitAPIException):
            return st.data_editor(df, **kwargs)


def page_header(titulo, subtitulo=""):
    sub = f'<div class="sub">{subtitulo}</div>' if subtitulo else ""
    st.markdown(
        f'<div class="omp-header"><h1>{titulo}</h1>{sub}</div>', unsafe_allow_html=True)


def badge(texto, tipo="info"):
    return f'<span class="omp-badge {tipo}">{texto}</span>'


def kpi_row(items):
    """items: [(valor, etiqueta, sub, color_css)]"""
    cells = "".join(
        f'<div class="omp-card"><div class="kpi-label">{etq}</div>'
        f'<div class="kpi-value" style="color:{color}">{val}</div>'
        f'<div class="kpi-sub">{sub}</div></div>'
        for val, etq, sub, color in items)
    st.markdown(f'<div class="omp-grid">{cells}</div>', unsafe_allow_html=True)


def section_title(texto):
    st.markdown(f'<div class="omp-section-title"><span class="bar"></span>{texto}</div>',
                unsafe_allow_html=True)


def resumen_box(texto, tipo="ok"):
    cls = {"ok": "omp-resumen", "warn": "omp-requiere", "error": "omp-error"}.get(tipo, "omp-resumen")
    st.markdown(f'<div class="{cls}">{texto}</div>', unsafe_allow_html=True)


def render_result(result):
    """Renderiza el resultado estandar de un agente o analisis."""
    estado = result.get("estado", "ok")
    if estado == "error":
        resumen_box(f"<b>Error:</b> {result.get('resumen', '')}", "error")
    elif estado == "requiere":
        resumen_box(f"<b>Requiere insumos:</b> {result.get('resumen', '')}", "warn")
    else:
        resumen_box(result.get("resumen", ""))
    tabs_names = []
    if result.get("tablas"):
        tabs_names.append("Tablas")
    if result.get("figuras"):
        tabs_names.append("Figuras")
    if result.get("pasos"):
        tabs_names.append("Siguientes pasos")
    if result.get("refs"):
        tabs_names.append("Referencias")
    if not tabs_names:
        return
    tabs = st.tabs(tabs_names)
    i = 0
    if result.get("tablas"):
        with tabs[i]:
            for nombre, df in result["tablas"].items():
                st.caption(nombre)
                tabla(df)
        i += 1
    if result.get("figuras"):
        with tabs[i]:
            for nombre, fig in result["figuras"].items():
                st.caption(nombre)
                try:
                    chart(fig)
                except Exception as e:
                    st.warning(f"Figura '{nombre}' no renderizable: {e}")
        i += 1
    if result.get("pasos"):
        with tabs[i]:
            for paso in result["pasos"]:
                st.markdown(f"- {paso}")
        i += 1
    if result.get("refs"):
        with tabs[i]:
            st.markdown('<div class="omp-refs">' +
                        "".join(f"<span>· {r}</span>" for r in result["refs"]) +
                        "</div>", unsafe_allow_html=True)


def footer():
    st.markdown(
        '<div class="omp-footer"><b>OMNIPLAST</b> v' + VERSION + ' · Plataforma doctoral de '
        'análisis integral de microplásticos · Proyecto demo con <b>DATOS SINTÉTICOS</b> '
        'etiquetados (Ciénaga de La Virgen, Cartagena) · Índices y umbrales referenciales '
        'citados en docs/METODOLOGIA.md · API REST: puerto 8010/docs</div>',
        unsafe_allow_html=True)


def selector_matriz(df, key="matriz"):
    opciones = sorted(df["matriz"].unique()) if "matriz" in df.columns else []
    return st.selectbox("Matriz", opciones, key=key) if opciones else None


def exportar_dfs(nombre, dfs: dict):
    from omniplast.config import EXPORT_DIR
    EXPORT_DIR.mkdir(parents=True, exist_ok=True)
    import io
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as w:
        for hoja, df in dfs.items():
            if isinstance(df, pd.DataFrame) and not df.empty:
                df.to_excel(w, sheet_name=str(hoja)[:31], index=False)
    st.download_button("Descargar resultados (Excel)", buf.getvalue(),
                       file_name=f"{nombre}.xlsx", mime="application/vnd.ms-excel")
