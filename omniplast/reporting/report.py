"""Generador de informes academicos DOCX con figuras, tablas y referencias."""
from datetime import datetime
from pathlib import Path

import pandas as pd

from ..config import APP_NAME, EXPORT_DIR, TAGLINE, VERSION


def _fig_to_image_bytes(fig):
    try:
        return fig.to_image(format="png", width=900, scale=1.6)
    except Exception:
        return None


def _style_doc(doc):
    from docx.shared import Pt, RGBColor
    style = doc.styles["Normal"]
    style.font.name = "Calibri"
    style.font.size = Pt(10.5)
    for i in range(1, 4):
        h = doc.styles[f"Heading {i}"]
        h.font.color.rgb = RGBColor(0x0A, 0x25, 0x40)
    return doc


def _add_table(doc, df, max_rows=25):
    from docx.shared import Pt, RGBColor
    df = df.head(max_rows)
    df = df.copy()
    df = df.astype(object).where(df.notna(), "—")
    table = doc.add_table(rows=1, cols=len(df.columns))
    table.style = "Light Grid Accent 1"
    for j, col in enumerate(df.columns):
        cell = table.rows[0].cells[j]
        cell.text = str(col)
        for p in cell.paragraphs:
            for r in p.runs:
                r.font.bold = True
                r.font.size = Pt(8.5)
    for _, row in df.iterrows():
        cells = table.add_row().cells
        for j, v in enumerate(row):
            if isinstance(v, float):
                v = f"{v:,.4g}"
            cells[j].text = str(v)[:90]
            for p in cells[j].paragraphs:
                for r in p.runs:
                    r.font.size = Pt(8.5)
    return table


def generate_docx(titulo, secciones, filename=None, meta=None):
    """secciones: lista de dicts {titulo, parrafos=[...], tablas={nombre: df},
    figuras={nombre: plotly.Fig}, notas=[...]}."""
    from docx import Document
    from docx.shared import Inches, Pt

    doc = Document()
    _style_doc(doc)
    meta = meta or {}
    doc.add_heading(f"{APP_NAME} — {TAGLINE}", level=3)
    doc.add_heading(titulo, level=1)
    sub = (f"Proyecto: {meta.get('proyecto', '—')} · Región: {meta.get('region', '—')} · "
           f"Generado: {datetime.now():%Y-%m-%d %H:%M} · v{VERSION}")
    p = doc.add_paragraph(sub)
    p.runs[0].font.size = Pt(9)
    if meta.get("advertencia"):
        doc.add_paragraph(meta["advertencia"])

    for sec in secciones:
        doc.add_heading(sec.get("titulo", "Sección"), level=2)
        for par in sec.get("parrafos", []):
            doc.add_paragraph(par)
        for nombre, df in sec.get("tablas", {}).items():
            if isinstance(df, pd.DataFrame) and not df.empty:
                doc.add_paragraph(nombre).runs[0].font.bold = True
                _add_table(doc, df)
        for nombre, fig in sec.get("figuras", {}).items():
            img = _fig_to_image_bytes(fig)
            if img:
                doc.add_paragraph(nombre).runs[0].font.bold = True
                tmp = EXPORT_DIR / f"_tmp_{nombre.replace(' ', '_')}.png"
                tmp.write_bytes(img)
                doc.add_picture(str(tmp), width=Inches(6.2))
                tmp.unlink(missing_ok=True)
        for nota in sec.get("notas", []):
            p = doc.add_paragraph(f"Nota: {nota}")
            for r in p.runs:
                r.font.size = Pt(8.5)
                r.font.italic = True

    refs = meta.get("referencias", [])
    if refs:
        doc.add_heading("Referencias citadas por el sistema", level=2)
        for ref in refs:
            p = doc.add_paragraph(ref, style="List Bullet")
            for r in p.runs:
                r.font.size = Pt(9)
    doc.add_heading("Descargo metodológico", level=2)
    doc.add_paragraph(
        "Los índices y umbrales usados por OMNIPLAST son referenciales y citados en "
        "docs/METODOLOGIA.md. Los resultados generados con datos sintéticos (proyecto "
        "demo) NO representan mediciones reales y no deben citarse como evidencia "
        "ambiental. Valide los puntos de corte con la literatura vigente antes de uso "
        "regulatorio o de publicación científica.")
    EXPORT_DIR.mkdir(parents=True, exist_ok=True)
    filename = filename or f"informe_{datetime.now():%Y%m%d_%H%M}.docx"
    out = EXPORT_DIR / filename
    doc.save(out)
    return out


def informe_calidad(payload, proyecto="Proyecto activo", region="—"):
    """Informe ejecutivo de calidad del ecosistema a partir de un payload de indices."""
    secciones = []
    res = payload.get("resumen_general")
    if res:
        secciones.append({"titulo": "Resumen ejecutivo", "parrafos": [res]})
    for sitio, datos in payload.get("sitios", {}).items():
        sec = {"titulo": f"Estación {sitio}",
               "parrafos": [datos.get("interpretacion", "")],
               "tablas": {k: v for k, v in datos.items()
                          if isinstance(v, pd.DataFrame)}}
        secciones.append(sec)
    if payload.get("recomendaciones") is not None:
        secciones.append({"titulo": "Recomendaciones de gestión",
                          "tablas": {"prioridades": payload["recomendaciones"]}})
    refs = payload.get("referencias", [
        "Hakanson (1980) Hydrobiologia 91:832-848",
        "Lithner et al. (2011) ES&T 45:4727-4735",
        "Brown et al. (1970) WQI-NSF",
        "Carlson (1977) Limnol. Oceanogr. 22:361-369",
    ])
    return generate_docx("Informe integrado de calidad frente a microplásticos",
                         secciones,
                         meta={"proyecto": proyecto, "region": region,
                               "referencias": refs})
