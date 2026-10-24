"""Capa de acceso a datos: sesiones, dataframes planos y persistencia de resultados."""
from functools import lru_cache

import pandas as pd
from sqlalchemy import desc, select
from sqlalchemy.orm import Session

from ..config import DB_PATH
from .models import (AgentRun, Base, BioParam, Campaign, EnvParam, GeneticResult,
                     GeneticRun, LiteratureRef, MPRecord, Project, ResultStore,
                     Sample, Site, create_all, init_engine)


@lru_cache(maxsize=1)
def get_engine():
    engine = init_engine(DB_PATH)
    create_all(engine)
    return engine


def get_session():
    return Session(get_engine())


PARAM_META = {
    "temperatura_c": ("Temperatura", "°C", "fisicoquimico"),
    "ph": ("pH", "—", "fisicoquimico"),
    "od_mg_l": ("Oxígeno disuelto", "mg/L", "fisicoquimico"),
    "dbo_mg_l": ("DBO₅", "mg/L", "fisicoquimico"),
    "dqo_mg_l": ("DQO", "mg/L", "fisicoquimico"),
    "turbidez_ntu": ("Turbidez", "NTU", "fisicoquimico"),
    "salinidad_psu": ("Salinidad", "PSU", "fisicoquimico"),
    "conductividad_us": ("Conductividad", "µS/cm", "fisicoquimico"),
    "tss_mg_l": ("Sólidos suspendidos", "mg/L", "fisicoquimico"),
    "nitrato_mg_l": ("Nitrato", "mg/L", "fisicoquimico"),
    "nitrito_mg_l": ("Nitrito", "mg/L", "fisicoquimico"),
    "fosfato_mg_l": ("Fosfato", "mg/L", "fisicoquimico"),
    "amonio_mg_l": ("Amonio", "mg/L", "fisicoquimico"),
    "clorofila_a_ug_l": ("Clorofila a", "µg/L", "biologico"),
    "toc_mg_l": ("COT agua", "mg/L", "fisicoquimico"),
    "toc_pct": ("Materia orgánica sed.", "%", "sedimento"),
    "arena_pct": ("Arena", "%", "sedimento"),
    "limo_arcilla_pct": ("Limo + arcilla", "%", "sedimento"),
    "coliformes_totales_nmp": ("Coliformes totales", "NMP/100 mL", "microbiologico"),
    "e_coli_nmp": ("E. coli", "NMP/100 mL", "microbiologico"),
    "sod_u_mg_prot": ("SOD", "U/mg prot", "biomarcador"),
    "cat_u_mg_prot": ("CAT", "U/mg prot", "biomarcador"),
    "gst_u_mg_prot": ("GST", "U/mg prot", "biomarcador"),
    "lpo_nmol_tba_mg": ("LPO (TBARS)", "nmol/mg", "biomarcador"),
    "factor_condicion": ("Factor de condición", "—", "biologico"),
    "hsi_pct": ("Índice hepatosomático", "%", "biologico"),
    "gsi_pct": ("Índice gonadosomático", "%", "biologico"),
}


def _q(stmt):
    with get_session() as s:
        rows = s.execute(stmt).mappings().all()
        return pd.DataFrame([dict(r) for r in rows])


def sites_df():
    return _q(select(Site.id, Site.codigo, Site.nombre, Site.lat, Site.lon,
                     Site.zona, Site.descripcion))


def campaigns_df():
    return _q(select(Campaign.id, Campaign.codigo, Campaign.nombre, Campaign.temporada,
                     Campaign.fecha_inicio, Campaign.fecha_fin, Campaign.objetivo))


def samples_df():
    stmt = (select(Sample.id, Sample.codigo, Sample.fecha, Sample.matriz, Sample.subtipo,
                  Sample.profundidad_m, Sample.replicas, Sample.metodo, Sample.colector,
                  Sample.qa_qc,
                  Site.id.label("site_id"), Site.codigo.label("sitio"), Site.nombre.label("sitio_nombre"),
                  Site.lat, Site.lon, Site.zona,
                  Campaign.id.label("campaign_id"), Campaign.codigo.label("campana"),
                  Campaign.temporada)
            .join(Site, Sample.site_id == Site.id)
            .join(Campaign, Sample.campaign_id == Campaign.id))
    df = _q(stmt)
    if not df.empty:
        df["fecha"] = pd.to_datetime(df["fecha"])
        df["anio"] = df["fecha"].dt.year
    return df


def _samples_meta():
    smp = samples_df()
    if smp.empty:
        return smp
    return smp.rename(columns={"id": "sample_id"}).drop(
        columns=["site_id", "campaign_id"], errors="ignore")


def mp_df():
    stmt = (select(MPRecord.id, MPRecord.sample_id, MPRecord.polimero, MPRecord.morfologia,
                   MPRecord.color, MPRecord.tamano_min_um, MPRecord.tamano_max_um,
                   MPRecord.abundancia, MPRecord.unidad, MPRecord.metodo_id, MPRecord.peso_mg)
            .join(Sample, MPRecord.sample_id == Sample.id))
    df = _q(stmt)
    if df.empty:
        return df
    smp = _samples_meta()
    df = df.merge(smp, on="sample_id", how="left")
    df["total_sample"] = df.groupby("sample_id")["abundancia"].transform("sum")
    return df


def env_df():
    stmt = (select(EnvParam.sample_id, EnvParam.parametro, EnvParam.valor, EnvParam.unidad)
            .join(Sample, EnvParam.sample_id == Sample.id))
    df = _q(stmt)
    if df.empty:
        return df
    smp = _samples_meta()
    return df.merge(smp, on="sample_id", how="left")


def env_wide():
    df = env_df()
    if df.empty:
        return pd.DataFrame()
    wide = df.pivot_table(index="sample_id", columns="parametro", values="valor", aggfunc="mean")
    meta = df.drop_duplicates("sample_id")[[
        "sample_id", "codigo", "sitio", "sitio_nombre", "zona", "lat", "lon",
        "campana", "temporada", "fecha", "matriz", "subtipo", "profundidad_m", "anio"]].set_index("sample_id")
    return meta.join(wide).reset_index()


def bio_df():
    stmt = (select(BioParam.sample_id, BioParam.taxon, BioParam.tipo, BioParam.parametro,
                   BioParam.valor, BioParam.unidad, BioParam.tejido)
            .join(Sample, BioParam.sample_id == Sample.id))
    df = _q(stmt)
    if df.empty:
        return df
    smp = _samples_meta()
    return df.merge(smp, on="sample_id", how="left")


def qpcr_df():
    stmt = (select(GeneticResult.sample_id, GeneticResult.target, GeneticResult.metrica,
                  GeneticResult.valor, GeneticResult.ct, GeneticResult.detectado,
                  GeneticRun.tipo, GeneticRun.objetivo)
            .join(GeneticRun, GeneticResult.run_id == GeneticRun.id))
    df = _q(stmt)
    if df.empty:
        return df
    if df["tipo"].eq("qpcr").any():
        df = df[df["tipo"] == "qpcr"]
    smp = _samples_meta()
    return df.merge(smp, on="sample_id", how="left")


def microbiome_df():
    stmt = (select(GeneticResult.sample_id, GeneticResult.target, GeneticResult.valor)
            .join(GeneticRun, GeneticResult.run_id == GeneticRun.id)
            .where(GeneticRun.tipo == "16S"))
    df = _q(stmt)
    if df.empty:
        return pd.DataFrame()
    smp = _samples_meta()
    return df.merge(smp, on="sample_id", how="left")


def microbiome_wide():
    df = microbiome_df()
    if df.empty:
        return pd.DataFrame()
    wide = df.pivot_table(index="sample_id", columns="target", values="valor", aggfunc="sum")
    meta = df.drop_duplicates("sample_id")[[
        "sample_id", "codigo", "sitio", "zona", "campana", "temporada", "matriz", "anio"]].set_index("sample_id")
    return meta.join(wide).reset_index()


def save_result(pagina, analisis, resumen, payload):
    with get_session() as s:
        s.add(ResultStore(pagina=pagina, analisis=analisis, resumen=resumen, payload=payload))
        s.commit()


def list_results(limit=50):
    return _q(select(ResultStore.id, ResultStore.pagina, ResultStore.analisis,
                     ResultStore.resumen, ResultStore.creado).order_by(desc(ResultStore.id)).limit(limit))


def log_agent_run(agent_id, agent_name, status, resumen, payload, duracion):
    with get_session() as s:
        s.add(AgentRun(agent_id=agent_id, agent_name=agent_name, status=status,
                       resumen=resumen, payload=payload, duracion_s=duracion))
        s.commit()


def agent_runs_df(limit=100):
    return _q(select(AgentRun.id, AgentRun.agent_id, AgentRun.agent_name,
                     AgentRun.status, AgentRun.duracion_s, AgentRun.creado)
              .order_by(desc(AgentRun.id)).limit(limit))


def save_literature(refs):
    with get_session() as s:
        for r in refs:
            s.add(LiteratureRef(**r))
        s.commit()


def literature_df():
    return _q(select(LiteratureRef.id, LiteratureRef.titulo, LiteratureRef.autores,
                     LiteratureRef.anio, LiteratureRef.revista, LiteratureRef.doi,
                     LiteratureRef.fuente))


def project_summary():
    with get_session() as s:
        return {
            "proyectos": s.query(Project).count(),
            "sitios": s.query(Site).count(),
            "campanas": s.query(Campaign).count(),
            "muestras": s.query(Sample).count(),
            "mp_registros": s.query(MPRecord).count(),
            "env_registros": s.query(EnvParam).count(),
            "geneticos": s.query(GeneticResult).count(),
        }
