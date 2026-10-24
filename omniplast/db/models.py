"""Modelos ORM (SQLite) para proyectos, sitios, campanas, muestras y resultados."""
from datetime import datetime

from sqlalchemy import (Boolean, Column, DateTime, Float, ForeignKey, Integer,
                        JSON, String, Text, create_engine)
from sqlalchemy.orm import declarative_base, relationship

Base = declarative_base()


class Project(Base):
    __tablename__ = "projects"
    id = Column(Integer, primary_key=True)
    codigo = Column(String(40), unique=True)
    nombre = Column(String(200))
    descripcion = Column(Text, default="")
    ecosistema = Column(String(60))
    region = Column(String(200))
    responsable = Column(String(200))
    fecha_inicio = Column(DateTime)
    fecha_fin = Column(DateTime)
    creado = Column(DateTime, default=datetime.utcnow)


class Site(Base):
    __tablename__ = "sites"
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey("projects.id"))
    codigo = Column(String(30))
    nombre = Column(String(200))
    lat = Column(Float)
    lon = Column(Float)
    zona = Column(String(40))
    descripcion = Column(Text, default="")
    sites = relationship("Project")


class Campaign(Base):
    __tablename__ = "campaigns"
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey("projects.id"))
    codigo = Column(String(30))
    nombre = Column(String(200))
    temporada = Column(String(20))
    fecha_inicio = Column(DateTime)
    fecha_fin = Column(DateTime)
    objetivo = Column(Text, default="")


class Sample(Base):
    __tablename__ = "samples"
    id = Column(Integer, primary_key=True)
    campaign_id = Column(Integer, ForeignKey("campaigns.id"))
    site_id = Column(Integer, ForeignKey("sites.id"))
    codigo = Column(String(50))
    fecha = Column(DateTime)
    matriz = Column(String(20))
    subtipo = Column(String(60))
    profundidad_m = Column(Float, default=0.0)
    replicas = Column(Integer, default=3)
    metodo = Column(String(120))
    colector = Column(String(120))
    qa_qc = Column(JSON, default=dict)
    notas = Column(Text, default="")


class MPRecord(Base):
    __tablename__ = "mp_records"
    id = Column(Integer, primary_key=True)
    sample_id = Column(Integer, ForeignKey("samples.id"))
    polimero = Column(String(20))
    morfologia = Column(String(20))
    color = Column(String(20))
    tamano_min_um = Column(Float)
    tamano_max_um = Column(Float)
    abundancia = Column(Float)
    unidad = Column(String(20))
    peso_mg = Column(Float, default=0.0)
    metodo_id = Column(String(60))
    notas = Column(Text, default="")


class EnvParam(Base):
    __tablename__ = "env_params"
    id = Column(Integer, primary_key=True)
    sample_id = Column(Integer, ForeignKey("samples.id"))
    parametro = Column(String(40))
    valor = Column(Float)
    unidad = Column(String(30))
    metodo = Column(String(120), default="")
    detectado = Column(Boolean, default=True)


class BioParam(Base):
    __tablename__ = "bio_params"
    id = Column(Integer, primary_key=True)
    sample_id = Column(Integer, ForeignKey("samples.id"))
    taxon = Column(String(120))
    tipo = Column(String(40))
    parametro = Column(String(40))
    valor = Column(Float)
    unidad = Column(String(40))
    tejido = Column(String(60), default="")
    notas = Column(Text, default="")


class GeneticRun(Base):
    __tablename__ = "genetic_runs"
    id = Column(Integer, primary_key=True)
    sample_id = Column(Integer, ForeignKey("samples.id"))
    tipo = Column(String(20))
    objetivo = Column(String(200))
    plataforma = Column(String(80), default="")
    protocolo = Column(Text, default="")


class GeneticResult(Base):
    __tablename__ = "genetic_results"
    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, ForeignKey("genetic_runs.id"))
    sample_id = Column(Integer, ForeignKey("samples.id"))
    target = Column(String(120))
    metrica = Column(String(30))
    valor = Column(Float)
    ct = Column(Float)
    unidad = Column(String(40))
    detectado = Column(Boolean, default=True)


class ResultStore(Base):
    __tablename__ = "results_store"
    id = Column(Integer, primary_key=True)
    pagina = Column(String(80))
    analisis = Column(String(200))
    resumen = Column(Text, default="")
    payload = Column(JSON)
    creado = Column(DateTime, default=datetime.utcnow)


class AgentRun(Base):
    __tablename__ = "agent_runs"
    id = Column(Integer, primary_key=True)
    agent_id = Column(String(20))
    agent_name = Column(String(200))
    status = Column(String(20))
    resumen = Column(Text, default="")
    payload = Column(JSON)
    duracion_s = Column(Float, default=0.0)
    creado = Column(DateTime, default=datetime.utcnow)


class LiteratureRef(Base):
    __tablename__ = "literature_refs"
    id = Column(Integer, primary_key=True)
    titulo = Column(Text)
    autores = Column(Text, default="")
    anio = Column(Integer)
    revista = Column(String(300), default="")
    doi = Column(String(120), default="")
    fuente = Column(String(60), default="")
    resumen = Column(Text, default="")
    etiquetas = Column(String(300), default="")
    creado = Column(DateTime, default=datetime.utcnow)


def init_engine(db_path):
    return create_engine(f"sqlite:///{db_path}", connect_args={"check_same_thread": False})


def create_all(engine):
    Base.metadata.create_all(engine)
