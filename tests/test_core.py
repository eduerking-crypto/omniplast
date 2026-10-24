"""Suite de pruebas OMNIPLAST: nucleo cientifico, indices, genetica, geo y agentes."""
import numpy as np
import pandas as pd


def _core():
    from omniplast.core import indices, sources, uncertainty
    return indices, sources, uncertainty


def _stats():
    from omniplast.core.statistics import inferential, timeseries, multivariate
    return inferential, timeseries, multivariate


def test_indices_basicos():
    indices, _, _ = _core()
    assert indices.contamination_factor(200, 50) == 4.0
    pli = indices.pli_site([4.0, 1.0, 1.0, 1.0])
    assert pli > 1.0
    assert indices.cf_classify(4.0) == "Considerable"
    assert indices.cf_classify(0.5) == "Baja"


def test_phi():
    indices, _, _ = _core()
    df = pd.DataFrame({"polimero": ["PE", "PVC", "UNKNOWN"],
                       "abundancia": [50, 25, 25]})
    r = indices.polymer_hazard_index(df)
    esperado = (11 * 50 + 5000 * 25) / 100
    assert abs(r["phi"] - esperado) < 0.5
    assert r["desconocidos_pct"] == 25.0
    assert r["clase"] == "Extremadamente peligroso"


def test_wqi_nsf():
    indices, _, _ = _core()
    r = indices.wqi_nsf({"od_mg_l": 8.0, "ph": 7.2, "dbo_mg_l": 2.0,
                         "temperatura_c": 25.0, "turbidez_ntu": 5.0,
                         "tss_mg_l": 10.0, "nitrato_mg_l": 0.5,
                         "nitrito_mg_l": 0.01, "fosfato_mg_l": 0.02,
                         "coliformes_totales_nmp": 10})
    assert 70 <= r["wqi"] <= 95
    assert r["clase"] in ("Excelente", "Buena")
    assert indices.wqi_nsf({})["wqi"] != indices.wqi_nsf({})["wqi"]


def test_tsi_carlson():
    indices, _, _ = _core()
    r = indices.tsi_carlson(chla=80)
    assert 60 < r["tsi"] < 80
    assert "Eutró" in r["clase"] or "eutró" in r["clase"]


def test_mp_tier():
    indices, _, _ = _core()
    assert indices.mp_tier(50, "agua") == "Baja"
    assert indices.mp_tier(10000, "agua") == "Extremadamente alta"
    assert indices.mp_tier(100, "sedimento") == "Baja"


def test_iicam_extremos():
    indices, _, _ = _core()
    serie = pd.Series([100, 500, 1000, 2000, 4000])
    mejor = indices.iicam({"mp_total": 100}, serie, phi_site=1, wqi_site=95)
    peor = indices.iicam({"mp_total": 4000}, serie, phi_site=5000, wqi_site=25)
    assert mejor["iicam"] > 80
    assert peor["iicam"] < 40
    assert mejor["clase"] == "Excelente"
    assert peor["clase"] == "Crítica"


def test_mann_kendall_y_sen():
    _, timeseries, _ = _stats()
    creciente = pd.Series(np.arange(10) * 2.0 + np.random.default_rng(1).normal(0, .1, 10))
    r = timeseries.mann_kendall(creciente)
    assert r["Z"] > 1.96 and r["tendencia"] == "creciente"
    s = timeseries.sens_slope(pd.Series([10, 12, 11, 14, 15, 17]))
    assert abs(s["pendiente"] - 1.4) < 0.5


def test_dunn_y_kruskal():
    inferential, _, _ = _stats()
    rng = np.random.default_rng(2)
    df = pd.DataFrame({"valor": np.concatenate([rng.normal(10, 1, 30),
                                                rng.normal(13, 1, 30),
                                                rng.normal(16, 1, 30)]),
                       "grupo": ["A"] * 30 + ["B"] * 30 + ["C"] * 30})
    r = inferential.kruskal_wallis(df, "valor", "grupo")
    assert r["p"] < 0.001
    dunn = inferential.dunn_posthoc(df, "valor", "grupo")
    assert (dunn["p_holm"] < 0.001).all()
    d = inferential.cohen_d(df[df["grupo"] == "A"]["valor"],
                            df[df["grupo"] == "C"]["valor"])
    assert d < -2


def test_permanova_separa_grupos():
    _, _, multivariate = _stats()
    rng = np.random.default_rng(3)
    a = np.concatenate([rng.normal(0, 1, (10, 5))])
    b = np.concatenate([rng.normal(4, 1, (10, 5))])
    X = np.vstack([a, b])
    df = pd.DataFrame(X)
    dist = multivariate.distance_matrix(df, method="euclidean")
    grupos = pd.Series(["g1"] * 10 + ["g2"] * 10)
    r = multivariate.permanova(dist, grupos, permutations=199)
    assert r["pseudo_F"] > 10 and r["p"] < 0.05
    nmds = multivariate.nmds(dist, n_init=3)
    assert nmds["stress"] < 0.15


def test_qpcr_ddct_y_eficiencia():
    from omniplast.core.genetics import qpcr
    eff = qpcr.efficiency_from_slope(-3.32)
    assert 99 <= eff["eficiencia_pct"] <= 101
    dd = pd.DataFrame({
        "sample_id": [1, 2, 3, 4, 5, 6],
        "codigo": ["m1", "m2", "m3", "m4", "m5", "m6"],
        "zona": ["ref", "ref", "ref", "imp", "imp", "imp"],
        "target": ["petase"] * 6, "ct": [30, 30.2, 29.8, 26, 26.1, 25.9]})
    dd2 = dd.copy()
    dd2["target"] = "16s"
    dd2["ct"] = [16, 16.1, 15.9, 16.0, 16.2, 15.8]
    r = qpcr.ddct_table(pd.concat([dd, dd2]), target="petase", ref="16s",
                        control_group="ref")
    assert r["tabla"]["fc"].iloc[3] > 5  # sobreexpresado en impacto (~8x)


def test_diversidad_alfa():
    from omniplast.core.genetics import community
    mb = pd.DataFrame({"codigo": ["a", "b"], "zona": ["z1", "z2"],
                       "Pseudomonas": [100, 50], "Vibrio": [100, 150],
                       "Bacillus": [0, 100]})
    alpha = community.alpha_diversity(mb)
    assert alpha.loc[0, "riqueza"] == 2
    assert abs(alpha.loc[0, "shannon"] - np.log(2)) < 1e-6  # H' natural = ln 2
    assert alpha.loc[1, "riqueza"] == 3


def pytest_approx(x):
    class A:
        def __eq__(self, other):
            return abs(other - x) < 1e-6
    return A()


def test_hwe_y_fst():
    from omniplast.core.genetics import genotyping
    r = genotyping.hardy_weinberg_test([50, 100, 50])  # p=q=0.5 exactamente esperado
    assert r["p"] > 0.05
    fst = genotyping.fst_pairwise([0, 0, 1, 1], [0, 1, 1, 1])
    assert -0.2 < fst["fst"] < 0.4


def test_stokes_y_apportionment():
    from omniplast.core.geo import trajectory
    indices, sources, _ = _core()
    v_pe = trajectory.stokes_settling_velocity(0.5, 960)
    v_pvc = trajectory.stokes_settling_velocity(0.5, 1400)
    assert v_pe == 0.0 and v_pvc > 0
    mp = pd.DataFrame({"polimero": ["PET", "PET", "PA", "PA", "PE", "PP"],
                       "morfologia": ["fibra"] * 4 + ["fragmento"] * 2,
                       "abundancia": [10, 10, 5, 5, 10, 10]})
    r = sources.apportion_sources(mp)
    assert r["tabla"].iloc[0]["fuente"] == "Textiles y lavado doméstico"


def test_monte_carlo_y_hq():
    from omniplast.core.eco import toxicology
    indices, sources, uncertainty = _core()
    mc = uncertainty.monte_carlo(lambda a, b: a * b,
                                 {"a": ("normal", 10, 1), "b": ("normal", 5, 1)}, n=4000)
    assert abs(mc["salida"]["p50"] - 50) < 5
    r = toxicology.hazard_quotient_table({"PE": 2.0}, {"PE": 2.0})
    assert r["tabla"].iloc[0]["hq"] == 1.0
    assert r["tabla"].iloc[0]["clase"] == "Riesgo medio"


def test_plantillas_excel():
    from omniplast.reporting import plantillas
    vacia = plantillas.generar_plantilla_bytes()
    llena = plantillas.generar_ejemplo_bytes()
    assert len(vacia) > 10_000
    assert len(llena) > 80_000
    from openpyxl import load_workbook
    import io
    wb = load_workbook(io.BytesIO(vacia))
    esperadas = {"LEEME", "Sitios", "Campanas", "Muestras", "Microplasticos",
                 "Parametros_Ambientales", "Biota", "qPCR", "Microbioma_16S", "Catalogos"}
    assert esperadas <= set(wb.sheetnames)
    wb2 = load_workbook(io.BytesIO(llena))
    assert wb2["Microplasticos"].max_row > 1000


def test_registro_100_agentes():
    from omniplast.agents import REGISTRY, run_agent
    from omniplast.agents.engine import HANDLERS
    assert len(REGISTRY) == 100
    faltantes = set(REGISTRY) - set(HANDLERS)
    assert not faltantes
    for aid in ["AG-007", "AG-019", "AG-063", "AG-100"]:
        r = run_agent(aid)
        assert r["estado"] in ("ok", "requiere")
        assert r["resumen"]


def test_db_y_dao():
    from omniplast.db import dao
    s = dao.project_summary()
    assert s["sitios"] == 10 and s["campanas"] == 4
    assert dao.mp_df()["matriz"].nunique() == 3
    assert len(dao.microbiome_wide()) == 40
