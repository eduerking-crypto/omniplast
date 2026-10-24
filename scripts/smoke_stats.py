import sys
sys.path.insert(0, ".")
import numpy as np
import pandas as pd
from omniplast.db import dao
from omniplast.core.statistics import (descriptive, inferential, multivariate,
                                       timeseries, ml, multilevel, bayes)
from omniplast.core import indices

mp = dao.mp_df()
water = mp[(mp.matriz == "agua")]
tot = water.groupby("sitio")["abundancia"].sum()
print("MP agua/sitio:", tot.round(0).to_dict())

env = dao.env_wide()
envw = env[env.matriz == "agua"]

r = inferential.kruskal_wallis(water.assign(t=water.abundancia), "t", "zona")
print("KW zona:", round(r["estadistico"], 3), r["p"])
r2 = inferential.anova_one(water, "abundancia", "temporada")
print("ANOVA temporada:", round(r2["estadistico"], 3), r2["p"], r2.get("eta2"))
d = inferential.dunn_posthoc(water, "abundancia", "zona")
print(d.head(3).to_string())

mb = dao.microbiome_wide()
if not mb.empty:
    taxa = mb.columns[7:]
    dist = multivariate.distance_matrix(mb[taxa], method="bray")
    perm = multivariate.permanova(dist, mb["zona"], permutations=99)
    print("PERMANOVA microbioma:", round(perm["pseudo_F"], 3), perm["p"], round(perm["R2"], 3))
    nmds = multivariate.nmds(dist, n_init=4)
    print("NMDS stress:", round(nmds["stress"], 4), nmds["calidad_ajuste"])
    sim = multivariate.simper(mb[taxa], mb["zona"])
    first_par = list(sim["por_par"].keys())[0]
    print("SIMPER top:", sim["por_par"][first_par].head(3)[["variable", "contribucion_pct"]].to_string())

ph = dao.qpcr_df()
merged = ph[ph.target == "petase"].merge(
    mp.groupby("sample_id")["abundancia"].sum().rename("mp_total"),
    left_on="sample_id", right_index=True)
glm = multilevel.glm_counts(merged.assign(mp_total=merged.mp_total), "valor", ["mp_total"], family="negativebinomial")
print("GLM NB petase~mp:", glm.get("interpretacion"), "alpha:", glm.get("alpha_sobredispersion"))

sites_mp = mp[mp.matriz == "sedimento"].groupby("sitio")["abundancia"].sum()
bg = sites_mp.min()
plis = {s: indices.pli_site([indices.contamination_factor(v, bg) for _ in [0]]) for s, v in sites_mp.items()}
print("PLI JPC-01:", plis.get("JPC-01"))
phi = indices.polymer_hazard_index(mp[mp.sitio == "JPC-01"])
print("PHI JPC-01:", round(phi["phi"], 1), phi["clase"], "unk:", round(phi["desconocidos_pct"], 1))

comb = dao.env_wide()
comb_s = comb[comb.matriz == "sedimento"].merge(
    mp[mp.matriz == "sedimento"].groupby("sample_id")["abundancia"].sum().rename("mp_total"),
    on="sample_id", how="left").merge(
    ph[ph.target == "petase"].pivot_table(index="sample_id", values="valor", aggfunc="max"),
    on="sample_id", how="left").rename(columns={"valor": "petase"})
comb_s = comb_s.dropna(subset=["mp_total", "petase", "toc_pct", "arena_pct"])
rda = multivariate.rda(comb_s[["toc_pct", "arena_pct", "salinidad_psu", "nitrato_mg_l"]],
                       comb_s[["mp_total", "petase"]], permutations=99)
print("RDA:", rda["interpretacion"])

comb_w2 = comb[comb.matriz == "agua"].dropna(subset=["dbo_mg_l", "od_mg_l"])
sample_dict = comb_w2.iloc[0].to_dict()
wqi = indices.wqi_nsf(sample_dict)
print("WQI:", wqi["wqi"], wqi["clase"])
site_row = {"mp_total": sites_mp.max()}
iicam = indices.iicam(site_row, sites_mp, phi_site=phi["phi"], wqi_site=wqi["wqi"])
print("IICAM JPC-01:", iicam["iicam"], iicam["clase"], iicam["ic95_mc"])

mk = timeseries.mann_kendall(pd.Series([10, 12, 11, 14, 15, 17, 18, 20]))
print("MK:", mk["Z"], mk["interpretacion"])
b = bayes.mcmc_two_means(np.random.normal(5, 1, 30), np.random.normal(6, 1, 30), n_iter=3000)
print("MCMC:", b["interpretacion"])
rf = ml.random_forest_regressor(comb_s,
    ["toc_pct", "arena_pct", "salinidad_psu", "nitrato_mg_l", "fosfato_mg_l"], "mp_total")
print("RF:", rf.get("interpretacion", rf.get("error")))
print("OK: motor estadistico completo verificado")
