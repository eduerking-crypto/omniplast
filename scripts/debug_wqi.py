import sys
sys.path.insert(0, ".")
from omniplast.db import dao
from omniplast.core.indices import _wqi_rating, wqi_nsf

env = dao.env_wide()
row = env[env.matriz == "agua"].dropna(subset=["dbo_mg_l"]).iloc[0].to_dict()
for p in ["od_mg_l", "coliformes_totales_nmp", "ph", "dbo_mg_l", "temperatura_c",
          "fosfato_mg_l", "nitrito_mg_l", "turbidez_ntu", "tss_mg_l"]:
    print(p, "valor:", row.get(p), "q:", _wqi_rating(p, row.get(p)))
res = wqi_nsf(row)
print(res["wqi"], res["clase"], res["ratings"])
