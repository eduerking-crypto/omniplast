"""Motor de incertidumbre: Monte Carlo generico y analisis de sensibilidad."""
import numpy as np
import pandas as pd
from scipy import stats


def monte_carlo(func, param_specs, n=20000, seed=42, percentiles=(2.5, 25, 50, 75, 97.5)):
    """param_specs: {nombre: ("normal", mu, sd) | ("uniform", a, b) | ("triangular", a, c, b)
    | ("lognormal", mu, sigma) | ("beta", a, b) | ("poisson", lam)}."""
    rng = np.random.default_rng(seed)
    draws = {}
    for name, spec in param_specs.items():
        kind = spec[0].lower()
        if kind == "normal":
            draws[name] = rng.normal(spec[1], spec[2], n)
        elif kind == "uniform":
            draws[name] = rng.uniform(spec[1], spec[2], n)
        elif kind == "triangular":
            draws[name] = rng.triangular(spec[1], spec[2], spec[3], n)
        elif kind == "lognormal":
            draws[name] = rng.lognormal(spec[1], spec[2], n)
        elif kind == "beta":
            draws[name] = rng.beta(spec[1], spec[2], n)
        elif kind == "poisson":
            draws[name] = rng.poisson(spec[1], n).astype(float)
        else:
            raise ValueError(f"Distribución {kind} no soportada")
    out = np.array([func(**{k: v[i] for k, v in draws.items()}) for i in range(n)])
    result = {f"p{p}": float(np.percentile(out, p)) for p in percentiles}
    result["media"] = float(np.mean(out))
    result["sd"] = float(np.std(out))
    sens = {}
    for k, v in draws.items():
        try:
            rho, p = stats.spearmanr(v, out)
            sens[k] = {"rho": float(rho), "p": float(p)}
        except Exception:
            pass
    sens_df = pd.DataFrame(sens).T.sort_values("rho", key=abs, ascending=False)
    return {"test": f"Monte Carlo (n={n})", "salida": result,
            "sensibilidad": sens_df, "muestra": draws,
            "interpretacion": f"Mediana={np.median(out):.4g}, banda P2.5-P97.5=["
                              f"{np.percentile(out, 2.5):.4g}, {np.percentile(out, 97.5):.4g}]"}


def tornado_sensitivity(base_func, param_specs, n=5000, delta=0.10, seed=42):
    """Analisis tornado: efecto de ±delta en cada parametro sobre la salida media."""
    rng = np.random.default_rng(seed)
    draws = {name: rng.normal(spec[1], spec[2], n) if spec[0] == "normal" else
             rng.uniform(spec[1], spec[2], n) for name, spec in param_specs.items()}
    base = np.mean([base_func(**{k: v[i] for k, v in draws.items()}) for i in range(n)])
    rows = []
    for name, spec in param_specs.items():
        up = dict(draws)
        up[name] = np.clip(draws[name] * (1 + delta), 1e-12, None)
        dn = dict(draws)
        dn[name] = np.clip(draws[name] * (1 - delta), 1e-12, None)
        val_up = np.mean([base_func(**{k: v[i] for k, v in up.items()}) for i in range(n)])
        val_dn = np.mean([base_func(**{k: v[i] for k, v in dn.items()}) for i in range(n)])
        rows.append({"parametro": name, "menor": val_dn, "base": base, "mayor": val_up,
                     "rango": abs(val_up - val_dn)})
    tab = pd.DataFrame(rows).sort_values("rango", ascending=False)
    top = tab.iloc[0]
    return {"test": f"Análisis tornado (±{delta:.0%})", "tabla": tab.round(4),
            "interpretacion": f"Parámetro más influyente: {top['parametro']} "
                              f"(rango {top['rango']:.4g} en la salida)"}
