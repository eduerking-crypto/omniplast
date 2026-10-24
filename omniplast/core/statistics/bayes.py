"""Estadistica bayesiana: conjugados, MCMC simple, factores de Bayes por BIC."""
import numpy as np
from scipy import stats


def normal_posterior(x, prior_mu=0.0, prior_sd=10.0, known_sd=None):
    """Posterior normal-normal para la media con varianza conocida (o estimada)."""
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    n = len(x)
    if n < 2:
        return {"test": "Posterior normal", "error": "n≥2"}
    sd = known_sd if known_sd is not None else np.std(x, ddof=1)
    post_prec = 1 / prior_sd ** 2 + n / sd ** 2
    post_var = 1 / post_prec
    post_mu = post_var * (prior_mu / prior_sd ** 2 + x.sum() / sd ** 2)
    ci = (post_mu - 1.96 * np.sqrt(post_var), post_mu + 1.96 * np.sqrt(post_var))
    return {"test": "Conjugado normal-normal", "media_prior": prior_mu, "media_posterior": float(post_mu),
            "sd_posterior": float(np.sqrt(post_var)),
            "ic95_credible": (float(ci[0]), float(ci[1])),
            "p_mayor_que_cero": float(1 - stats.norm.cdf(0, post_mu, np.sqrt(post_var))),
            "interpretacion": f"Media posterior={post_mu:.4g} (SD={np.sqrt(post_var):.4g}); "
                              f"probabilidad de que μ>0: {1 - stats.norm.cdf(0, post_mu, np.sqrt(post_var)):.1%}"}


def beta_binomial(k, n, a=1, b=1):
    """Posterior Beta para una proporcion (deteccion de un polimero, etc.)."""
    if n <= 0:
        return {"test": "Beta-binomial", "error": "n>0"}
    a_post, b_post = a + k, b + (n - k)
    mean = a_post / (a_post + b_post)
    lo, hi = stats.beta.ppf([.025, .975], a_post, b_post)
    return {"test": "Beta-binomial conjugado", "media": float(mean),
            "ic95_credible": (float(lo), float(hi)),
            "a_post": float(a_post), "b_post": float(b_post),
            "interpretacion": f"Proporción estimada={mean:.3f}, IC creíble 95%=[{lo:.3f}, {hi:.3f}]"}


def _log_posterior_diff(m, a, b, sd_a, sd_b, prior_mu=0, prior_sd=10):
    if len(a) == 0 or len(b) == 0:
        return -np.inf
    ll = ((-0.5 * ((a - m[0]) / sd_a) ** 2).sum() + (-0.5 * ((b - m[1]) / sd_b) ** 2).sum())
    lp = -0.5 * ((m - prior_mu) / prior_sd) ** 2 @ np.eye(2) @ (m - prior_mu) / 2
    return ll + lp


def mcmc_two_means(a, b, n_iter=12000, burn=2000, seed=42):
    """Metropolis-Hastolds para la diferencia de medias (estilo BEST simplificado)."""
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) < 3 or len(b) < 3:
        return {"test": "MCMC", "error": "n≥3 por grupo"}
    sd_a, sd_b = np.std(a, ddof=1), np.std(b, ddof=1)
    current = np.array([np.mean(a), np.mean(b)])
    chain = []
    for _ in range(n_iter):
        prop = current + rng.normal(0, [sd_a * .2, sd_b * .2])
        lp_cur = _log_posterior_diff(current, a, b, sd_a, sd_b)
        lp_prop = _log_posterior_diff(prop, a, b, sd_a, sd_b)
        if np.log(rng.random()) < lp_prop - lp_cur:
            current = prop
        chain.append(current.copy())
    chain = np.array(chain[burn:])
    diff = chain[:, 0] - chain[:, 1]
    lo, hi = np.percentile(diff, [2.5, 97.5])
    p_sup = float((diff > 0).mean())
    return {"test": "Metropolis-Hastings (2 medias)", "diferencia_media": float(diff.mean()),
            "ic95_credible": (float(lo), float(hi)), "p_dif_mayor_0": p_sup,
            "n_muestras": int(len(diff)),
            "interpretacion": f"Δ={diff.mean():.4g} (IC creíble 95% [{lo:.4g}, {hi:.4g}]); "
                              f"P(Δ>0)={p_sup:.1%} → {'diferencia creíble' if p_sup < .05 or p_sup > .95 else 'diferencia no concluyente'}"}


def bayes_factor_bic(bic1, bic2):
    """Aproximacion del factor de Bayes a partir de BIC (Kass & Raftery 1995)."""
    bf = np.exp((bic2 - bic1) / 2)
    if bf < 1:
        strength = "evidencia débil o nula"
    elif bf < 3:
        strength = "evidencia positiva"
    elif bf < 20:
        strength = "evidencia fuerte"
    elif bf < 150:
        strength = "evidencia muy fuerte"
    else:
        strength = "evidencia decisiva"
    return {"test": "Factor de Bayes (BIC)", "BF": float(bf),
            "interpretacion": f"BF≈{bf:.2g} ({strength}) a favor del modelo con menor BIC"}


def power_analysis(effect_d=0.5, alpha=0.05, power=0.8, n=None, groups=2):
    """Tamano de muestra o potencia alcanzada para t-test / ANOVA aproximado."""
    from scipy import stats as st
    from scipy.optimize import brentq
    if n is None:
        lam = (groups + (groups > 2)) * 0
        if groups == 2:
            z_a = st.norm.ppf(1 - alpha / 2)
            z_b = st.norm.ppf(power)
            n = 2 * ((z_a + z_b) / effect_d) ** 2
            return {"test": "Análisis de potencia (aprox. normal)", "n_por_grupo": int(np.ceil(n)),
                    "n_total": int(np.ceil(n) * 2), "d_esperado": effect_d,
                    "interpretacion": f"Se requieren ≈{int(np.ceil(n))} observaciones por grupo "
                                      f"para detectar d={effect_d} con potencia {power:.0%}"}
    else:
        z_a = st.norm.ppf(1 - alpha / 2)
        se = np.sqrt(2 / n)
        power_calc = st.norm.cdf(effect_d / se - z_a) + st.norm.cdf(-effect_d / se - z_a)
        return {"test": "Análisis de potencia", "potencia_estimada": float(power_calc),
                "interpretacion": f"Con n={n} por grupo y d={effect_d}, potencia≈{power_calc:.0%} "
                                  f"({'suficiente' if power_calc >= .8 else 'insuficiente (riesgo de error tipo II)'})"}
