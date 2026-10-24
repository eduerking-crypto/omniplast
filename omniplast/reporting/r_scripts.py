"""Generador de scripts R ejecutables en RStudio (vegan, DESeq2, lavaan, glmmTMB).

Los scripts se escriven en data/exports/*.R con lectura de CSV parametrizada para
que el usuario ejecute el analisis canonico en R cuando el flujo doctoral lo exija.
"""
from pathlib import Path

from ..config import EXPORT_DIR, R_AVAILABLE_HINT


def _write(name, code):
    EXPORT_DIR.mkdir(parents=True, exist_ok=True)
    path = EXPORT_DIR / name
    path.write_text(code, encoding="utf-8")
    return path


VEGAN_TEMPLATE = '''# =============================================================
# OMNIPLAST :: Analisis multivariante canónico (vegan) en RStudio
# Requiere: install.packages(c("vegan", "ggplot2", "readr"))
# Uso: Rscript omp_vegan.R < abundancias.csv> <grupos.csv> <columna_grupo>
# =============================================================
library(vegan); library(ggplot2); library(readr)
args <- commandArgs(trailingOnly = TRUE)
abund  <- read.csv(args[1], row.names = 1)
grupos <- read.csv(args[2])
gcol   <- args[3]

mat <- as.matrix(abund)
mat_hel <- decostand(mat, "hellinger")          # transformación Hellinger
dist <- vegdist(mat, method = "bray")           # Bray-Curtis (raíz ya aplicada por hellinger)

# NMDS con 250 intentos (estándar metaMDS)
nmds <- metaMDS(mat, distance = "bray", trymax = 250)
print(paste("Stress:", nmds$stress))

# PERMANOVA + betadisper (SIEMPRE reportar ambos)
perm <- adonis2(mat ~ grupos[[gcol]], data = grupos, method = "bray", permutations = 999)
print(perm)
disp <- betadisper(vegdist(mat, "bray"), grupos[[gcol]])
print(permutest(disp))

# ANOSIM y SIMPER
print(anosim(dist, grupos[[gcol]]))
print(summary(simper(mat, grupos[[gcol]])))

# CAP (constrained Analysis of Principal Coordinates) opcional con variables ambientales
# cap <- capscale(mat ~ ., data = env, distance = "bray")

g <- data.frame(NMDS1 = scores(nmds)[, 1], NMDS2 = scores(nmds)[, 2],
                grupo = grupos[[gcol]])
ggplot(g, aes(NMDS1, NMDS2, color = grupo)) +
  geom_point(size = 3) + stat_ellipse(level = 0.95) +
  theme_minimal() + ggtitle("NMDS Bray-Curtis")
ggsave("omp_nmds.png", width = 7, height = 5)
'''

DESEQ2_TEMPLATE = '''# =============================================================
# OMNIPLAST :: Abundancia diferencial canónica (DESeq2) en RStudio
# Requiere: BiocManager::install("DESeq2")
# Uso: Rscript omp_deseq2.R <conteos.csv> <grupos.csv> <columna_grupo>
# =============================================================
library(DESeq2)
args <- commandArgs(trailingOnly = TRUE)
counts <- read.csv(args[1], row.names = 1)
meta   <- read.csv(args[2])
gcol   <- args[3]

dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = meta,
                              design = as.formula(paste("~", gcol)))
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
resOrdered <- res[order(res$padj), ]
print(head(resOrdered, 25))
write.csv(as.data.frame(resOrdered), "omp_deseq2_resultados.csv")
'''

LAVAAN_TEMPLATE = '''# =============================================================
# OMNIPLAST :: Modelos de ecuaciones estructurales (SEM) en RStudio
# Requiere: install.packages(c("lavaan", "semPlot"))
# El modelo ejemplo: DBO -> MP -> (LPO) con biomarcador como mediador.
# Ajuste las rutas según su teoría causal antes de interpretar.
# =============================================================
library(lavaan)
datos <- read.csv("datos_sem.csv")

modelo <- '
  # rutas estructurales
  mp     ~ dbo + tss + salinidad
  lpo    ~ mp + dbo
  # covarianzas
  dbo ~~ tss
'
fit <- sem(modelo, data = datos, estimator = "MLR")
print(summary(fit, standardized = TRUE, fit.measures = TRUE, ci = TRUE))
# Criterios doctorales: CFI > 0.95, RMSEA < 0.06, SRMR < 0.08 (Hu & Bentler 1999)
'''

GLMMTMB_TEMPLATE = '''# =============================================================
# OMNIPLAST :: GLMM de conteos de MP (glmmTMB) en RStudio
# Requiere: install.packages(c("glmmTMB", "DHARMa"))
# Binomial negativa con intercepto aleatorio por sitio.
# =============================================================
library(glmmTMB); library(DHARMa)
datos <- read.csv("datos_glmm.csv")

m <- glmmTMB(abundancia ~ dbo + tss + temporada + (1 | sitio),
             family = nbinom2, data = datos)
print(summary(m))

# Diagnóstico de residuos simulados (obligatorio para publicación)
sim <- simulateResiduals(m)
plot(sim)
testDispersion(sim)
'''


def generate_all(csv_abundancia=None, csv_grupos=None, grupo_col="zona"):
    ab = csv_abundancia or "abundancias.csv"
    gr = csv_grupos or "grupos.csv"
    paths = {
        "omp_vegan.R": VEGAN_TEMPLATE.replace("args[1]", f"'{ab}'").replace("args[2]", f"'{gr}'").replace("args[3]", f"'{grupo_col}'"),
        "omp_deseq2.R": DESEQ2_TEMPLATE.replace("args[1]", f"'{ab}'").replace("args[2]", f"'{gr}'").replace("args[3]", f"'{grupo_col}'"),
        "omp_sem_lavaan.R": LAVAAN_TEMPLATE,
        "omp_glmm_glmmTMB.R": GLMMTMB_TEMPLATE,
    }
    out = {name: _write(name, code) for name, code in paths.items()}
    return {"scripts": out, "nota": R_AVAILABLE_HINT}


def run_rscript(path, timeout=90):
    """Ejecuta Rscript si R está en PATH; devuelve stdout o None."""
    import shutil
    import subprocess
    if shutil.which("Rscript") is None:
        return {"ok": False, "error": "Rscript no está en PATH. " + R_AVAILABLE_HINT}
    try:
        res = subprocess.run(["Rscript", str(path)], capture_output=True,
                             text=True, timeout=timeout)
        return {"ok": res.returncode == 0, "salida": res.stdout[-4000:],
                "error": res.stderr[-2000:] if res.returncode else None}
    except Exception as e:
        return {"ok": False, "error": str(e)}
