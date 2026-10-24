# Metodología y fórmulas — OMNIPLAST 1.0

Cada índice/prueba implementado, su fórmula, sus clases de interpretación y su cita
primaria. Los umbrales marcados **[referencial]** no tienen valor normativo oficial para
microplásticos: son editables en `omniplast/config.py`.

## 1. Definiciones y protocolos

- **Microplástico**: 1 µm – 5 mm; **nanoplástico**: 1–1000 nm (Hartmann et al. 2019, ES&T 53:1039-1047).
- **Muestreo agua**: red neustón/manta 300–330 µm (GESAMP 2019); registrar volumen filtrado; unidades part/m³.
- **Muestreo sedimento**: núcleo/cuadrante, capa 0–5 cm, separación por densidad, unidades part/kg peso seco (Masura et al. 2015, NOAA).
- **Muestreo biota**: tracto digestivo aislado, n≥10, reportar % con MP y part/individuo (Lusher et al. 2013).
- **QA/QC**: blancos <5%, duplicados ≥10%, recuperación 70–130%, identificación espectral ≥10% (Koelmans et al. 2019, ES&T 53:9615-9623).
- **Identificación**: FTIR ≥20 µm; Raman ≥1 µm; Nile Red solo cribado (Erni-Cassola 2017).
- **Normativa**: Reglamento (UE) 2024/1448 (metodología agua de consumo); ISO/TR 21960:2020; ISO 24187:2023. Sin valor límite ambiental internacional vigente **[referencial]**.

## 2. Índices de contaminación por MP

### 2.1 Factor de contaminación (CF)
CF = C_muestra / C_fondo — Hakanson (1980) Hydrobiologia 91:832-848.
Clases: <1 baja · 1–3 moderada · 3–6 considerable · ≥6 muy alta.

### 2.2 Pollution Load Index (PLI)
PLI = (∏ CFᵢ)^(1/n) — Tomlinson et al. (1980) Helgol. Meeresunters. 33:566-575.
PLI>1 = deterioro de la calidad.

### 2.3 Polymer Hazard Index (PHI)
PHI = Σ (Hₚ × Aₚ)/100, con H = puntaje de riesgo químico de Lithner et al. (2011,
ES&T 45:4727-4735): PE 11 · PP 11 · PS 30 · PET 51 · PC 68 · PA 2056 · PU 1156 ·
PVC 5000 · ABS 8552 (PMMA/PLA/PVA/POM **[referencial]**). Clases: <1 bajo · 1–10
moderado · 10–100 alto · 100–1000 peligroso · >1000 extremadamente peligroso.
Polímeros sin puntaje se excluyen y se reporta % UNKNOWN.

### 2.4 Clases de carga de MP **[referencial]**
Agua (part/m³): <100 baja · 100–1000 moderada · 1000–5000 alta · >5000 extremadamente alta.
Sedimento (part/kg): <500 · 500–2000 · 2000–5000 · >5000. Biota (part/ind): <1 · 1–5 · 5–15 · >15.
Rangos construidos de estudios costeros tropicales; ajustar por región antes de publicar.

### 2.5 IICAM (propuesto en OMNIPLAST)
IICAM = Σ wᵢ·Sᵢ / Σ wᵢ con subíndices 0–100:
- S_mp = 100 − percentil de la carga en la distribución del sistema
- S_phi = clip(100 − 20·log₁₀(PHI), 0, 100)
- S_wqi = WQI-NSF
- S_bio, S_tox = opcionales (Shannon/IBR escalados)
Pesos por defecto: mp .30, phi .20, wqi .20, bio .15, tox .15. Incertidumbre: Monte Carlo
(σ=5 por subíndice, 2000 iteraciones). Clases: ≥80 Excelente · 60–79 Buena · 45–59
Aceptable · 30–44 Deteriorada · <30 Crítica.
**Estado**: propuesta metodológica del sistema; validación con panel de expertos pendiente (v2).

## 3. Calidad del agua

### 3.1 WQI-NSF (Brown et al. 1970)
WQI = Σ(qᵢ·wᵢ)/Σwᵢ; 9 parámetros con pesos: OD .17 · coliformes .16 · pH .11 ·
DBO .11 · T .10 · PO₄ .10 · NO₂ .10 · turbidez .08 · SST .07. Curvas q(·) aproximadas
por interpolación de puntos publicados. Clases: >90 excelente · 71–90 buena · 51–70
media · 26–50 mala · <26 muy mala.

### 3.2 Saturación de oxígeno
C_s(T) = 14.652 − 0.41022·T + 0.007991·T² − 7.7774e-5·T³ (APHA/Benson-Krause),
corrección salina empírica ×(1 − 0.0097·S). Hipoxia <2 mg/L.

### 3.3 TSI de Carlson (1977)
TSI(chla) = 30.6 + 9.81·ln(chla) · TSI(PT) = 14.42·ln(PT_µg/L) + 4.15.
Clases: <30 oligotrófico · 40–50 mesotrófico · 60–70 eutrófico · >70 hipereutrófico.

## 4. Estadística (decisiones y supuestos)

- **Conteos de MP**: GLM Binomial Negativa cuando α de sobredispersión > 0 (Cameron &
  Trivedi 2013); nunca ANOVA directo. Offset log-volumen si difiere.
- **Pseudorreplicación**: sitio como intercepto aleatorio (Zuur et al. 2009) o cluster GEE
  (Liang & Zeger 1986).
- **PERMANOVA** (Anderson 2001): pseudo-F por permutación; **siempre acompañado de
  betadisper** (la diferencia puede ser de dispersión).
- **NMDS**: stress de Kruskal <0.05 excelente · <0.10 buena · <0.20 con cautela.
- **Corrección múltiple**: Holm (comparaciones) y Benjamini-Hochberg (ómica).
- **Tamaños de efecto**: d de Cohen (0.2/0.5/0.8), η² (0.06/0.14), ε², δ de Cliff, V de Cramer.
- **Potencia**: n por grupo ≈ 2((z_α+z_β)/d)² para d esperado 0.8, potencia .80.

## 5. Molecular

- **ΔΔCt** (Livak & Schmittgen 2001; MIQE Bustin 2009): FC = (1+E)^(−ΔΔCt); E de la curva
  estándar: E = 10^(−1/pendiente) − 1, válida 90–110% con R²≥0.99. Ct>38 = no detección.
- **Panel de genes**: petase/mhetase (Yoshida 2016 Science 351:1196-1203), alkB (van Beilen
  & Funhoff 2007), laccasa, phaC (Tokiwa 2009), 16S como normalizador. Cebadores
  funcionales: diseñar desde secuencias de referencia (módulo Cebadores); solo los
  universales clásicos van pre-cargados (515F-Y/806R-A: Apprill/Parada 2015; 27F/1492R:
  Lane/Turner; ITS1F/ITS4: White 1990).
- **Diversidad**: Shannon H′=−Σp·ln p (natural); Chao1 = S + F₁²/2F₂; Hill q=0/1/2.
- **DA**: DESeq2-lite (size factors mediana de razones + MWU + BH FDR 5%).
- **HWE**: χ² de genotype counts; Fst Weir & Cockerham (1984) simplificado bialélico.

## 6. Ecotoxicología

- **IBR/2** (Beliaeff & Burgeot 2002; Lam & Gray 2003): zᵢ = (M_tratado − M_ref)/SD_ref;
  IBR/2 = Σ|zᵢ|/2. Perfil por biomarcador (SOD, CAT, GST, LPO).
- **Histopatología** (Bernet et al. 1999): IAL ∝ Σ(prev×severidad×importancia).
- **AMBI** (Borja 2000): (0·%I + 1.5·%II + 3·%III + 4.5·%IV + 6·%V)/10 con lista AZTI.
- **BMWP-Col** (Roldán 2003): >120 muy buena · 101–120 buena · 61–100 aceptable · 36–60
  dudosa · 16–35 crítica · <16 muy crítica.

## 7. Riesgo y One Health

- **HQ = PEC/PNEC** (EU TGD R.10 2003); PNEC referenciales por polímero **[editable]**.
- **SSD lognormal**: HC5 = exp(μ + 1.645σ) de NOEC/EC10 (≥5 especies); PAF(x) = Φ((ln x−μ)/σ).
- **Exposición pesca**: p/kg·pc/día por Monte Carlo (concentración lognormal de biota,
  porción triangular, frecuencia, peso corporal); OMS 2019: evidencia insuficiente para
  límites → comparar escenarios, no "seguridad".

## 8. Geo y remoto

- **Kriging ordinario** con variograma esférico/exponencial/gaussiano (Isaaks &
  Srivastava 1989); reportar nugget/meseta/alcance y validación leave-one-out.
- **Moran I** global (k-NN, permutaciones) y **LISA** HH/LL/HL/LH (Anselin 1995);
  **hotspots Gi*** con confianza 90/95/99 (Getis & Ord 1992).
- **FDI flotantes** (Biermann et al. 2020 Sci. Rep. 10:19836): NIR − (RED + (SWIR−RED)·Δλ);
  valida en campo: detecta escombros genéricos.
- **Drift lagrangiano**: advección viento/corriente + difusión + Stokes con corrección de
  arrastre (Waldschläger & Schüttrumpf 2019) — orden de magnitud, no transporte formal
  (usar OpenDrift/ROMS para publicar).

## 9. Fuentes y recomendaciones

- **Apportionment**: similitud coseno de la huella polímero×morfología contra 8 fuentes
  tipificadas + softmax; resultado = hipótesis contrastable, no prueba.
- **Motor de recomendaciones (AG-100)**: reglas trazables a umbrales (PLI>1, PHI≥100,
  WQI<50, fibra>40%, E. coli>STV 410, HQ>1) con evidencia numérica por fila.
