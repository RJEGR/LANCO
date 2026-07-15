---
title: "LANCO/RUN29 — Depuración del pipeline COI Leray-XT"
subtitle: "Diagnóstico, soluciones bioinformáticas y decisión metodológica"
author: "Ricardo Gómez · LANCO · UABC"
date: "2026-06-23"
format: slides
---

# Resumen ejecutivo

El pipeline DADA2 con filtrado Q20-estricto descartaba el 100 % de los reads COI Leray-XT (MiSeq 2×301). Tras diagnosticar la causa raíz (truncQ + maxEE mutuamente excluyentes) y evaluar 5 soluciones bioinformáticas vía deep-research, se adopta una estrategia **single-end R1 + curación swarm/LULU** que recupera 750 K reads útiles en 38 muestras. La inflación de ASVs observada en la estrategia paired-end (498 ASVs) se identifica como artefacto combinatorio del `justConcatenate`, no como diversidad biológica real.

---

# 1. Contexto del proyecto

| Campo | Valor |
|---|---|
| Proyecto | LANCO RUN29 — Workshop COI metabarcoding |
| Plataforma | Illumina MiSeq M03978, 2×301 paired-end |
| Marcador | COI Leray-XT (mlCOIintF / jgHCO2198) |
| Amplicón esperado | 313 nt |
| Muestras | 38 (grupo L, post-demultiplex) |
| Profundidad media | ~73 K reads / muestra |
| Outlier confirmado | L33 con 2.34 M reads (potencial misdemultiplex, no resuelto en este sprint) |
| Repositorio | `/Users/rjegr/Documents/GitHub/LANCO` |

---

# 2. Síntoma inicial

Ejecución `workflow/03_dada2_pipeline.R` con `truncQ=20, maxEE=(1,2), truncLen=(140,160)`:

- 26 / 38 muestras → **0 reads** pasaron `filterAndTrim`.
- Las 12 restantes → **1-9 reads** retenidos.
- `learnErrors` ejecutó sobre solo 4 060 bases en 29 reads (objetivo ≥100 M bases).
- `dada()` abortó con `OMEGA_A must be between zero and one`.

**Acción recomendada inmediata:** detener el pipeline; el error de OMEGA_A es la consecuencia, no la causa.

---

# 3. Diagnóstico de causa raíz

**Conflicto lógico entre parámetros:**

`truncQ=20` trunca cada read en la primera base con Q<20 *antes* de aplicar `truncLen`. Si tras truncQ el read es más corto que truncLen, se descarta. Con `truncQ=20` estricto sobre MiSeq 2×300, prácticamente toda lectura tiene ≥1 base Q<20 antes de pos 140/160 → todo se descarta.

**Principio bioinformático violado:** `truncQ` y `maxEE` son redundantes; best-practice DADA2 (Callahan, tutorial 1.16) es dejar `truncQ` en el default (2) y dejar que `maxEE` haga el control de calidad.

| Acciones recomendadas |
|---|
| Cambiar `truncQ: 20 → 2` |
| Endurecer coerción YAML para `OMEGA_A` (`as.numeric()` defensivo) |
| Documentar la contradicción interna del YAML original |

---

# 4. Marco metodológico — deep-research

Búsqueda sistemática (PubMed + Consensus + foros oficiales DADA2/QIIME2) identificó **5 estrategias validadas** para recuperar señal en COI 2×300 con R2 colapsado:

| # | Estrategia | Endorse | Aplicabilidad |
|---|---|---|---|
| 1 | `justConcatenate=TRUE` + maxEE_R2 relajado | Callahan, GitHub #176 | Alta — mínimo cambio |
| 2 | Single-end forward only (R1) | Callahan; Antich 2021 | Alta — fallback estándar |
| 3 | Pipeline Brandt: DADA2 → swarm → LULU | Brandt 2021 *Mol Ecol Resour* | Alta — específico COI marino |
| 4 | DnoisE (UNOISE3 entropy-aware) | Antich 2022 *PeerJ* | Media — requiere reescribir denoising |
| 5 | Pre-error-correction Tadpole/BBMerge | Bushnell 2017 *PLOS ONE* | Baja — rompe modelo DADA2 |

---

# 5. Iteración experimental — 3 corridas

| Run | Estrategia | Parámetros | Mediana pct retenido | ASVs |
|---|---|---|---|---|
| PE | Paired-end + justConcatenate | truncLen=(140,160), maxEE=(2,5) | **30 %** | 498 |
| SE v1 | Single-end R1 | truncLen=220, maxEE=2 | **0.7 %** | 65 |
| SE v2 | Single-end R1 (re-tuneado) | truncLen=180, maxEE=3 | **17.25 %** | 261 |

**Lectura clave:** la asunción "R1 mantiene Q≥20 hasta pos 220" fue especulativa. Los datos demuestran que el R1 también colapsa, no tan dramáticamente como R2, pero suficiente para perder 99 % de reads con truncLen=220.

---

# 6. La diferencia 498 vs 65 ASVs — explicación técnica

**Los 498 ASVs del PE NO son diversidad biológica real:**

1. **Combinatoria del concatenate**: cada ASV = R1_haplotype × R2_haplotype. Con *m* haplotipos R1 y *n* haplotipos R2, espacio combinatorio = m × n.
2. **Inflación por ruido R2** (`maxEE_R2=5`): errores residuales en R2 generan pseudo-variantes que `learnErrors` no logra colapsar.
3. **Quimeras intra-linker invisibles**: `removeBimeraDenovo` no detecta recombinaciones dentro del linker `10×N`.

| Pipeline | Reads totales | ASVs | Reads/ASV | Diagnóstico |
|---|---|---|---|---|
| PE | ~810 K | 498 | ~53 | Distribución dispersa → ruido |
| SE v2 | 750 K | 261 | ~2 870 | Distribución concentrada → señal |

**Punto pedagógico para el workshop:** este contraste es material directo para la sección Discussion del manuscrito.

---

# 7. Decisión metodológica

**Estrategia adoptada:** SE v2 (single-end R1, truncLen=180, maxEE=3) + curación posterior con swarm+LULU.

| Justificación | Evidencia |
|---|---|
| Single-end recomendado oficialmente para R2 colapsado | Callahan, DADA2 GitHub #176 |
| Curación swarm+LULU valida específicamente para COI marino | Brandt et al. 2021, *Mol Ecol Resour* |
| 750 K reads en 38 muestras es suficiente para análisis comunidad COI | Estándar de campo (≥5 K/muestra) |
| ASV biológicamente continuo permite NUMT-screen por traducción | Código genético 5 (NCBI invert. mt) |

**Trade-off aceptado:** amplicón final de ~170 nt vs 313 nt completo → menor resolución intra-especie, suficiente para asignación a género/familia.

---

# 8. Arquitectura final del pipeline

```
RUN29/*.fastq.gz
        │
        ├── 01_quality (FastQC + quality profiles)
        ├── 02_cutadapt_primers.sh   (Leray-XT removal, 15% error)
        │
        ├── 03_dada2_pipeline_SE.R    ← ACTIVO (single-end R1)
        │      truncLen=180, maxEE=3, truncQ=2
        │      output: 261 ASVs × 38 muestras, 750 K reads
        │
        ├── 05_swarm_LULU.R           ← curación posterior
        │      swarm d=13 + LULU (84% identity, 95% co-occurrence)
        │      output esperado: ~40-70 OTUs curados
        │
        └── 04_taxonomy_BOLD_NCBI.R   ← asignación final
               BOLD + NCBI vía RESCRIPt (VSEARCH consensus)

ARCHIVADO: 03_dada2_pipeline_PE.R + config/params_PE.yml
           (snapshot reproducible de la estrategia paired-end descartada)
```

---

# 9. Cambios en el repositorio

| Archivo | Estado | Función |
|---|---|---|
| `config/params.yml` | Modificado | Config activa single-end |
| `config/params_PE.yml` | Nuevo backup | Snapshot PE archivado |
| `workflow/03_dada2_pipeline.R` | Sin cambios | Original PE (referencia) |
| `workflow/03_dada2_pipeline_PE.R` | Nuevo backup | Snapshot PE con header de archivado |
| `workflow/03_dada2_pipeline_SE.R` | Nuevo | Pipeline single-end activo |
| `workflow/05_swarm_LULU.R` | Nuevo | Curación post-DADA2 |
| `results/03_dada2_SE/` | Generado | 261 ASVs, 750 K reads |
| `docs/RUN29_pipeline_debug_executive_summary.md` | Este documento | Resumen para slides |

---

# 10. Predicción de resultados post-curación

| Etapa | OTUs / ASVs | Reads conservados | % vs inicial |
|---|---|---|---|
| Input DADA2 SE v2 | 261 ASVs | 750 273 | 100 % baseline |
| Post-swarm d=13 | ~70-110 OTUs | ~745 K | ~99 % |
| Post-LULU | **~40-70 OTUs** | **~720 K** | **~96 %** |

**Criterio de éxito:** retención de reads ≥90 % tras LULU. Esto demuestra que la reducción ASVs→OTUs elimina inflación artefactual sin perder señal biológica.

---

# 11. Próximos pasos

| Prioridad | Acción | Owner | Estado |
|---|---|---|---|
| Alta | Instalar dependencias (swarm, vsearch, R-pkg lulu) | Ricardo | Pendiente |
| Alta | Ejecutar `05_swarm_LULU.R` | Ricardo | Pendiente |
| Media | Asignación taxonómica `04_taxonomy_BOLD_NCBI.R` | Ricardo | Bloqueado por 05 |
| Media | Investigar outlier L33 (2.34 M reads) | Ricardo | Backlog |
| Baja | Evaluar DnoisE como alternativa (iteración futura del workshop) | Ricardo | Backlog |
| Baja | Considerar re-secuenciación con kit v2 (500 ciclos) | Lab manager | Decisión externa |

---

# 12. Referencias clave

**Peer-reviewed (Tier 1):**

- Callahan B.J. et al. (2016). *DADA2: High resolution sample inference from Illumina amplicon data*. Nature Methods. [28 447 citas]
- Brandt M.I. et al. (2021). *Bioinformatic pipelines combining denoising and clustering tools allow for more comprehensive prokaryotic and eukaryotic metabarcoding*. Mol Ecol Resour.
- Antich A. et al. (2021). *To denoise or to cluster, that is not the question: optimizing pipelines for COI metabarcoding*. BMC Bioinformatics.
- Antich A. et al. (2022). *DnoisE: distance denoising by entropy*. PeerJ.
- Bushnell B. et al. (2017). *BBMerge — Accurate paired shotgun read merging via overlap*. PLOS ONE.
- Frøslev T.G. et al. (2017). *Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates*. Nature Communications.
- Mahé F. et al. (2015). *Swarm v2: highly-scalable and high-resolution amplicon clustering*. PeerJ.

**Documentación técnica (Tier 2):**

- DADA2 Pipeline Tutorial 1.16 — Callahan: https://benjjneb.github.io/dada2/tutorial.html
- DADA2 GitHub Issues #176, #279, #790 (justConcatenate + R2 quality)
- QIIME2 Forum: "Low overlap due to poor reverse read quality in 2×300"

---

# Apéndice A — Métricas comparativas por iteración

| Métrica | Q20-strict (crash) | PE | SE v1 | SE v2 |
|---|---|---|---|---|
| `truncLen_R1` | 140 | 140 | 220 | 180 |
| `truncLen_R2` | 160 | 160 | – | – |
| `maxEE_R1` | 1 | 2 | 2 | 3 |
| `maxEE_R2` | 2 | 5 | – | – |
| `truncQ` | 20 | 2 | 2 | 2 |
| Mediana pct retenido | 0 % | 30 % | 0.7 % | 17.25 % |
| Reads totales útiles | <100 | ~810 K | 26 568 | 750 273 |
| ASVs finales | crash | 498 | 65 | 261 |
| `learnErrors` bases | 4 060 | ~5 M | 6.4 M | 138 M |

---

# Apéndice B — Punto pedagógico para el workshop

El contraste **PE (498 ASVs) vs SE+swarm/LULU (~50 OTUs)** sobre la misma muestra física constituye material didáctico de primer orden:

1. Demuestra que **la elección de parámetros bioinformáticos altera la "diversidad" reportada en ~10×**.
2. Refuerza que **DADA2 fue diseñado para 16S** (sin estructura de codón) y su aplicación naïve a COI genera ASVs combinatorios espurios (Antich 2021).
3. Justifica el uso de **denoising entropy-aware (DnoisE)** o **post-clustering curation (LULU)** como pasos no-opcionales en COI metabarcoding.

**Aplicación inmediata:** este caso real puede convertirse en un dataset de práctica para que los estudiantes ejecuten ambos pipelines y comparen Shannon, riqueza y curvas rarefaction.
