# LANCO / RUN29 — Pipeline COI metabarcoding
### Resumen técnico para *slides ejecutivos* · Workshop "Microbiome amplicon workshop design"

---

## Resumen ejecutivo

Construimos un workflow reproducible DADA2 + RESCRIPt para procesar 38 muestras COI Leray-XT del run MiSeq M03978 de LANCO, con trimming Q20-reality data-driven y clasificación taxonómica híbrida BOLD/NCBI. Diagnóstico empírico sobre los 76 FASTQ reveló dos sorpresas: el run contiene tres librerías multiplexadas (no solo COI) y la calidad colapsa antes de lo asumido, lo que llevó a usar `mergePairs(justConcatenate=TRUE)` como decisión informada por datos. El workflow queda listo para ejecución local con dos entornos conda parametrizados desde un único `params.yml`.

---

## Slide 1 — Contexto y objetivo

- **Proyecto:** LANCO/RUN29 (Laboratorio Nacional CONAHCYT de Oceanografía)
- **Doble propósito:** producción + material reproducible para taller
- **Input:** 60 muestras paired-end, MiSeq M03978, 2×301 nt
- **Marker declarado:** COI Leray-XT (mlCOIintF / jgHCO2198, ~313 pb)
- **Output esperado:** ASV table + taxonomía + diagnósticos de calidad

---

## Slide 2 — Composición real del run (hallazgo crítico)

Inspección automática por regex de primer sobre 2 000 reads/muestra reveló **tres librerías multiplexadas** en el mismo lane:

| Prefijo | n | Amplicón detectado | % reads con primer | Acciones recomendadas |
|---|---|---|---|---|
| L1–L38 | 38 | COI Leray-XT (mlCOIintF) | 92–95 % | **Procesar en este workflow** |
| T1, T2 | 2 | ITS1F (hongos) | ~64 % | Excluir — n insuficiente |
| T3–T31 | 20 | 16S V3–V4 (341F) | 75–95 % | Excluir — workflow paralelo futuro |

Decisión: **alcance restringido a COI**. Mapping filtrado en `metadata/samples_RUN29_COI.tsv` (38 filas); original respaldado.

---

## Slide 3 — Arquitectura del pipeline

```
RUN29/L*.fastq.gz
    │
    ├─ 00_setup_envs.sh ── crea lanco_coi + qiime2-rescript
    │
    ├─ 01_inspect_reads.R ──── quality profiles + primer hits
    │
    ├─ 02_cutadapt_primers.sh ─ remueve Leray-XT (linked adapters)
    │                           --discard-untrimmed
    │
    ├─ 03_dada2_pipeline.R ──── filterAndTrim → learnErrors → dada
    │                           → mergePairs(justConcatenate=T)
    │                           → removeBimeraDenovo
    │                           → length filter (305–315 nt)
    │
    ├─ 04a_build_BOLD_db_RESCRIPt.sh  (UNA vez, 2–8 h)
    ├─ 04b_build_NCBI_db_RESCRIPt.sh  (UNA vez, 1–4 h)
    │
    └─ 04_taxonomy_BOLD_NCBI.R ── VSEARCH consensus
                                  BOLD primario + NCBI rescate
                                  → tax_consensus.tsv
```

Parametrización centralizada en `config/params.yml` (leído por bash con `yq` y por R con `yaml::read_yaml`).

---

## Slide 4 — Diagnóstico de calidad (n=38, sandbox Python)

Métricas agregadas sobre los 76 FASTQ:

| Métrica | Valor |
|---|---|
| Reads/muestra (mediana) | 62 764 |
| Min reads | 28 046 (L29) |
| Max reads | **2 368 881 (L33, 37.7× la mediana)** |
| Mediana % primer fwd presente en R1 | 93.0 % |
| Mediana % primer rev presente en R2 | 93.2 % |
| Mediana Q≥20 R1 sostenida hasta | pos ~130 (raw) |
| Mediana Q≥20 R2 sostenida hasta | pos ~145 (raw) |

Outputs persistidos en `results/01_quality/` (PDFs + TSVs).

---

## Slide 5 — Conflicto Q20 strict vs merge biológico

| Condición | truncLen R1 + R2 | Mínimo merge (313 + 12) | Margen |
|---|---|---|---|
| Q20 strict según mediana (104+119) | **223 nt** | 325 nt | **−102 nt** ❌ |
| Configuración inicial (220+150) | 370 nt | 325 nt | +45 nt ✓ |
| **Decisión Q20-reality (140+160)** | **300 nt** | 325 nt | **−25 nt** → justConcatenate |

**Decisión Ricardo (2026-06-22):** honrar Q20 real aunque eso impida merge biológico. Usar `mergePairs(justConcatenate=TRUE)` → ASV = R1 + 10×N + revComp(R2) = 310 nt.

---

## Slide 6 — Parámetros finales (`params.yml`)

| Bloque | Parámetro | Valor | Justificación |
|---|---|---|---|
| cutadapt | error_rate | 0.15 | Degeneración alta de COI |
| cutadapt | discard_untrimmed | true | Off-target = ruido para DADA2 |
| filterAndTrim | truncLen R1 | 140 (post-cutadapt) | Q20-reality empírico |
| filterAndTrim | truncLen R2 | 160 (post-cutadapt) | Q20-reality empírico |
| filterAndTrim | truncQ | 20 | Decisión Ricardo |
| filterAndTrim | maxEE | (1, 2) | Estricto |
| dada2 | pool | pseudo | Balance sensibilidad/costo |
| mergePairs | justConcatenate | true | Inevitable con Q20-reality |
| coi_filtering | min/max_length | 305 / 315 nt | 140 + 10N + 160 ± 5 |
| coi_filtering | NUMT translation check | **false** | Linker N rompe marco |

---

## Slide 7 — Taxonomía: BOLD + NCBI vía RESCRIPt

Reemplaza MIDORI2 con DBs construidas explícitamente:

| DB | Filtros aplicados | Output |
|---|---|---|
| **BOLD** (tutorial qiime2 forum 16129) | Marker COI-5P; taxa Metazoa+Algae+Fungi; longitud 250–1600; max 5 N | `bold-coi-{seqs,tax,classifier}.qza` + FASTA DADA2-format |
| **NCBI** (tutorial qiime2 forum 16500) | Entrez query con SLEN[250:1600]; excluye env/unclassified; 7 ranks | `ncbi-coi-{seqs,tax,classifier}.qza` + FASTA DADA2-format |

**Estrategia de clasificación:** VSEARCH `--usearch_global` (default) → tolera el gap N central del ASV concatenado mejor que Naive Bayes. BOLD primario; ASVs sin Género → NCBI rescate. Consenso por rango con `min_consensus = 0.51`.

---

## Slide 8 — Entregables

| Artefacto | Ubicación | Propósito |
|---|---|---|
| Workflow scripts | `workflow/0[0-4]*.{R,sh}` | 7 scripts numerados, ejecución secuencial |
| Configuración | `config/params.yml` | Única fuente de parámetros |
| Mapping COI | `metadata/samples_RUN29_COI.tsv` | 38 muestras + slots para vars ambientales |
| Diagnósticos | `results/01_quality/*.{pdf,tsv}` | Quality profiles + primer hits + read counts |
| Documentación | `workflow/README.md` | Tabla de orden de ejecución + acciones |
| Memoria persistente | `memory/*.md` | 5 entradas (user, project×2, reference, feedback) |

---

## Slide 9 — Riesgos identificados

| Riesgo | Severidad | Acciones recomendadas |
|---|---|---|
| L33 con 37.7× la mediana de reads | Medio | Documentar en `track_reads.tsv`; revisar si distorsiona learnErrors |
| Retención post-Q20 desconocida | Medio | Primer test real validará si maxEE descarta demasiado |
| ASV concatenado (310 nt con linker N) | Bajo | VSEARCH tolera gaps; NUMT screen desactivado por diseño |
| Bottom 5 muestras con 28–35K reads | Bajo | Suficiente para inferencia DADA2; monitorear |
| QIIME2 no nativo en Apple Silicon | Bajo | Setup script ya detecta y fuerza Rosetta 2 |

---

## Slide 10 — Próximos pasos

| # | Acción | Tiempo estimado |
|---|---|---|
| 1 | `bash workflow/00_setup_envs.sh` | 20–40 min |
| 2 | Llenar columnas ambientales en `samples_RUN29_COI.tsv` | manual |
| 3 | `conda activate lanco_coi && bash workflow/02_cutadapt_primers.sh` | 5–15 min |
| 4 | `Rscript workflow/03_dada2_pipeline.R` (primer test real Q20-reality) | 30–60 min |
| 5 | En paralelo: `conda activate qiime2-rescript && bash workflow/04a_build_BOLD_db_RESCRIPt.sh` | 2–8 h |
| 6 | `bash workflow/04b_build_NCBI_db_RESCRIPt.sh` (configurar `NCBI_API_KEY`) | 1–4 h |
| 7 | `Rscript workflow/04_taxonomy_BOLD_NCBI.R` | 30–90 min |
| 8 | Construir phyloseq + diversidad/ordenación (fuera de scope actual) | siguiente sprint |

---

## Slide 11 — Métricas del proyecto

| Métrica | Valor |
|---|---|
| Scripts producidos | 7 |
| Líneas de código (R + bash) | ~1 400 |
| Parámetros centralizados | 35+ en `params.yml` |
| Tareas completadas | 21 (todas) |
| FASTQ inspeccionados | 76 |
| Decisiones técnicas documentadas | 7 (Q20-reality, justConcatenate, BOLD-first, etc.) |
| Memorias persistentes guardadas | 5 |
| Tiempo estimado de ejecución full pipeline (post-DB) | 1–2 h por run |

---

## Apéndice — Decisiones de diseño clave

1. **Solo COI por ahora** — n=38 vs alternativas multi-amplicón (3× trabajo).
2. **DADA2 en R sobre nf-core** — control pedagógico paso a paso para workshop.
3. **RESCRIPt sobre MIDORI2** — versionado explícito de scope taxonómico + query Entrez.
4. **Q20-reality + justConcatenate** sobre Q20-nominal + merge — honrar criterios de calidad por encima de conveniencia técnica.
5. **VSEARCH consensus** sobre Naive Bayes — robustez a gaps N centrales.
6. **BOLD primario + NCBI rescate** — mejor cobertura metazoos costeros + rescate de raros.
7. **Outliers documentados, no preprocesados automáticamente** — transparencia sobre intervención.
