# LANCO / RUN29 — Workflow COI metabarcoding

**Resumen ejecutivo.** Procesamiento de 60 muestras COI Leray-XT (MiSeq 2×301) de los grupos L (n=38) y T (n=22). El pipeline cubre cinco etapas: inspección de calidad, remoción de primers con cutadapt, inferencia de ASVs con DADA2 (con filtro NUMT por traducción), clasificación taxonómica con MIDORI2 y verificación contra BOLD. Diseñado como material reproducible para el *Microbiome amplicon workshop*.

## Estructura del repositorio

```
LANCO/
├── RUN29/                      # FASTQ crudos (no tocar)
├── config/
│   └── params.yml              # única fuente de parámetros
├── metadata/
│   └── samples_RUN29.tsv       # mapping file (rellenar columnas ambientales)
├── workflow/
│   ├── 01_inspect_reads.R
│   ├── 02_cutadapt_primers.sh
│   ├── 03_dada2_pipeline.R
│   ├── 04_taxonomy_BOLD_NCBI.R
│   └── README.md               # este archivo
├── results/
│   ├── 01_quality/
│   ├── 02_cutadapt/
│   ├── 03_dada2/
│   └── 04_taxonomy/
└── logs/
```

## Alcance de RUN29

La inspección automática de los FASTQ reveló que RUN29 contiene tres librerías multiplexadas:

- **L1–L38 (n=38)**: COI Leray-XT  ← procesadas por este workflow
- **T1, T2**: ITS1F (hongos)         ← excluidas
- **T3–T31 (n=20)**: 16S V3–V4 (341F) ← excluidas

El mapping file `metadata/samples_RUN29_COI.tsv` ya está filtrado a las 38 muestras COI. El original con los 60 sampleos está respaldado en `samples_RUN29.tsv.bak_YYYYMMDD`.

## Orden de ejecución y acciones recomendadas

| Paso | Script | Entrada | Salida clave | Acciones recomendadas |
|------|--------|---------|--------------|----------------------|
| 1 | `01_inspect_reads.R` | `RUN29/L*.fastq.gz` | `quality_profile_R{1,2}.pdf`, `primer_hits.tsv` | Confirmar Q drop y presencia de primer (>90 %) en muestras L |
| 2 | `02_cutadapt_primers.sh` | `RUN29/L*.fastq.gz` | `results/02_cutadapt/L*.fastq.gz`, `cutadapt_summary.tsv` | Verificar `pct_passed > 80 %` |
| 3 | `03_dada2_pipeline.R` (Q20 strict) | `results/02_cutadapt/` | `asvs.fasta`, `asv_table.tsv`, `track_reads.tsv`, `quality_postfilter_R*.pdf` | Si `pct_q20 < 40 %` en >5 muestras: relajar `truncQ` a 15 |
| 4a | `04a_build_BOLD_db_RESCRIPt.sh` | (red) | `db/rescript/bold-coi-{seqs,tax,classifier}.qza`, `.fasta`, `.tsv` | Ejecutar UNA SOLA VEZ. 2–8 h |
| 4b | `04b_build_NCBI_db_RESCRIPt.sh` | (red) | `db/rescript/ncbi-coi-{seqs,tax,classifier}.qza`, `.fasta`, `.tsv` | Ejecutar UNA SOLA VEZ. 1–4 h. Usar `NCBI_API_KEY` |
| 4 | `04_taxonomy_BOLD_NCBI.R` | `asvs.fasta` + DBs | `tax_BOLD.tsv`, `tax_NCBI.tsv`, `tax_consensus.tsv`, `assignment_summary.tsv` | Comparar % asignación BOLD vs NCBI por rango — material de discusión del workshop |

## Dependencias

R (≥ 4.3):

```r
install.packages(c("BiocManager", "yaml", "dplyr", "readr", "tibble",
                   "ggplot2", "rprojroot"))
BiocManager::install(c("dada2", "Biostrings", "ShortRead"))
```

Dos entornos conda separados (cutadapt no necesita QIIME2):

```bash
# Entorno 1 — pipeline DADA2 + cutadapt + vsearch
mamba create -n lanco_coi -c bioconda -c conda-forge \
  cutadapt=4.* vsearch yq biopython python=3.11
mamba activate lanco_coi

# Entorno 2 — RESCRIPt para construir DBs (solo pasos 04a/04b)
mamba create -n qiime2-rescript \
  -c qiime2 -c conda-forge -c bioconda \
  qiime2 q2-rescript python=3.10
mamba activate qiime2-rescript
pip install git+https://github.com/bokulich-lab/RESCRIPt.git
```

Bases de datos: construidas con `04a_build_BOLD_db_RESCRIPt.sh` y `04b_build_NCBI_db_RESCRIPt.sh`. No requieren descarga manual previa más allá de un `NCBI_API_KEY` opcional para acelerar Entrez.

## Decisiones del workshop

1. **Trimming Q20 estricto** (`truncQ=20`, `maxEE=(1,2)`): cada read se trunca en la primera base con Q<20. Los `truncLen` (R1=220, R2=150 post-cutadapt) se eligieron empíricamente del colapso de calidad observado en muestras L.
2. **Pseudo-pooling en DADA2** para balance sensibilidad/costo en eDNA marino heterogéneo.
3. **Filtro NUMT por traducción** (código 5, Invertebrate Mitochondrial) contra pseudogenes.
4. **Filtro de longitud 299–320 nt** centrado en Leray (313 nt).
5. **DBs vía RESCRIPt** (BOLD + NCBI) en lugar de MIDORI2: control total del scope taxonómico (Metazoa/Algae/Fungi para LANCO) y query Entrez explícito para NCBI. Reproducible y citable.
6. **VSEARCH consensus** como clasificador por defecto: robusto a gaps de COI, no requiere reentrenamiento. Naive Bayes (sklearn) y `assignTaxonomy()` disponibles como alternativas.
7. **BOLD primero, NCBI para rescate**: BOLD tiene mejor cobertura de metazoos costeros; NCBI rescata clados raros (microeucariotas, parásitos).

## Estado actual

- [x] Estructura del proyecto creada
- [x] Mapping file filtrado a COI (`samples_RUN29_COI.tsv`, n=38)
- [x] `params.yml` con valores Leray-XT y Q20 strict
- [x] Scripts 01, 02, 03 (Q20), 04a, 04b, 04 listos
- [ ] Rellenar columnas `site`, `depth_m`, `date_sampled`, `lat`, `lon` en mapping file
- [ ] Configurar entornos conda (lanco_coi + qiime2-rescript)
- [ ] Ejecutar 04a + 04b una vez para generar DBs
- [ ] Primera corrida 01→03→04 y ajustar `truncLen` con datos reales

## Backups

Sigue la regla del proyecto: antes de sobrescribir cualquier resultado, copia con sufijo `.bak_YYYYMMDD`. Los FASTQ crudos en `RUN29/` se tratan como read-only.
