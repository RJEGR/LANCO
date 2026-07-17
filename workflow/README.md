# LANCO / RUN29 — Workflow COI metabarcoding

**Resumen ejecutivo.** Pipeline reproducible para 38 muestras COI Leray-XT (MiSeq 2×301) del proyecto LANCO/CICESE, ejecutado en el cluster LUSTRE con `qiime2-2026`. Procesa desde FASTQ crudos hasta tabla de OTUs curados (swarm+LULU) con asignación taxonómica BOLD+NCBI, y sirve como material del *Microbiome amplicon workshop*.

## Estructura del repositorio

```
LANCO/
├── RUN29/                             # FASTQ crudos (read-only)
├── config/
│   ├── params.yml                     # única fuente de parámetros (SE activo)
│   └── params_PE.yml                  # snapshot PE (referencia histórica)
├── metadata/
│   ├── samples_RUN29_COI.tsv          # mapping activo (38 muestras L)
│   └── samples_RUN29.tsv              # mapping completo con target_detected
├── workflow/
│   ├── 00_setup_envs.sh
│   ├── 01_inspect_reads.R
│   ├── 02_cutadapt_primers.sh
│   ├── 03_dada2_pipeline_SE.R         # estrategia ACTIVA
│   ├── 03_dada2_pipeline.R            # PE con justConcatenate (deprecado)
│   ├── 03_dada2_pipeline_PE.R         # snapshot PE
│   ├── 04a_build_BOLD_db_SciLifeLab.sh   # ★ nuevo (reemplaza RESCRIPt BOLD)
│   ├── 04a_build_BOLD_db_RESCRIPt.sh     # obsoleto (get-bold-data ya no existe)
│   ├── 04b_build_NCBI_db_RESCRIPt.sh
│   ├── 04_build_dbs_RESCRIPt.slurm       # ★ wrapper SLURM: BOLD + NCBI
│   ├── 04_taxonomy_BOLD_NCBI.R
│   ├── 05_classify_ASVs_BOLD.slurm       # ★ finaliza DB + clasifica ASVs
│   ├── 05_swarm_LULU.R
│   └── README.md                      # este archivo
├── db/rescript/                       # DBs .qza + FASTA/TSV exportados
├── results/
│   ├── 01_quality/
│   ├── 02_cutadapt/
│   ├── 03_dada2_SE/                   # output activo
│   ├── 04_taxonomy/
│   └── 05_swarm_lulu/
└── logs/
```

## Alcance de RUN29

RUN29 contiene tres librerías multiplexadas:

- **L1–L38 (n=38)**: COI Leray-XT ← procesadas por este workflow
- **T1, T2**: ITS1F (hongos) ← excluidas
- **T3–T31 (n=20)**: 16S V3–V4 (341F) ← excluidas

Mapping activo: `metadata/samples_RUN29_COI.tsv`. Original respaldado en `samples_RUN29.tsv.bak_YYYYMMDD`.

## Orden de ejecución

| Paso | Script | Entrada | Salida clave | Acciones recomendadas |
|------|--------|---------|--------------|----------------------|
| 0 | `00_setup_envs.sh` | (red) | env `qiime2-2026` + RESCRIPt + R packages | Ejecutar UNA vez. `verify` para re-chequear |
| 1 | `01_inspect_reads.R` | `RUN29/L*.fastq.gz` | `quality_profile_R{1,2}.pdf`, `primer_hits.tsv` | Confirmar Q drop y presencia de primer (>90 %) |
| 2 | `02_cutadapt_primers.sh` | `RUN29/L*.fastq.gz` | `results/02_cutadapt/L*.fastq.gz` | Verificar `pct_passed > 80 %` |
| 3 | `03_dada2_pipeline_SE.R` | `results/02_cutadapt/` | `asvs.fasta`, `asv_table.tsv`, `track_reads.tsv` | Retención esperada ≥70 %; si <40 %, ajustar `truncLen_R1` |
| 4a | `04a_build_BOLD_db_SciLifeLab.sh` | (red, 416 MB) | `bold-coi-{seqs,tax,classifier}.qza` + FASTA/TSV | UNA vez. 30–90 min. Descarga en nodo login si sin DNS |
| 4b | `04b_build_NCBI_db_RESCRIPt.sh` | (red Entrez) | `ncbi-coi-{seqs,tax,classifier}.qza` + FASTA/TSV | UNA vez. 15–40 min con `NCBI_API_KEY`, 1–4 h sin ella |
| 4 (SLURM) | `04_build_dbs_RESCRIPt.slurm` | scripts 4a+4b | Todos los `.qza` finales | Wrapper para lanzar 4a+4b en cluster: `sbatch workflow/04_build_dbs_RESCRIPt.slurm` |
| 5 (SLURM) | `05_classify_ASVs_BOLD.slurm` | `asvs.fasta` + BOLD DB | `results/04_taxonomy/tax_BOLD.tsv` | Finaliza DB (dereplicate) + clasifica con VSEARCH consensus |
| 5R | `04_taxonomy_BOLD_NCBI.R` | ASVs + FASTAs exportados | `tax_consensus.tsv`, `assignment_summary.tsv` | Consenso híbrido BOLD (primario) + NCBI (rescate). Material de discusión workshop |
| 6 | `05_swarm_LULU.R` | ASVs + tabla | `curated_otu_table.tsv` | Curación post-DADA2 (Brandt 2021) |

## Ejecución en cluster LUSTRE (CICESE)

**Env conda:** `qiime2-2026` en `/LUSTRE/apps/Anaconda/2026/miniconda3/envs/`. Los SLURM wrappers exportan el PATH y activan el env automáticamente. RESCRIPt debe estar instalado dentro de `qiime2-2026` (`pip install git+https://github.com/bokulich-lab/RESCRIPt.git`).

**yq (mikefarah v4.18+):** requerido para leer `params.yml` desde bash. Instalación sin sudo:

```bash
mkdir -p ~/.local/bin
curl -fsSL https://github.com/mikefarah/yq/releases/download/v4.44.3/yq_linux_amd64 -o ~/.local/bin/yq
chmod +x ~/.local/bin/yq
```

**Lanzamiento SLURM:**

```bash
cd /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/LANCO

# Construir DBs (BOLD SciLifeLab + NCBI vía RESCRIPt)
sbatch --export=ALL,NCBI_API_KEY=xxxxxxxx workflow/04_build_dbs_RESCRIPt.slurm

# Clasificar ASVs vs BOLD
sbatch workflow/05_classify_ASVs_BOLD.slurm
# Cambiar input de ASVs:
sbatch --export=ALL,ASV_FASTA=results/03_dada2/asvs.fasta workflow/05_classify_ASVs_BOLD.slurm
```

**Restricciones típicas HPC LUSTRE**: los nodos de cómputo no resuelven DNS externo (Figshare bloqueado, NCBI parcialmente permitido). Para `04a` descarga los 2 archivos de SciLifeLab en el nodo login primero — el script tiene `download_check` con MD5 y hará `skip` cuando corra en el nodo:

```bash
mkdir -p db/rescript/_work_bold_scilifelab
cd db/rescript/_work_bold_scilifelab
curl -fSL -o coidb.clustered.fasta.gz        https://ndownloader.figshare.com/files/65310486
curl -fSL -o coidb.qiime2.info.exclNA.tsv.gz https://ndownloader.figshare.com/files/65310780
md5sum coidb.*   # verificar contra los MD5 en 04a_build_BOLD_db_SciLifeLab.sh
```

## Dependencias

**R (≥ 4.3):**

```r
install.packages(c("BiocManager", "yaml", "dplyr", "readr", "tibble",
                   "ggplot2", "rprojroot"))
BiocManager::install(c("dada2", "Biostrings", "ShortRead"))
devtools::install_github("tobiasgf/lulu")   # para paso 05
```

**Env conda** (creado por `00_setup_envs.sh`, o usar el `qiime2-2026` compartido del cluster):

```bash
mamba create -n qiime2-2026 -c qiime2 -c conda-forge -c bioconda \
  qiime2 python=3.12 cutadapt=4.* vsearch swarm biopython
mamba activate qiime2-2026
pip install git+https://github.com/bokulich-lab/RESCRIPt.git
```

## Decisiones del workshop

1. **Trimming SINGLE-END R1** (activo desde 2026-06-23). El R2 del MiSeq 2×301 colapsa <Q20 después del ciclo ~180 (kit v3 endémico) y un intento previo con `justConcatenate` retuvo solo ~30 % de reads. La estrategia actual usa solo R1 con `truncLen=180`, `maxEE=3`. Trade-off: perdemos los ~100 nt 3' del amplicón (313 nt full-length) pero el ASV es biológicamente continuo → **NUMT-screen por traducción reactivado** (código genético 5, invertebrado).
2. **Pseudo-pooling en DADA2** para balance sensibilidad/costo en eDNA marino heterogéneo.
3. **Filtro de longitud 170–185 nt** centrado en el ASV single-end (~180 nt post-cutadapt).
4. **DB BOLD vía SciLifeLab v6** (DOI [10.17044/scilifelab.20514192.v6](https://doi.org/10.17044/scilifelab.20514192.v6), derivada del BOLD Data Package 15-May-2026, ~24 M secuencias). Sustituye la ruta anterior con `qiime rescript get-bold-data`, que fue eliminada en releases recientes de RESCRIPt tras romperse la API pública de BOLD. Ventajas: pre-curada, clusterizada al 100 % por BOLD BIN, formato QIIME2 nativo, reproducible por DOI+MD5.
5. **DB NCBI vía RESCRIPt** (`get-ncbi-data`) con query Entrez explícito filtrado por longitud 250–1600 nt, excluyendo `environmental samples` y `unclassified`. Reproducible y citable.
6. **VSEARCH consensus** como clasificador por defecto (params: `perc_identity=0.85`, `maxaccepts=10`, `min_consensus=0.51`): robusto a gaps de COI, no requiere reentrenar. Naive Bayes (sklearn) disponible como alternativa.
7. **BOLD primero, NCBI para rescate**: BOLD tiene mejor cobertura de metazoos costeros; NCBI rescata clados raros (microeucariotas, parásitos).
8. **Curación post-DADA2 con swarm + LULU** (Brandt 2021): DADA2 puede sobre-llamar variantes en marcadores codificantes; swarm clusteriza localmente y LULU elimina ASVs espurios por co-ocurrencia.

## Troubleshooting

| Síntoma | Causa | Solución |
|---|---|---|
| `yq: command not found` en SLURM | Los jobs no cargan `~/.bashrc` | El SLURM ya exporta `~/.local/bin` al PATH; asegúrate de tener `yq` instalado ahí |
| `Error: unknown command ".foo" for yq` | yq v4 < 4.18 requiere `eval` | Actualiza a v4.18+: `curl ...yq_linux_amd64 -o ~/.local/bin/yq && chmod +x ~/.local/bin/yq` |
| `mkdir: Permission denied 'logs'` + `workflow/04a...: No such file or directory` | SLURM copia el script al spool; `$BASH_SOURCE` NO apunta al repo | Ya arreglado: los SLURM usan `$SLURM_SUBMIT_DIR` |
| `gzip: invalid option -- 'k'` | gzip 1.5 no soporta `-k` | Ya arreglado: `04a` usa `zcat foo.gz > foo` en su lugar |
| `Could not resolve host: ndownloader.figshare.com` | Nodos de cómputo sin DNS externo | Descarga manual en nodo login (ver sección "Ejecución en cluster") |
| `Error: QIIME 2 plugin 'rescript' has no action 'get-bold-data'` | RESCRIPt eliminó la acción; API BOLD cambió | Ya migrado a `04a_build_BOLD_db_SciLifeLab.sh` |
| `ChunkedEncodingError: Response ended prematurely` en NCBI | Fallo transitorio de red durante Entrez fetch | Relanza el job; NCBI típicamente funciona al 2º intento. Corre entre 21:00–05:00 ET |

## Recursos SLURM

Configuración por defecto en los wrappers (ajustable con `#SBATCH` directives):

- `--partition=cicese`
- `--cpus-per-task=8` (empata `params.yml → compute.cores_total`)
- `--mem=64G`
- `--time=06:00:00` para `04`, `04:00:00` para `05`
- `--mail-user=rgomez41@uabc.edu.mx` — notifica BEGIN/END/FAIL

## Estado actual (2026-07-15)

- [x] Estructura del proyecto creada
- [x] Mapping filtrado a COI (`samples_RUN29_COI.tsv`, n=38)
- [x] `params.yml` con estrategia SINGLE-END R1 activa
- [x] Scripts 00–05 listos (incluye SLURM wrappers)
- [x] DB BOLD SciLifeLab intermedios generados (`bold-len-filt-*.qza`, ~340 MB c/u)
- [ ] Finalizar DB BOLD (dereplicate + export) — corre `05_classify_ASVs_BOLD.slurm`
- [ ] Ejecutar `04b` NCBI (fallo transitorio previo por ChunkedEncodingError)
- [ ] Clasificación completa de ASVs con VSEARCH consensus
- [ ] Rellenar columnas `site`, `depth_m`, `date_sampled`, `lat`, `lon` en mapping
- [ ] Rescate NCBI para Unassigned (04_taxonomy_BOLD_NCBI.R)

## Backups

Regla del proyecto: antes de sobrescribir cualquier resultado, copia con sufijo `.bak_YYYYMMDD`. Los FASTQ crudos en `RUN29/` se tratan como read-only.

## Referencias

- BOLD SciLifeLab v6 — Sundh J. (2026). *COI reference sequences from BOLD DB*. Swedish Museum of Natural History. DOI: [10.17044/scilifelab.20514192.v6](https://doi.org/10.17044/scilifelab.20514192.v6).
- BOLD Data Packages — [boldsystems.org/data/data-packages](https://boldsystems.org/data/data-packages/).
- RESCRIPt — Robeson MS et al. (2021). *RESCRIPt: Reproducible sequence taxonomy reference database management*. PLoS Comput Biol. 17(11):e1009581.
- Brandt MI et al. (2021). *Bioinformatic pipelines combining denoising and clustering tools allow accurate representation of biodiversity in eDNA*. Mol Ecol Resour. 21(6):1904–1921.
- Wangensteen OS et al. (2018). *DNA metabarcoding of littoral hard-bottom communities: high diversity and database gaps revealed by two molecular markers*. PeerJ 6:e4705. (primers Leray-XT)
