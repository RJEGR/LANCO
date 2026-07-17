#!/usr/bin/env bash
# =============================================================================
# 05_classify_ASVs_BOLD.slurm
# Paso 1: Finaliza la DB BOLD SciLifeLab (dereplicate + export FASTA/TSV)
# Paso 2: Clasifica los ASVs de DADA2 vs BOLD usando classify-consensus-vsearch
# -----------------------------------------------------------------------------
# Uso:
#   cd /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/LANCO
#   sbatch workflow/05_classify_ASVs_BOLD.slurm
#
# Override del FASTA de ASVs (por defecto SE):
#   sbatch --export=ALL,ASV_FASTA=results/03_dada2/asvs.fasta \
#          workflow/05_classify_ASVs_BOLD.slurm
# =============================================================================
#SBATCH --job-name=lanco_classify_bold
#SBATCH --output=logs/05_classify_%x_%j.out
#SBATCH --error=logs/05_classify_%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=cicese
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rgomez41@uabc.edu.mx

set -euo pipefail

# ---- Env --------------------------------------------------------------------
export PATH=$HOME/.local/bin:/LUSTRE/apps/Anaconda/2026/miniconda3/condabin:/LUSTRE/apps/Anaconda/2026/miniconda3/bin:$PATH
# shellcheck disable=SC1091
source activate qiime2-2026

ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$ROOT" || { echo "[FATAL] no puedo entrar a $ROOT"; exit 1; }
mkdir -p logs results/04_taxonomy

echo "[$(date -Iseconds)] Job arrancó en $(hostname):$(pwd)"
echo " JobID:  ${SLURM_JOB_ID:-N/A}"
echo " qiime:  $(command -v qiime)"
echo " yq:     $(command -v yq)"

# ---- Parámetros -------------------------------------------------------------
PARAMS="config/params.yml"
PID=$(yq       '.taxonomy.classification.perc_identity'  "$PARAMS")
MAXACC=$(yq    '.taxonomy.classification.maxaccepts'     "$PARAMS")
MINCONS=$(yq   '.taxonomy.classification.min_consensus'  "$PARAMS")
THREADS=$(yq   '.taxonomy.classification.threads'        "$PARAMS")
DEREP_MODE=$(yq '.taxonomy.bold.dereplicate_mode'        "$PARAMS")
CORES=$(yq     '.compute.cores_total'                    "$PARAMS")

SEQS_QZA=$(yq  '.taxonomy.bold.seqs_qza'   "$PARAMS")
TAX_QZA=$(yq   '.taxonomy.bold.tax_qza'    "$PARAMS")
SEQS_FA=$(yq   '.taxonomy.bold.seqs_fasta' "$PARAMS")
TAX_TSV=$(yq   '.taxonomy.bold.tax_tsv'    "$PARAMS")

# FASTA de ASVs a clasificar. Default: single-end (estrategia activa).
ASV_FASTA="${ASV_FASTA:-results/03_dada2_SE/asvs.fasta}"
[[ -s "$ASV_FASTA" ]] || { echo "[FATAL] no existe $ASV_FASTA"; exit 1; }
echo " ASVs a clasificar: $ASV_FASTA ($(grep -c '^>' "$ASV_FASTA") secuencias)"

WORK="db/rescript/_work_bold_scilifelab"
TAXO_OUT="results/04_taxonomy"
mkdir -p "$WORK" "$TAXO_OUT"

echo "=========================================================="
echo " PASO 1 — Finalizar BOLD DB (dereplicate + export)"
echo "=========================================================="

# ---- 1a) Dereplicate (retoma desde bold-len-filt-*.qza) --------------------
if [[ -s "$SEQS_QZA" && -s "$TAX_QZA" ]]; then
  echo "[skip] $SEQS_QZA y $TAX_QZA ya existen"
else
  echo "[run] dereplicate (mode=$DEREP_MODE, threads=$CORES)"
  time qiime rescript dereplicate \
    --i-sequences "$WORK/bold-len-filt-seqs.qza" \
    --i-taxa      "$WORK/bold-len-filt-tax.qza" \
    --p-mode      "$DEREP_MODE" \
    --p-threads   "$CORES" \
    --o-dereplicated-sequences "$SEQS_QZA" \
    --o-dereplicated-taxa      "$TAX_QZA"
fi

# ---- 1b) Export FASTA (DADA2 format) + TSV (para el R downstream) ---------
if [[ -s "$SEQS_FA" && -s "$TAX_TSV" ]]; then
  echo "[skip] $SEQS_FA y $TAX_TSV ya existen"
else
  echo "[run] export FASTA + TSV"
  TMPEXP=$(mktemp -d)
  qiime tools export --input-path "$SEQS_QZA" --output-path "$TMPEXP/seqs"
  qiime tools export --input-path "$TAX_QZA"  --output-path "$TMPEXP/tax"

  python3 - <<PY
import csv, os, shutil
tmap = {}
with open(os.path.join("$TMPEXP","tax","taxonomy.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        tmap[row["Feature ID"]] = row["Taxon"]
with open(os.path.join("$TMPEXP","seqs","dna-sequences.fasta")) as fi, \
     open("$SEQS_FA", "w") as fo:
    for line in fi:
        if line.startswith(">"):
            sid = line[1:].strip().split()[0]
            fo.write(f">{tmap.get(sid,'Unassigned')}\n")
        else:
            fo.write(line)
shutil.copy(os.path.join("$TMPEXP","tax","taxonomy.tsv"), "$TAX_TSV")
print(f"FASTA DADA2 : $SEQS_FA")
print(f"TSV taxonomy: $TAX_TSV")
PY
  rm -rf "$TMPEXP"
fi

echo "=========================================================="
echo " PASO 2 — Clasificar ASVs (classify-consensus-vsearch)"
echo "=========================================================="
echo " perc_identity=$PID  maxaccepts=$MAXACC  min_consensus=$MINCONS  threads=$THREADS"

# ---- 2a) Importar ASVs a QIIME2 --------------------------------------------
ASV_QZA="$TAXO_OUT/asvs.qza"
if [[ -s "$ASV_QZA" ]]; then
  echo "[skip] $ASV_QZA ya existe"
else
  echo "[run] import ASVs"
  qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path "$ASV_FASTA" \
    --output-path "$ASV_QZA"
fi

# ---- 2b) classify-consensus-vsearch vs BOLD --------------------------------
CLASSIF_QZA="$TAXO_OUT/tax_BOLD_classification.qza"
SEARCH_QZA="$TAXO_OUT/tax_BOLD_search.qza"
if [[ -s "$CLASSIF_QZA" ]]; then
  echo "[skip] $CLASSIF_QZA ya existe"
else
  echo "[run] classify-consensus-vsearch"
  time qiime feature-classifier classify-consensus-vsearch \
    --i-query          "$ASV_QZA" \
    --i-reference-reads "$SEQS_QZA" \
    --i-reference-taxonomy "$TAX_QZA" \
    --p-perc-identity  "$PID" \
    --p-maxaccepts     "$MAXACC" \
    --p-min-consensus  "$MINCONS" \
    --p-strand both \
    --p-threads        "$THREADS" \
    --o-classification "$CLASSIF_QZA" \
    --o-search-results "$SEARCH_QZA" \
    --verbose 2>&1 | tee "logs/05_classify_vsearch.log"
fi

# ---- 2c) Export TSV para R downstream --------------------------------------
TAX_BOLD_TSV="$TAXO_OUT/tax_BOLD.tsv"
TMPEXP=$(mktemp -d)
qiime tools export --input-path "$CLASSIF_QZA" --output-path "$TMPEXP"
cp "$TMPEXP/taxonomy.tsv" "$TAX_BOLD_TSV"
rm -rf "$TMPEXP"

# ---- Resumen ---------------------------------------------------------------
N_TOTAL=$(($(wc -l < "$TAX_BOLD_TSV") - 1))
N_UNASS=$(grep -c $'\tUnassigned\t' "$TAX_BOLD_TSV" 2>/dev/null || echo 0)
N_ASS=$((N_TOTAL - N_UNASS))

echo
echo "=========================================================="
echo " Clasificación BOLD completa — $(date -Iseconds)"
echo " Total ASVs:      $N_TOTAL"
echo " Asignados:       $N_ASS"
echo " Unassigned:      $N_UNASS"
echo "----------------------------------------------------------"
echo " Salidas:"
echo "  DB BOLD final:  $SEQS_QZA / $TAX_QZA"
echo "  DB FASTA:       $SEQS_FA"
echo "  DB TSV:         $TAX_TSV"
echo "  Clasificación:  $CLASSIF_QZA (.qza) + $TAX_BOLD_TSV"
echo "----------------------------------------------------------"
echo " Acciones recomendadas:"
echo "  - Inspeccionar tax_BOLD.tsv y calcular % asignación por rango:"
echo "      awk -F'\\t' 'NR>1{n=split(\$2,a,\";\"); print n}' $TAX_BOLD_TSV | sort | uniq -c"
echo "  - Rescate NCBI para Unassigned: workflow/04b_build_NCBI_db_RESCRIPt.sh"
echo "  - Consenso híbrido BOLD+NCBI: workflow/04_taxonomy_BOLD_NCBI.R"
echo "=========================================================="
