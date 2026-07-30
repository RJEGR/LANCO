#!/bin/bash
# =============================================================================
# 06a_classify_OTUs_BOLD.slurm  — VSEARCH usearch_global vs BOLD COI DB
# LANCO/RUN29 — clasifica OTUs curados (paso 05) contra BOLD.
# Salida: OTUs asignados + FASTA de unassigned para paso 06b (NCBI).
# =============================================================================
#SBATCH --job-name=06a_bold
#SBATCH --output=logs/06a_classify_OTUs_BOLD_%j.out
#SBATCH --error=logs/06a_classify_OTUs_BOLD_%j.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --time=04:00:00
# --- AJUSTA a tu cluster: partition, account, mail ---
# #SBATCH --partition=<cluster_partition>
# #SBATCH --account=<cluster_account>
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=rgomez41@uabc.edu.mx
# =============================================================================

set -euo pipefail

# ---- Entorno --------------------------------------------------------------
# Ajustar path del conda init según cluster (típicamente /lustre/apps/miniconda/...)
#source ~/miniconda3/etc/profile.d/conda.sh
#conda activate coi_pipeline    # env con swarm+vsearch (docs/INSTALL_dependencies.md)

# ---- Env --------------------------------------------------------------------
export PATH=$HOME/.local/bin:/LUSTRE/apps/Anaconda/2026/miniconda3/condabin:/LUSTRE/apps/Anaconda/2026/miniconda3/bin:$PATH
# shellcheck disable=SC1091
source activate qiime2-2026

# ---- Paths (ROOT = raíz del repo LANCO en LUSTRE) -------------------------
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$ROOT"

OTUS_FA="results/05_swarm_lulu/lulu/curated_otus.fasta"
BOLD_DB="db/rescript/bold-coi-seqs.fasta"
OUT_DIR="results/06_taxonomy/bold"
BLAST6_OUT="${OUT_DIR}/otus_vs_bold.b6"

mkdir -p "$OUT_DIR" logs

# ---- Preflight ------------------------------------------------------------
for f in "$OTUS_FA" "$BOLD_DB"; do
  [[ -s "$f" ]] || { echo "ERROR: falta $f"; exit 1; }
done
command -v vsearch >/dev/null || { echo "ERROR: vsearch no en PATH"; exit 1; }

echo "[06a] host=$(hostname) date=$(date -Iseconds)"
echo "[06a] OTUs input:   $OTUS_FA  ($(grep -c '^>' "$OTUS_FA") seqs)"
echo "[06a] BOLD DB:      $BOLD_DB  ($(grep -c '^>' "$BOLD_DB") seqs)"
echo "[06a] threads=${SLURM_CPUS_PER_TASK:-24}  RAM=${SLURM_MEM_PER_NODE:-100G}"

# ---- VSEARCH usearch_global (parámetros de config/params.yml) -------------
vsearch \
  --usearch_global   "$OTUS_FA" \
  --db               "$BOLD_DB" \
  --id               0.85 \
  --maxaccepts       10 \
  --maxrejects       32 \
  --top_hits_only \
  --output_no_hits \
  --query_cov        0.9 \
  --blast6out        "$BLAST6_OUT" \
  --threads          "${SLURM_CPUS_PER_TASK:-24}" \
  --strand           both \
  --log              "${OUT_DIR}/vsearch_bold.log"

echo "[06a] VSEARCH terminado. hits=$(wc -l < "$BLAST6_OUT") líneas."

# ---- Consensus taxonomy en R ----------------------------------------------
Rscript workflow/06a_classify_OTUs_BOLD.R

echo "[06a] Done. Siguiente: sbatch workflow/06b_classify_OTUs_NCBI.slurm"
