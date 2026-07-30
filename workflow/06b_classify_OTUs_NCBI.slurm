#!/bin/bash
# =============================================================================
# 06b_classify_OTUs_NCBI.slurm  — VSEARCH vs NCBI COI, solo OTUs huérfanos BOLD
# LANCO/RUN29 — segunda fase de clasificación taxonómica.
# Input:  results/06_taxonomy/bold/otus_unassigned_BOLD.fasta
# Output: results/06_taxonomy/tax_consensus.tsv  (BOLD + NCBI + unassigned)
# =============================================================================
#SBATCH --job-name=06b_ncbi
#SBATCH --output=logs/06b_classify_OTUs_NCBI_%j.out
#SBATCH --error=logs/06b_classify_OTUs_NCBI_%j.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --time=03:00:00
# --- AJUSTA a tu cluster: partition, account, mail ---
# #SBATCH --partition=<cluster_partition>
# #SBATCH --account=<cluster_account>
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=rgomez41@uabc.edu.mx
# #SBATCH --dependency=afterok:<06a_jobid>   # opcional: encadena tras 06a
# =============================================================================

set -euo pipefail

# ---- Entorno --------------------------------------------------------------
#source ~/miniconda3/etc/profile.d/conda.sh
#conda activate coi_pipeline

# ---- Env --------------------------------------------------------------------
export PATH=$HOME/.local/bin:/LUSTRE/apps/Anaconda/2026/miniconda3/condabin:/LUSTRE/apps/Anaconda/2026/miniconda3/bin:$PATH
# shellcheck disable=SC1091
source activate qiime2-2026

# ---- Paths ---------------------------------------------------------------
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$ROOT"

IN_FA="results/06_taxonomy/bold/otus_unassigned_BOLD.fasta"
NCBI_DB="db/rescript/ncbi-coi-seqs.fasta"
OUT_DIR="results/06_taxonomy/ncbi"
BLAST6_OUT="${OUT_DIR}/otus_vs_ncbi.b6"

mkdir -p "$OUT_DIR" logs

# ---- Preflight ------------------------------------------------------------
for f in "$IN_FA" "$NCBI_DB"; do
  [[ -s "$f" ]] || { echo "ERROR: falta $f"; exit 1; }
done
command -v vsearch >/dev/null || { echo "ERROR: vsearch no en PATH"; exit 1; }

# Edge case: si BOLD asignó todo, no hay nada que hacer
N_IN=$(grep -c '^>' "$IN_FA" || true)
if [[ "$N_IN" -eq 0 ]]; then
  echo "[06b] Sin OTUs unassigned tras BOLD. Nada que hacer."
  # crear blast6 vacío para que el R merge funcione
  touch "$BLAST6_OUT"
  Rscript workflow/06b_classify_OTUs_NCBI.R
  exit 0
fi

echo "[06b] host=$(hostname) date=$(date -Iseconds)"
echo "[06b] OTUs unassigned BOLD: $N_IN"
echo "[06b] NCBI DB:  $NCBI_DB ($(grep -c '^>' "$NCBI_DB") seqs)"

# ---- VSEARCH usearch_global (mismos parámetros que 06a) -------------------
vsearch \
  --usearch_global   "$IN_FA" \
  --db               "$NCBI_DB" \
  --id               0.85 \
  --maxaccepts       10 \
  --maxrejects       32 \
  --top_hits_only \
  --output_no_hits \
  --query_cov        0.9 \
  --blast6out        "$BLAST6_OUT" \
  --threads          "${SLURM_CPUS_PER_TASK:-24}" \
  --strand           both \
  --log              "${OUT_DIR}/vsearch_ncbi.log"

echo "[06b] VSEARCH terminado. hits=$(wc -l < "$BLAST6_OUT") líneas."

# ---- Consensus + merge final ---------------------------------------------
Rscript workflow/06b_classify_OTUs_NCBI.R

echo "[06b] Pipeline taxonomía completo."
echo "[06b] Output final: results/06_taxonomy/tax_consensus.tsv"
