#!/usr/bin/env bash
# =============================================================================
# 04b_build_NCBI_db_RESCRIPt.sh
# Construcción de una base de datos COI curada desde NCBI usando RESCRIPt
# -----------------------------------------------------------------------------
# Adaptado del tutorial oficial:
#   https://forum.qiime2.org/t/building-a-coi-database-from-ncbi-references/16500
#
# Estrategia:
#   1) Descargar secuencias COI desde NCBI Nucleotide vía Entrez query
#      (config: taxonomy.ncbi.entrez_query)
#   2) Recortar headers / asignar taxonomía completa a 7 ranks
#   3) Filtros: longitud, ambigüedad, homopolímeros
#   4) Dereplicar
#   5) Entrenar Naive Bayes classifier (opcional)
#   6) Exportar a FASTA + TSV para consumo en R/DADA2/VSEARCH
# -----------------------------------------------------------------------------
# Requisitos (mismo entorno que 04a):
#   source activate qiime2-2026       # env del cluster LUSTRE
# Variable de entorno opcional:
#   export NCBI_API_KEY="tu_api_key"      # acelera ~10x el rate limit
# -----------------------------------------------------------------------------
# Tiempo esperado: 1–4 h (limitado por rate limit de Entrez).
# Espacio en disco: ~5–15 GB.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT="$( cd "${SCRIPT_DIR}/.." && pwd )"
PARAMS="${ROOT}/config/params.yml"

command -v qiime >/dev/null 2>&1 || {
  echo "[ERROR] Activa el entorno qiime2-2026." >&2; exit 1; }
command -v yq >/dev/null 2>&1 || { echo "[ERROR] yq requerido." >&2; exit 1; }
qiime rescript --help >/dev/null 2>&1 || {
  echo "[ERROR] Plugin RESCRIPt no instalado." >&2; exit 1; }

# ---- Leer parámetros --------------------------------------------------------
QUERY=$(yq    '.taxonomy.ncbi.entrez_query' "$PARAMS")
MIN_LEN=$(yq  '.taxonomy.ncbi.min_length'   "$PARAMS")
MAX_LEN=$(yq  '.taxonomy.ncbi.max_length'   "$PARAMS")
MAX_AMBIG=$(yq '.taxonomy.ncbi.max_ambig'   "$PARAMS")
DEREP_MODE=$(yq '.taxonomy.ncbi.dereplicate_mode' "$PARAMS")
SEQS_QZA=$(yq '.taxonomy.ncbi.seqs_qza'     "$PARAMS")
TAX_QZA=$(yq  '.taxonomy.ncbi.tax_qza'      "$PARAMS")
CLS_QZA=$(yq  '.taxonomy.ncbi.classifier_qza' "$PARAMS")
SEQS_FA=$(yq  '.taxonomy.ncbi.seqs_fasta'   "$PARAMS")
TAX_TSV=$(yq  '.taxonomy.ncbi.tax_tsv'      "$PARAMS")
CORES=$(yq    '.compute.cores_total'        "$PARAMS")

mapfile -t RANKS < <(yq '.taxonomy.ncbi.ranks[]' "$PARAMS")
RANKS_JOIN=$(IFS=' '; echo "${RANKS[*]}")

# ---- Directorios -----------------------------------------------------------
OUT_DIR="${ROOT}/db/rescript"
WORK_DIR="${ROOT}/db/rescript/_work_ncbi"
LOG_DIR="${ROOT}/logs/04b_ncbi"
mkdir -p "$OUT_DIR" "$WORK_DIR" "$LOG_DIR"
cd "$WORK_DIR"

echo "[04b] Query:    $QUERY"
echo "[04b] Ranks:    $RANKS_JOIN"
echo "[04b] Filtros:  length [$MIN_LEN, $MAX_LEN], max ambig $MAX_AMBIG"
echo "[04b] Output:   ${ROOT}/${SEQS_QZA}"
echo

# ---- 1) Descargar de NCBI vía RESCRIPt --------------------------------------
RAW_SEQ="ncbi-coi-raw-seqs.qza"
RAW_TAX="ncbi-coi-raw-tax.qza"

if [[ -s "$RAW_SEQ" && -s "$RAW_TAX" ]]; then
  echo "[skip] Descarga NCBI ya hecha"
else
  echo "[run]  get-ncbi-data (esto puede tardar 1–4 h)..."
  qiime rescript get-ncbi-data \
    --p-query "$QUERY" \
    --p-ranks $RANKS_JOIN \
    --p-n-jobs "$CORES" \
    --o-sequences "$RAW_SEQ" \
    --o-taxonomy  "$RAW_TAX" \
    --verbose 2>&1 | tee "${LOG_DIR}/get-ncbi-data.log"
fi

# ---- 2) Cull-seqs: descartar seqs con muchas degeneradas/homopolímeros -----
echo "[04b] cull-seqs..."
qiime rescript cull-seqs \
  --i-sequences "$RAW_SEQ" \
  --p-num-degenerates    "$MAX_AMBIG" \
  --p-homopolymer-length 12 \
  --p-n-jobs "$CORES" \
  --o-clean-sequences ncbi-coi-culled.qza

# ---- 3) Filtro de longitud -------------------------------------------------
echo "[04b] filter-seqs-length..."
qiime rescript filter-seqs-length \
  --i-sequences   ncbi-coi-culled.qza \
  --p-global-min  "$MIN_LEN" \
  --p-global-max  "$MAX_LEN" \
  --o-filtered-seqs  ncbi-coi-len-filt.qza \
  --o-discarded-seqs ncbi-coi-len-discarded.qza

qiime rescript filter-taxa \
  --i-taxonomy "$RAW_TAX" \
  --m-ids-to-keep-file ncbi-coi-len-filt.qza \
  --o-filtered-taxonomy ncbi-coi-len-filt-tax.qza

# ---- 4) Dereplicar ---------------------------------------------------------
echo "[04b] dereplicate (mode=$DEREP_MODE)..."
qiime rescript dereplicate \
  --i-sequences ncbi-coi-len-filt.qza \
  --i-taxa      ncbi-coi-len-filt-tax.qza \
  --p-mode      "$DEREP_MODE" \
  --p-threads   "$CORES" \
  --o-dereplicated-sequences "../$(basename $SEQS_QZA)" \
  --o-dereplicated-taxa      "../$(basename $TAX_QZA)"

# ---- 5) Entrenar Naive Bayes (opcional) ------------------------------------
echo "[04b] fit-classifier-naive-bayes..."
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads    "../$(basename $SEQS_QZA)" \
  --i-reference-taxonomy "../$(basename $TAX_QZA)" \
  --o-classifier         "../$(basename $CLS_QZA)" \
  --verbose 2>&1 | tee "${LOG_DIR}/fit-classifier.log"

# ---- 6) Exportar a FASTA + TSV para DADA2/VSEARCH --------------------------
echo "[04b] Exportando a FASTA + TSV..."
TMPEXP=$(mktemp -d)
qiime tools export --input-path "../$(basename $SEQS_QZA)" --output-path "$TMPEXP/seqs"
qiime tools export --input-path "../$(basename $TAX_QZA)" --output-path "$TMPEXP/tax"

python3 - <<PY
import csv, os, shutil
tax_map = {}
with open(os.path.join("$TMPEXP", "tax", "taxonomy.tsv")) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        tax_map[row["Feature ID"]] = row["Taxon"]

in_fa  = os.path.join("$TMPEXP", "seqs", "dna-sequences.fasta")
out_fa = "${ROOT}/${SEQS_FA}"
with open(in_fa) as fi, open(out_fa, "w") as fo:
    for line in fi:
        if line.startswith(">"):
            sid = line[1:].strip().split()[0]
            tax = tax_map.get(sid, "Unassigned")
            fo.write(f">{tax}\n")
        else:
            fo.write(line)
print("FASTA DADA2-format:", out_fa)
shutil.copy(os.path.join("$TMPEXP", "tax", "taxonomy.tsv"), "${ROOT}/${TAX_TSV}")
print("TSV taxonomy:      ${ROOT}/${TAX_TSV}")
PY

rm -rf "$TMPEXP"

echo
echo "=== Resumen 04b (NCBI via RESCRIPt) ==="
echo "QZA seqs:        ${ROOT}/${SEQS_QZA}"
echo "QZA tax:         ${ROOT}/${TAX_QZA}"
echo "QZA classifier:  ${ROOT}/${CLS_QZA}"
echo "FASTA DADA2:     ${ROOT}/${SEQS_FA}"
echo "TSV tax:         ${ROOT}/${TAX_TSV}"
echo
echo "Acciones recomendadas:"
echo " - Verificar: grep -c '^>' ${ROOT}/${SEQS_FA}"
echo " - Continuar con workflow/04_taxonomy_BOLD_NCBI.R (hybrid assignment)."
