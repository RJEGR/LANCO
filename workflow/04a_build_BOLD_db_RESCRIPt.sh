#!/usr/bin/env bash
# =============================================================================
# 04a_build_BOLD_db_RESCRIPt.sh
# Construcción de una base de datos COI curada desde BOLD usando RESCRIPt
# -----------------------------------------------------------------------------
# Adaptado del tutorial oficial:
#   https://forum.qiime2.org/t/building-a-coi-database-from-bold-references/16129
#
# Estrategia:
#   1) Descargar secuencias COI-5P desde BOLD para los taxa de interés
#      (config: taxonomy.bold.taxa  →  Metazoa, Algae, Fungi para LANCO)
#   2) Limpiar headers, normalizar taxonomía, filtrar por longitud y ambigüedad
#   3) Dereplicar para colapsar redundancia
#   4) Entrenar (opcional) un Naive Bayes classifier
#   5) Exportar a FASTA + TSV consumible por R/DADA2 y VSEARCH
# -----------------------------------------------------------------------------
# Requisitos:
#   bash workflow/00_setup_envs.sh rescript   # crea qiime2-rescript + RESCRIPt + yq v4
#   mamba activate qiime2-rescript
# (yq v4 -mikefarah/yq- no existe en conda; 00_setup_envs.sh lo instala como
#  binario directo desde GitHub releases en el bin del entorno activo.)
# -----------------------------------------------------------------------------
# Tiempo esperado: 2–8 h dependiendo del taxon (Metazoa es el más pesado).
# Espacio en disco: ~10–25 GB para BOLD Metazoa completo.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT="$( cd "${SCRIPT_DIR}/.." && pwd )"
PARAMS="${ROOT}/config/params.yml"

# ---- Verificar dependencias -------------------------------------------------
command -v qiime >/dev/null 2>&1 || {
  echo "[ERROR] qiime no encontrado. Activa el entorno qiime2-rescript." >&2; exit 1; }
command -v yq    >/dev/null 2>&1 || { echo "[ERROR] yq requerido." >&2; exit 1; }

qiime rescript --help >/dev/null 2>&1 || {
  echo "[ERROR] Plugin RESCRIPt no instalado. Instala con:" >&2
  echo "  pip install git+https://github.com/bokulich-lab/RESCRIPt.git" >&2
  exit 1; }

# ---- Leer parámetros --------------------------------------------------------
MARKER=$(yq '.taxonomy.bold.rescript_marker' "$PARAMS")
MIN_LEN=$(yq '.taxonomy.bold.min_length'     "$PARAMS")
MAX_LEN=$(yq '.taxonomy.bold.max_length'     "$PARAMS")
MAX_AMBIG=$(yq '.taxonomy.bold.max_ambig'    "$PARAMS")
DEREP_MODE=$(yq '.taxonomy.bold.dereplicate_mode' "$PARAMS")
SEQS_QZA=$(yq '.taxonomy.bold.seqs_qza'      "$PARAMS")
TAX_QZA=$(yq  '.taxonomy.bold.tax_qza'       "$PARAMS")
CLS_QZA=$(yq  '.taxonomy.bold.classifier_qza' "$PARAMS")
SEQS_FA=$(yq  '.taxonomy.bold.seqs_fasta'    "$PARAMS")
TAX_TSV=$(yq  '.taxonomy.bold.tax_tsv'       "$PARAMS")
CORES=$(yq    '.compute.cores_total'         "$PARAMS")

# Taxa a descargar (array)
mapfile -t TAXA < <(yq '.taxonomy.bold.taxa[]' "$PARAMS")

# ---- Preparar directorios ---------------------------------------------------
OUT_DIR="${ROOT}/db/rescript"
WORK_DIR="${ROOT}/db/rescript/_work_bold"
LOG_DIR="${ROOT}/logs/04a_bold"
mkdir -p "$OUT_DIR" "$WORK_DIR" "$LOG_DIR"

cd "$WORK_DIR"

echo "[04a] Marker:    $MARKER"
echo "[04a] Taxa:      ${TAXA[*]}"
echo "[04a] Filtros:   length [$MIN_LEN, $MAX_LEN], max ambig $MAX_AMBIG"
echo "[04a] Output:    ${ROOT}/${SEQS_QZA}"
echo

# ---- 1) Descargar BOLD por cada taxón --------------------------------------
DOWNLOADED_SEQS=()
DOWNLOADED_TAX=()
for TAXON in "${TAXA[@]}"; do
  RAW_SEQ="bold-${TAXON,,}-raw-seqs.qza"
  RAW_TAX="bold-${TAXON,,}-raw-tax.qza"

  if [[ -s "$RAW_SEQ" && -s "$RAW_TAX" ]]; then
    echo "[skip] $TAXON ya descargado"
  else
    echo "[run]  Descargando BOLD: $TAXON ($MARKER) ..."
    qiime rescript get-bold-data \
      --p-ranks   phylum class order family subfamily genus species \
      --p-taxa    "$TAXON" \
      --p-marker  "$MARKER" \
      --o-sequences "$RAW_SEQ" \
      --o-taxonomy  "$RAW_TAX" \
      --verbose 2>&1 | tee "${LOG_DIR}/get-bold-${TAXON}.log"
  fi
  DOWNLOADED_SEQS+=("$RAW_SEQ")
  DOWNLOADED_TAX+=("$RAW_TAX")
done

# ---- 2) Mergear los taxa descargados ---------------------------------------
if [[ ${#DOWNLOADED_SEQS[@]} -gt 1 ]]; then
  echo "[04a] Mergeando ${#DOWNLOADED_SEQS[@]} taxa..."
  qiime feature-table merge-seqs \
    --i-data "${DOWNLOADED_SEQS[@]}" \
    --o-merged-data bold-merged-seqs.qza
  qiime feature-table merge-taxa \
    --i-data "${DOWNLOADED_TAX[@]}" \
    --o-merged-data bold-merged-tax.qza
else
  cp "${DOWNLOADED_SEQS[0]}" bold-merged-seqs.qza
  cp "${DOWNLOADED_TAX[0]}"  bold-merged-tax.qza
fi

# ---- 3) Limpiar headers y normalizar taxonomía -----------------------------
echo "[04a] edit-taxonomy: normalizando ranks..."
qiime rescript edit-taxonomy \
  --i-taxonomy bold-merged-tax.qza \
  --p-search-strings  ' ' ',' ';' \
  --p-replacement-strings '_' '' '_' \
  --o-edited-taxonomy  bold-edited-tax.qza

# ---- 4) Filtro por longitud y ambigüedad -----------------------------------
echo "[04a] cull-seqs: descartar seqs degeneradas..."
qiime rescript cull-seqs \
  --i-sequences  bold-merged-seqs.qza \
  --p-num-degenerates "$MAX_AMBIG" \
  --p-homopolymer-length 12 \
  --p-n-jobs "$CORES" \
  --o-clean-sequences bold-culled-seqs.qza

echo "[04a] filter-seqs-length-by-taxon (250–1600 nt para COI)..."
qiime rescript filter-seqs-length \
  --i-sequences   bold-culled-seqs.qza \
  --p-global-min  "$MIN_LEN" \
  --p-global-max  "$MAX_LEN" \
  --o-filtered-seqs    bold-len-filt-seqs.qza \
  --o-discarded-seqs   bold-len-discarded.qza

# Realinear taxonomía al subset filtrado
qiime rescript filter-taxa \
  --i-taxonomy bold-edited-tax.qza \
  --m-ids-to-keep-file bold-len-filt-seqs.qza \
  --o-filtered-taxonomy bold-len-filt-tax.qza

# ---- 5) Dereplicar ---------------------------------------------------------
echo "[04a] dereplicate (mode=$DEREP_MODE)..."
qiime rescript dereplicate \
  --i-sequences bold-len-filt-seqs.qza \
  --i-taxa      bold-len-filt-tax.qza \
  --p-mode      "$DEREP_MODE" \
  --p-threads   "$CORES" \
  --o-dereplicated-sequences "../$(basename $SEQS_QZA)" \
  --o-dereplicated-taxa      "../$(basename $TAX_QZA)"

# ---- 6) Entrenar Naive Bayes (opcional pero recomendado) ------------------
echo "[04a] fit-classifier-naive-bayes (~30–90 min, depende del tamaño)..."
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads    "../$(basename $SEQS_QZA)" \
  --i-reference-taxonomy "../$(basename $TAX_QZA)" \
  --o-classifier         "../$(basename $CLS_QZA)" \
  --verbose 2>&1 | tee "${LOG_DIR}/fit-classifier.log"

# ---- 7) Exportar a FASTA + TSV para DADA2/VSEARCH (en R) -------------------
echo "[04a] Exportando a FASTA + TSV..."
TMPEXP=$(mktemp -d)
qiime tools export --input-path "../$(basename $SEQS_QZA)" --output-path "$TMPEXP/seqs"
qiime tools export --input-path "../$(basename $TAX_QZA)" --output-path "$TMPEXP/tax"

# Headers FASTA en formato DADA2: "Kingdom;Phylum;Class;Order;Family;Genus;Species"
python3 - <<PY
import csv, os
tax_map = {}
with open(os.path.join("$TMPEXP", "tax", "taxonomy.tsv")) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        tax_map[row["Feature ID"]] = row["Taxon"]

in_fa  = os.path.join("$TMPEXP", "seqs", "dna-sequences.fasta")
out_fa = "${ROOT}/${SEQS_FA}"
with open(in_fa) as fi, open(out_fa, "w") as fo:
    sid = None
    for line in fi:
        if line.startswith(">"):
            sid = line[1:].strip().split()[0]
            tax = tax_map.get(sid, "Unassigned")
            fo.write(f">{tax}\n")
        else:
            fo.write(line)
print("FASTA DADA2-format:", out_fa)

# Tabla TSV plana
import shutil
shutil.copy(os.path.join("$TMPEXP", "tax", "taxonomy.tsv"), "${ROOT}/${TAX_TSV}")
print("TSV taxonomy:      ${ROOT}/${TAX_TSV}")
PY

rm -rf "$TMPEXP"

echo
echo "=== Resumen 04a (BOLD via RESCRIPt) ==="
echo "QZA seqs:        ${ROOT}/${SEQS_QZA}"
echo "QZA tax:         ${ROOT}/${TAX_QZA}"
echo "QZA classifier:  ${ROOT}/${CLS_QZA}"
echo "FASTA DADA2:     ${ROOT}/${SEQS_FA}"
echo "TSV tax:         ${ROOT}/${TAX_TSV}"
echo
echo "Acciones recomendadas:"
echo " - Verificar tamaño y nº de secuencias: grep -c '^>' ${ROOT}/${SEQS_FA}"
echo " - Construir DB NCBI: workflow/04b_build_NCBI_db_RESCRIPt.sh"
echo " - Asignar taxonomía:  workflow/04_taxonomy_BOLD_NCBI.R"
