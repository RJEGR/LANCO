#!/usr/bin/env bash
# =============================================================================
# 04a_build_BOLD_db_SciLifeLab.sh  —  ALTERNATIVA a RESCRIPt get-bold-data
# -----------------------------------------------------------------------------
# Descarga el dataset "COI reference sequences from BOLD DB" (SciLifeLab / NBIS,
# Sundh 2026-06-08 v6, DOI 10.17044/scilifelab.20514192), derivado del
# BOLD Data Package oficial del 15-May-2026. Los archivos ya vienen curados,
# clusterizados al 100 % por BOLD BIN y con taxonomía en formato QIIME2.
#
# Ventaja vs. RESCRIPt get-bold-data (que YA NO EXISTE en releases recientes):
#   - Fresco (BDP 15-May-2026, ~24 M secuencias BOLD)
#   - Un único download (2 archivos, ~416 MB)
#   - Sin rate-limits de la API pública de BOLD
#   - Reproducible: mismo DOI, checksums MD5/SHA en shasum.txt
#
# Trade-off:
#   - Alcance TOTAL de BOLD (no solo Metazoa+Algae+Fungi como en params.yml).
#     Si necesitas restringir a esos taxa: filter-taxa post-import (opcional).
# -----------------------------------------------------------------------------
# Requisitos:
#   source activate qiime2-2026
#   yq >= 4.18   (para leer params.yml)
#   ~5 GB libres en disco
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT="$( cd "${SCRIPT_DIR}/.." && pwd )"
PARAMS="${ROOT}/config/params.yml"

# ---- Verificaciones --------------------------------------------------------
command -v qiime >/dev/null || { echo "[ERROR] activa qiime2-2026" >&2; exit 1; }
command -v yq    >/dev/null || { echo "[ERROR] yq requerido"       >&2; exit 1; }
command -v curl  >/dev/null || { echo "[ERROR] curl requerido"     >&2; exit 1; }
qiime rescript --help >/dev/null 2>&1 || {
  echo "[WARN] RESCRIPt no encontrado — los filtros cull/filter/dereplicate se saltarán." >&2
  HAS_RESCRIPT=0
}
HAS_RESCRIPT=${HAS_RESCRIPT:-1}

# ---- Parámetros ------------------------------------------------------------
MIN_LEN=$(yq   '.taxonomy.bold.min_length'       "$PARAMS")
MAX_LEN=$(yq   '.taxonomy.bold.max_length'       "$PARAMS")
MAX_AMBIG=$(yq '.taxonomy.bold.max_ambig'        "$PARAMS")
DEREP_MODE=$(yq '.taxonomy.bold.dereplicate_mode' "$PARAMS")
SEQS_QZA=$(yq  '.taxonomy.bold.seqs_qza'         "$PARAMS")
TAX_QZA=$(yq   '.taxonomy.bold.tax_qza'          "$PARAMS")
CLS_QZA=$(yq   '.taxonomy.bold.classifier_qza'   "$PARAMS")
SEQS_FA=$(yq   '.taxonomy.bold.seqs_fasta'       "$PARAMS")
TAX_TSV=$(yq   '.taxonomy.bold.tax_tsv'          "$PARAMS")
CORES=$(yq     '.compute.cores_total'            "$PARAMS")

# ---- Config SciLifeLab -----------------------------------------------------
# figshare API IDs verificados vía https://api.figshare.com/v2/articles/20514192
# Versión v6 (2026-06-08), consenso 'exclNA' (mayor resolución de especies).
FASTA_URL="https://ndownloader.figshare.com/files/65310486"   # coidb.clustered.fasta.gz (~360 MB)
FASTA_MD5="7400215d82d34e1a4e6cdf9a60ee1d0d"
TAX_URL="https://ndownloader.figshare.com/files/65310780"     # coidb.qiime2.info.exclNA.tsv.gz (~56 MB)
TAX_MD5="ad4c11d1cd8e7af3e981cbe22b76817a"

# ---- Directorios -----------------------------------------------------------
OUT_DIR="${ROOT}/db/rescript"
WORK_DIR="${ROOT}/db/rescript/_work_bold_scilifelab"
LOG_DIR="${ROOT}/logs/04a_bold"
mkdir -p "$OUT_DIR" "$WORK_DIR" "$LOG_DIR"
cd "$WORK_DIR"

echo "=========================================================="
echo " Fuente:   SciLifeLab COI reference (v6, BDP 15-May-2026)"
echo " DOI:      10.17044/scilifelab.20514192.v6"
echo " Consenso: exclNA (mayor resolución especies)"
echo " Output:   ${ROOT}/${SEQS_QZA}  +  ${ROOT}/${TAX_QZA}"
echo "=========================================================="

# ---- 1) Descargar (idempotente + verifica MD5) -----------------------------
download_check () {
  local url="$1" out="$2" md5="$3"
  if [[ -s "$out" ]] && echo "$md5  $out" | md5sum -c --status 2>/dev/null; then
    echo "[skip] $out ya descargado y MD5 OK"
    return
  fi
  echo "[dl]   $out  ($url)"
  curl -fSL --retry 5 --retry-delay 15 --continue-at - "$url" -o "$out"
  echo "$md5  $out" | md5sum -c || { echo "[FATAL] MD5 mismatch: $out" >&2; exit 1; }
}
download_check "$FASTA_URL" "coidb.clustered.fasta.gz"          "$FASTA_MD5"
download_check "$TAX_URL"   "coidb.qiime2.info.exclNA.tsv.gz"   "$TAX_MD5"

# ---- 2) Descomprimir (mantiene .gz — compatible con gzip 1.5 sin -k) ------
# gzip < 1.6 no soporta -k; usamos redirección con zcat que es portable.
[[ -s coidb.clustered.fasta ]]        || zcat coidb.clustered.fasta.gz        > coidb.clustered.fasta
[[ -s coidb.qiime2.info.exclNA.tsv ]] || zcat coidb.qiime2.info.exclNA.tsv.gz > coidb.qiime2.info.exclNA.tsv

# ---- 3) Normalizar headers FASTA (solo processid, sin ' bin_uri:...') ------
# QIIME2 tools import necesita que los IDs FASTA coincidan con Feature ID del TSV.
# Formato original SciLifeLab: ">{processid} bin_uri:{BIN}"
# → nos quedamos solo con el primer campo (processid).
echo "[fix]  Normalizando headers FASTA → solo processid"
awk 'BEGIN{FS=" "} /^>/ {print ">"substr($1,2); next} {print}' \
    coidb.clustered.fasta > coidb.clustered.norm.fasta

# ---- 4) Importar a artifacts QIIME2 ----------------------------------------
echo "[qiime] Importando FASTA → FeatureData[Sequence]"
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path coidb.clustered.norm.fasta \
  --output-path bold-raw-seqs.qza

echo "[qiime] Importando taxonomía → FeatureData[Taxonomy]"
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path <(tail -n +2 coidb.qiime2.info.exclNA.tsv) \
  --output-path bold-raw-tax.qza \
  2>&1 | tee "${LOG_DIR}/import-tax.log" || {
    # Fallback: algunos releases traen header — reintentar con TSVTaxonomyFormat
    echo "[fallback] Reintentando con TSVTaxonomyFormat (con header)"
    qiime tools import \
      --type 'FeatureData[Taxonomy]' \
      --input-format TSVTaxonomyFormat \
      --input-path coidb.qiime2.info.exclNA.tsv \
      --output-path bold-raw-tax.qza
  }

# ---- 5) Filtros RESCRIPt (opcionales — la DB ya viene curada) --------------
if [[ "$HAS_RESCRIPT" -eq 1 ]]; then
  echo "[04a] cull-seqs (degeneradas / homopolímeros)"
  qiime rescript cull-seqs \
    --i-sequences bold-raw-seqs.qza \
    --p-num-degenerates    "$MAX_AMBIG" \
    --p-homopolymer-length 12 \
    --p-n-jobs "$CORES" \
    --o-clean-sequences bold-culled-seqs.qza

  echo "[04a] filter-seqs-length ($MIN_LEN–$MAX_LEN nt)"
  qiime rescript filter-seqs-length \
    --i-sequences   bold-culled-seqs.qza \
    --p-global-min  "$MIN_LEN" \
    --p-global-max  "$MAX_LEN" \
    --o-filtered-seqs  bold-len-filt-seqs.qza \
    --o-discarded-seqs bold-len-discarded.qza

  qiime rescript filter-taxa \
    --i-taxonomy bold-raw-tax.qza \
    --m-ids-to-keep-file bold-len-filt-seqs.qza \
    --o-filtered-taxonomy bold-len-filt-tax.qza

  echo "[04a] dereplicate (mode=$DEREP_MODE)"
  qiime rescript dereplicate \
    --i-sequences bold-len-filt-seqs.qza \
    --i-taxa      bold-len-filt-tax.qza \
    --p-mode      "$DEREP_MODE" \
    --p-threads   "$CORES" \
    --o-dereplicated-sequences "${ROOT}/${SEQS_QZA}" \
    --o-dereplicated-taxa      "${ROOT}/${TAX_QZA}"
else
  echo "[04a] RESCRIPt ausente — copiando raw directamente a destino"
  cp bold-raw-seqs.qza "${ROOT}/${SEQS_QZA}"
  cp bold-raw-tax.qza  "${ROOT}/${TAX_QZA}"
fi

# ---- 6) Naive Bayes (opcional pero recomendado, 30–90 min) -----------------
echo "[04a] fit-classifier-naive-bayes"
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads    "${ROOT}/${SEQS_QZA}" \
  --i-reference-taxonomy "${ROOT}/${TAX_QZA}" \
  --o-classifier         "${ROOT}/${CLS_QZA}" \
  --verbose 2>&1 | tee "${LOG_DIR}/fit-classifier.log"

# ---- 7) Exportar FASTA (DADA2 assignTaxonomy) + TSV -----------------------
echo "[04a] Exportando a FASTA + TSV"
TMPEXP=$(mktemp -d)
qiime tools export --input-path "${ROOT}/${SEQS_QZA}" --output-path "$TMPEXP/seqs"
qiime tools export --input-path "${ROOT}/${TAX_QZA}"  --output-path "$TMPEXP/tax"

python3 - <<PY
import csv, os, shutil
tmap = {}
with open(os.path.join("$TMPEXP","tax","taxonomy.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        tmap[row["Feature ID"]] = row["Taxon"]

with open(os.path.join("$TMPEXP","seqs","dna-sequences.fasta")) as fi, \
     open("${ROOT}/${SEQS_FA}", "w") as fo:
    for line in fi:
        if line.startswith(">"):
            sid = line[1:].strip().split()[0]
            fo.write(f">{tmap.get(sid,'Unassigned')}\n")
        else:
            fo.write(line)
shutil.copy(os.path.join("$TMPEXP","tax","taxonomy.tsv"), "${ROOT}/${TAX_TSV}")
print(f"FASTA DADA2 : ${ROOT}/${SEQS_FA}")
print(f"TSV taxonomy: ${ROOT}/${TAX_TSV}")
PY
rm -rf "$TMPEXP"

# ---- Resumen ---------------------------------------------------------------
echo
echo "=========================================================="
echo " 04a (BOLD via SciLifeLab) COMPLETO"
echo " QZA seqs:        ${ROOT}/${SEQS_QZA}"
echo " QZA tax:         ${ROOT}/${TAX_QZA}"
echo " QZA classifier:  ${ROOT}/${CLS_QZA}"
echo " FASTA DADA2:     ${ROOT}/${SEQS_FA}"
echo " TSV tax:         ${ROOT}/${TAX_TSV}"
echo "=========================================================="
echo
echo "Acciones recomendadas:"
echo "  - Verificar tamaño: grep -c '^>' ${ROOT}/${SEQS_FA}"
echo "  - Construir DB NCBI (rescate): workflow/04b_build_NCBI_db_RESCRIPt.sh"
echo "  - Asignar taxonomía:            workflow/04_taxonomy_BOLD_NCBI.R"
