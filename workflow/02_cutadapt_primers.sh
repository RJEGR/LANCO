#!/usr/bin/env bash
# =============================================================================
# 02_cutadapt_primers.sh
# Remoción de primers Leray-XT (mlCOIintF / jgHCO2198) con cutadapt
# -----------------------------------------------------------------------------
# Estrategia: "linked adapters" en ambas direcciones, con
# --discard-untrimmed (descarta reads donde no se haya encontrado el primer
# en ambos extremos esperados), porque en COI cualquier read sin primer es
# señal de off-target/lecturas mezcladas y degrada DADA2.
# -----------------------------------------------------------------------------
# Requisitos: cutadapt >= 4.0  (conda install -c bioconda cutadapt)
# Entrada:   RUN29/*_R{1,2}_001.fastq.gz
# Salida:    results/02_cutadapt/<sample>_R{1,2}.fastq.gz + log por muestra
# =============================================================================

set -euo pipefail

# ---- Localizar raíz del repo (donde vive RUN29/) ----------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT="$( cd "${SCRIPT_DIR}/.." && pwd )"

# ---- Parámetros (leemos params.yml con yq si está disponible, si no fallback)
PARAMS="${ROOT}/config/params.yml"
if command -v yq >/dev/null 2>&1; then
  FWD=$(yq '.primers.fwd_seq'        "$PARAMS")
  REV=$(yq '.primers.rev_seq'        "$PARAMS")
  MIN_LEN=$(yq '.cutadapt.min_length' "$PARAMS")
  ERR=$(yq '.cutadapt.error_rate'    "$PARAMS")
  CORES=$(yq '.cutadapt.cores'       "$PARAMS")
else
  echo "[WARN] yq no encontrado — uso valores por defecto" >&2
  FWD="GGWACWGGWTGAACWGTWTAYCCYCC"
  REV="TANACYTCNGGRTGNCCRAARAAYCA"
  MIN_LEN=50
  ERR=0.15
  CORES=4
fi

# ---- Complemento reverso de los primers (cutadapt necesita las 4 secuencias)-
# Usamos python (siempre disponible vía cutadapt) para revcomp seguro
revcomp () {
  python3 -c "
from Bio.Seq import Seq
import sys
print(str(Seq(sys.argv[1]).reverse_complement()))
" "$1" 2>/dev/null || python3 -c "
tbl = str.maketrans('ACGTRYSWKMBDHVN','TGCAYRWSMKVHDBN')
import sys; s = sys.argv[1].upper()[::-1]
print(s.translate(tbl))
" "$1"
}

FWD_RC=$(revcomp "$FWD")
REV_RC=$(revcomp "$REV")

echo "[02] Primers:"
echo "    FWD     = $FWD"
echo "    REV     = $REV"
echo "    FWD_RC  = $FWD_RC"
echo "    REV_RC  = $REV_RC"

# ---- Directorios ------------------------------------------------------------
IN_DIR="${ROOT}/RUN29"
OUT_DIR="${ROOT}/results/02_cutadapt"
LOG_DIR="${ROOT}/logs/02_cutadapt"
mkdir -p "$OUT_DIR" "$LOG_DIR"

# ---- Loop por muestra -------------------------------------------------------
shopt -s nullglob
for R1 in "${IN_DIR}"/*_R1_001.fastq.gz; do
  BN=$(basename "$R1")
  SAMPLE=$(echo "$BN" | sed -E 's/_S[0-9]+_L001_R1_001\.fastq\.gz$//')
  R2="${R1/_R1_/_R2_}"

  OUT_R1="${OUT_DIR}/${SAMPLE}_R1.fastq.gz"
  OUT_R2="${OUT_DIR}/${SAMPLE}_R2.fastq.gz"
  LOG="${LOG_DIR}/${SAMPLE}.log"

  if [[ -s "$OUT_R1" && -s "$OUT_R2" ]]; then
    echo "[skip] $SAMPLE ya procesado"
    continue
  fi

  echo "[run]  $SAMPLE"
  cutadapt \
    -g  "$FWD" \
    -G  "$REV" \
    -a  "$REV_RC" \
    -A  "$FWD_RC" \
    --discard-untrimmed \
    --minimum-length "$MIN_LEN" \
    -e  "$ERR" \
    -j  "$CORES" \
    -o  "$OUT_R1" \
    -p  "$OUT_R2" \
    "$R1" "$R2" \
    > "$LOG" 2>&1
done

# ---- Resumen ----------------------------------------------------------------
echo
echo "=== Resumen 02_cutadapt ==="
SUMMARY="${OUT_DIR}/cutadapt_summary.tsv"
{
  printf "sample_id\treads_in\treads_passed\tpct_passed\n"
  for L in "${LOG_DIR}"/*.log; do
    SAMPLE=$(basename "$L" .log)
    IN=$(grep -E "Total read pairs processed" "$L" | awk -F: '{gsub(/[^0-9]/,"",$2); print $2}')
    OK=$(grep -E "Pairs written \(passing filters\)" "$L" | awk -F: '{print $2}' | awk '{gsub(",","",$1); print $1}')
    PCT=$(awk -v a="$OK" -v b="$IN" 'BEGIN{if(b>0) printf "%.2f", 100*a/b; else print "NA"}')
    printf "%s\t%s\t%s\t%s\n" "$SAMPLE" "$IN" "$OK" "$PCT"
  done
} > "$SUMMARY"

echo "Summary: $SUMMARY"
echo
echo "Acciones recomendadas:"
echo " - Revisar $SUMMARY: pct_passed > 80 % es lo esperable para Leray-XT."
echo " - Si una muestra tiene pct_passed < 50 %, revisar log individual."
echo " - Continuar con workflow/03_dada2_pipeline.R."
