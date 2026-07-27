#!/usr/bin/env bash
# =============================================================================
# 05_classify_ASVs_NCBI.slurm
# Clasifica ASVs vs NCBI COI (classify-consensus-vsearch de QIIME2)
# -----------------------------------------------------------------------------
# Dos modos de operación (env var MODE):
#   MODE=rescate  (DEFAULT): reclasifica solo los ASVs marcados como
#     "Unassigned" en results/04_taxonomy/tax_BOLD.tsv. Produce
#     tax_NCBI_rescate.tsv + tax_consensus.tsv (BOLD primary + NCBI rescate).
#   MODE=all: clasifica TODOS los ASVs vs NCBI independientemente
#     (útil para comparar BOLD vs NCBI en el workshop). Produce tax_NCBI.tsv.
# -----------------------------------------------------------------------------
# Prerrequisitos:
#   - db/rescript/ncbi-coi-{seqs,tax}.qza  (generados por 04b)
#   - results/04_taxonomy/tax_BOLD.tsv     (generado por 05_classify_ASVs_BOLD.slurm)
#   - results/03_dada2_SE/asvs.fasta        (ASVs originales de DADA2)
# -----------------------------------------------------------------------------
# Uso:
#   cd /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/LANCO
#   # Rescate (default): solo Unassigned de BOLD
#   sbatch workflow/05_classify_ASVs_NCBI.slurm
#   # Todos los ASVs vs NCBI (para comparación):
#   sbatch --export=ALL,MODE=all workflow/05_classify_ASVs_NCBI.slurm
#   # Override del FASTA de ASVs:
#   sbatch --export=ALL,ASV_FASTA=results/03_dada2/asvs.fasta \
#          workflow/05_classify_ASVs_NCBI.slurm
# =============================================================================
#SBATCH --job-name=lanco_classify_ncbi
#SBATCH --output=logs/05_classify_%x_%j.out
#SBATCH --error=logs/05_classify_%x_%j.err
#SBATCH --time=06:00:00
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
echo " JobID: ${SLURM_JOB_ID:-N/A}"
echo " qiime: $(command -v qiime)"
echo " yq:    $(command -v yq)"

# ---- Parámetros -------------------------------------------------------------
PARAMS="config/params.yml"
PID=$(yq     '.taxonomy.classification.perc_identity'  "$PARAMS")
MAXACC=$(yq  '.taxonomy.classification.maxaccepts'     "$PARAMS")
MINCONS=$(yq '.taxonomy.classification.min_consensus'  "$PARAMS")
THREADS=$(yq '.taxonomy.classification.threads'        "$PARAMS")

NCBI_SEQS_QZA=$(yq '.taxonomy.ncbi.seqs_qza' "$PARAMS")
NCBI_TAX_QZA=$(yq  '.taxonomy.ncbi.tax_qza'  "$PARAMS")

# Modo de operación
MODE="${MODE:-rescate}"
ASV_FASTA="${ASV_FASTA:-results/03_dada2_SE/asvs.fasta}"
TAX_BOLD="results/04_taxonomy/tax_BOLD.tsv"
TAXO_OUT="results/04_taxonomy"

# ---- Sanity checks ----------------------------------------------------------
[[ -s "$NCBI_SEQS_QZA" ]] || { echo "[FATAL] no existe $NCBI_SEQS_QZA (corre 04b antes)"; exit 1; }
[[ -s "$NCBI_TAX_QZA"  ]] || { echo "[FATAL] no existe $NCBI_TAX_QZA (corre 04b antes)";  exit 1; }
[[ -s "$ASV_FASTA"     ]] || { echo "[FATAL] no existe $ASV_FASTA"; exit 1; }

case "$MODE" in
  rescate)
    [[ -s "$TAX_BOLD" ]] || { echo "[FATAL] MODE=rescate requiere $TAX_BOLD (corre 05_classify_ASVs_BOLD antes)"; exit 1; }
    ;;
  all)
    ;;
  *)
    echo "[FATAL] MODE debe ser 'rescate' o 'all' (recibido: $MODE)"; exit 1
    ;;
esac

echo "=========================================================="
echo " MODE:        $MODE"
echo " ASV FASTA:   $ASV_FASTA ($(grep -c '^>' "$ASV_FASTA") ASVs totales)"
echo " NCBI DB:     $NCBI_SEQS_QZA + $NCBI_TAX_QZA"
echo " VSEARCH:     pid=$PID  maxaccepts=$MAXACC  min_cons=$MINCONS  threads=$THREADS"
echo "=========================================================="

# ---- 1) Preparar FASTA de query --------------------------------------------
if [[ "$MODE" == "rescate" ]]; then
  # Extrae Feature IDs de ASVs marcados Unassigned en tax_BOLD.tsv
  QUERY_FA="$TAXO_OUT/asvs_unassigned_BOLD.fasta"
  UNASS_IDS="$TAXO_OUT/asvs_unassigned_BOLD.ids"

  echo "[rescate] Extrayendo Feature IDs Unassigned de $TAX_BOLD"
  # tax_BOLD.tsv formato: Feature ID<TAB>Taxon<TAB>Confidence
  awk -F'\t' 'NR>1 && $2=="Unassigned"{print $1}' "$TAX_BOLD" > "$UNASS_IDS"
  N_UNASS=$(wc -l < "$UNASS_IDS")
  echo "[rescate] $N_UNASS ASVs Unassigned a reclasificar"

  if [[ "$N_UNASS" -eq 0 ]]; then
    echo "[rescate] Nada que rescatar — todos los ASVs ya asignados por BOLD"
    cp "$TAX_BOLD" "$TAXO_OUT/tax_consensus.tsv"
    echo "  → tax_consensus.tsv = copia de tax_BOLD.tsv"
    exit 0
  fi

  # Subset FASTA por IDs (con seqkit si existe, si no con awk)
  if command -v seqkit >/dev/null 2>&1; then
    seqkit grep -f "$UNASS_IDS" "$ASV_FASTA" > "$QUERY_FA"
  else
    python3 <<PY
ids = set(open("$UNASS_IDS").read().split())
keep = False
with open("$ASV_FASTA") as fi, open("$QUERY_FA","w") as fo:
    for line in fi:
        if line.startswith(">"):
            sid = line[1:].strip().split()[0]
            keep = sid in ids
        if keep:
            fo.write(line)
PY
  fi
  echo "[rescate] FASTA query: $QUERY_FA ($(grep -c '^>' "$QUERY_FA") secuencias)"
else
  QUERY_FA="$ASV_FASTA"
  echo "[all] FASTA query: $QUERY_FA (todos los ASVs)"
fi

# ---- 2) Import query a QIIME2 ----------------------------------------------
QUERY_QZA="$TAXO_OUT/asvs_query_ncbi_${MODE}.qza"
if [[ -s "$QUERY_QZA" ]]; then
  echo "[skip] $QUERY_QZA ya existe"
else
  echo "[run] import query"
  qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path "$QUERY_FA" \
    --output-path "$QUERY_QZA"
fi

# ---- 3) classify-consensus-vsearch vs NCBI ---------------------------------
CLASSIF_QZA="$TAXO_OUT/tax_NCBI_${MODE}_classification.qza"
SEARCH_QZA="$TAXO_OUT/tax_NCBI_${MODE}_search.qza"
if [[ -s "$CLASSIF_QZA" ]]; then
  echo "[skip] $CLASSIF_QZA ya existe"
else
  echo "[run] classify-consensus-vsearch vs NCBI"
  time qiime feature-classifier classify-consensus-vsearch \
    --i-query              "$QUERY_QZA" \
    --i-reference-reads    "$NCBI_SEQS_QZA" \
    --i-reference-taxonomy "$NCBI_TAX_QZA" \
    --p-perc-identity      "$PID" \
    --p-maxaccepts         "$MAXACC" \
    --p-min-consensus      "$MINCONS" \
    --p-strand both \
    --p-threads            "$THREADS" \
    --o-classification     "$CLASSIF_QZA" \
    --o-search-results     "$SEARCH_QZA" \
    --verbose 2>&1 | tee "logs/05_classify_ncbi_${MODE}.log"
fi

# ---- 4) Export TSV ---------------------------------------------------------
if [[ "$MODE" == "rescate" ]]; then
  TAX_NCBI_TSV="$TAXO_OUT/tax_NCBI_rescate.tsv"
else
  TAX_NCBI_TSV="$TAXO_OUT/tax_NCBI.tsv"
fi

TMPEXP=$(mktemp -d)
qiime tools export --input-path "$CLASSIF_QZA" --output-path "$TMPEXP"
cp "$TMPEXP/taxonomy.tsv" "$TAX_NCBI_TSV"
rm -rf "$TMPEXP"

# ---- 5) Consenso híbrido BOLD primary + NCBI rescate (solo MODE=rescate) --
if [[ "$MODE" == "rescate" ]]; then
  TAX_CONS="$TAXO_OUT/tax_consensus.tsv"
  echo "[rescate] Combinando BOLD (primary) + NCBI (rescate) → $TAX_CONS"
  python3 <<PY
import csv
bold = {}
with open("$TAX_BOLD") as f:
    rdr = csv.DictReader(f, delimiter="\t")
    for row in rdr:
        bold[row["Feature ID"]] = (row["Taxon"], row.get("Confidence",""))

ncbi = {}
with open("$TAX_NCBI_TSV") as f:
    rdr = csv.DictReader(f, delimiter="\t")
    for row in rdr:
        ncbi[row["Feature ID"]] = (row["Taxon"], row.get("Confidence",""))

with open("$TAX_CONS", "w") as fo:
    fo.write("Feature ID\tTaxon\tConfidence\tSource\n")
    n_bold = n_ncbi = n_unass = 0
    for fid, (tax, conf) in bold.items():
        if tax != "Unassigned":
            fo.write(f"{fid}\t{tax}\t{conf}\tBOLD\n")
            n_bold += 1
        elif fid in ncbi and ncbi[fid][0] != "Unassigned":
            t2, c2 = ncbi[fid]
            fo.write(f"{fid}\t{t2}\t{c2}\tNCBI_rescate\n")
            n_ncbi += 1
        else:
            fo.write(f"{fid}\tUnassigned\t\tnone\n")
            n_unass += 1
    print(f"BOLD primary: {n_bold}")
    print(f"NCBI rescate: {n_ncbi}")
    print(f"Unassigned:   {n_unass}")
PY
fi

# ---- Resumen ---------------------------------------------------------------
N_TOTAL=$(($(wc -l < "$TAX_NCBI_TSV") - 1))
N_UNASS=$(grep -c $'\tUnassigned\t' "$TAX_NCBI_TSV" 2>/dev/null || echo 0)
N_ASS=$((N_TOTAL - N_UNASS))
PCT_ASS=$(awk -v a=$N_ASS -v t=$N_TOTAL 'BEGIN{if(t>0)printf "%.1f", 100*a/t; else print "0.0"}')

echo
echo "=========================================================="
echo " Clasificación NCBI completa (MODE=$MODE) — $(date -Iseconds)"
echo " Consultados:  $N_TOTAL ASVs"
echo " Asignados:    $N_ASS ($PCT_ASS %)"
echo " Unassigned:   $N_UNASS"
echo "----------------------------------------------------------"
echo " Salidas:"
echo "  Classification: $CLASSIF_QZA"
echo "  Taxonomía:      $TAX_NCBI_TSV"
if [[ "$MODE" == "rescate" ]]; then
  echo "  Consenso final: $TAXO_OUT/tax_consensus.tsv (BOLD + NCBI rescate)"
fi
echo "----------------------------------------------------------"
echo " Acciones recomendadas:"
echo "  - Inspeccionar profundidad de asignación por rango:"
echo "      awk -F'\\t' 'NR>1 && \$2!=\"Unassigned\"{n=split(\$2,a,\";\"); print n}' $TAX_NCBI_TSV | sort -n | uniq -c"
if [[ "$MODE" == "rescate" ]]; then
  echo "  - Comparar % asignación por fuente:"
  echo "      awk -F'\\t' 'NR>1{print \$NF}' $TAXO_OUT/tax_consensus.tsv | sort | uniq -c"
fi
echo "  - Continuar con curación swarm+LULU: workflow/05_swarm_LULU.R"
echo "=========================================================="
