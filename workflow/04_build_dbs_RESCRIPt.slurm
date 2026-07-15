#!/usr/bin/env bash
# =============================================================================
# 04_build_dbs_RESCRIPt.slurm — ESTRATEGIA ACTIVA: BOLD (SciLifeLab) + NCBI
# -----------------------------------------------------------------------------
# 2026-07-14: RESCRIPt actual NO expone `get-bold-data`. Sustituimos 04a por
# `04a_build_BOLD_db_SciLifeLab.sh` que descarga el dataset SciLifeLab v6
# (DOI 10.17044/scilifelab.20514192.v6, derivado del BOLD Data Package
# 15-May-2026, ~24 M secuencias, ya curado y en formato QIIME2).
# 04b (NCBI vía RESCRIPt) se mantiene sin cambios.
# -----------------------------------------------------------------------------
# Tiempo estimado:
#   04a (SciLifeLab): 20–90 min (descarga 416 MB + import + Naive Bayes)
#   04b (NCBI):       15–40 min con NCBI_API_KEY, 1–4 h sin ella
#   Walltime 06:00:00 sigue siendo suficiente.
# -----------------------------------------------------------------------------
# Uso:
#   cd /ruta/al/repo/LANCO
#   sbatch workflow/04_build_dbs_RESCRIPt.slurm
#
# Opcional (acelera NCBI ~10x, muy recomendado):
#   sbatch --export=ALL,NCBI_API_KEY=xxxxxxxx workflow/04_build_dbs_RESCRIPt.slurm
# =============================================================================
#SBATCH --job-name=lanco_rescript_dbs
# IMPORTANTE: 'logs/' YA existe en el repo; NO uses subcarpetas que no existan,
# porque SLURM abre estos archivos antes de ejecutar el script y el job muere
# en silencio si la ruta no es escribible.
# Revisar particiones: sinfo -s
#SBATCH --output=logs/04_rescript_%x_%j.out
#SBATCH --error=logs/04_rescript_%x_%j.err
#SBATCH --time=06:00:00
#SBATCH --partition=cicese
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rgomez41@uabc.edu.mx

set -euo pipefail

# ---- Entorno conda (cluster LUSTRE — Anaconda 2026) ------------------------
# ~/.local/bin va PRIMERO para que yq (instalado por el usuario, sin sudo) tenga
# prioridad. Los jobs SLURM NO cargan ~/.bashrc por defecto, por eso lo añadimos
# explícitamente aquí.
export PATH=$HOME/.local/bin:/LUSTRE/apps/Anaconda/2026/miniconda3/condabin:/LUSTRE/apps/Anaconda/2026/miniconda3/bin:$PATH
# shellcheck disable=SC1091
source activate qiime2-2026

# ---- Rutas ------------------------------------------------------------------
# CRÍTICO: SLURM COPIA el script al spool (/var/spool/slurm/job_XXX/slurm_script)
# y lo ejecuta desde allí. Por eso ${BASH_SOURCE[0]} apunta al spool, NO al repo,
# y `cd $ROOT` termina en /var/spool/slurm (sin workflow/, sin permisos).
# La única forma correcta de localizar el repo es $SLURM_SUBMIT_DIR: el dir
# desde donde se lanzó `sbatch`.
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$ROOT" || { echo "[FATAL] No puedo entrar a ROOT=$ROOT"; exit 1; }

# Sanity check: el repo debe contener workflow/04a y workflow/04b
if [[ ! -f workflow/04a_build_BOLD_db_RESCRIPt.sh || ! -f workflow/04b_build_NCBI_db_RESCRIPt.sh ]]; then
  echo "[FATAL] No encuentro workflow/04a o 04b desde ROOT=$ROOT"
  echo "        ¿lanzaste 'sbatch workflow/04_build_dbs_RESCRIPt.slurm' desde la raíz del repo?"
  echo "        Contenido actual:"; ls -la
  exit 1
fi

mkdir -p logs || { echo "[FATAL] No puedo crear logs/ en $ROOT (permisos)"; exit 1; }

# Marcador visible al inicio: si esta línea NO aparece en el .out, es que el
# script nunca arrancó (típicamente: directorio de output inexistente,
# partición/QoS mal, o el env conda no activó).
echo "[$(date -Iseconds)] SLURM job arrancó en $(hostname):$(pwd)"

# ---- Localizar yq (instalación local del usuario, sin sudo) ----------------
# yq NO viene con qiime2; es obligatorio porque 04a/04b leen params.yml con él.
# El usuario lo instaló en ~/.local/bin/yq (ya añadido al PATH arriba).
YQ_LOCAL="$HOME/.local/bin/yq"
if [[ -x "$YQ_LOCAL" ]]; then
  echo "[setup] yq local detectado: $YQ_LOCAL"
elif command -v yq >/dev/null 2>&1; then
  echo "[setup] yq encontrado en PATH: $(command -v yq)"
else
  echo "[FATAL] yq no está en \$HOME/.local/bin ni en PATH."           >&2
  echo "        Instálalo (sin sudo) con:"                              >&2
  echo "          mkdir -p ~/.local/bin"                                >&2
  echo "          curl -fsSL https://github.com/mikefarah/yq/releases/download/v4.44.3/yq_linux_amd64 -o ~/.local/bin/yq" >&2
  echo "          chmod +x ~/.local/bin/yq"                             >&2
  exit 1
fi

# ---- Verificaciones rápidas -------------------------------------------------
echo "=========================================================="
echo " JobID:       ${SLURM_JOB_ID:-N/A}"
echo " Node:        $(hostname)"
echo " Fecha:       $(date -Iseconds)"
echo " Conda env:   ${CONDA_DEFAULT_ENV:-unset}"
echo " which qiime: $(command -v qiime || echo 'NOT FOUND')"
echo " which yq:    $(command -v yq    || echo 'NOT FOUND')"
echo " yq version:  $(yq --version 2>/dev/null || echo 'N/A')"
echo " CPUs:        ${SLURM_CPUS_PER_TASK:-?}"
echo " Mem:         ${SLURM_MEM_PER_NODE:-?} MB"
echo " NCBI_API_KEY: $([[ -n "${NCBI_API_KEY:-}" ]] && echo 'set' || echo 'unset')"
echo "=========================================================="

command -v qiime >/dev/null || { echo "[FATAL] qiime no en PATH — ¿está activado qiime2-2026?"; exit 1; }
command -v yq    >/dev/null || { echo "[FATAL] yq no en PATH (falló autoinstalación)"; exit 1; }

# yq de mikefarah v4 requiere >= 4.18 para invocar expresiones sin 'eval' explícito
# (04a/04b usan sintaxis `yq '.path' file.yml`). Versiones antiguas fallan con
# "unknown command '.foo' for yq". Validamos que sea 4.18+.
YQ_VER_STR=$(yq --version 2>/dev/null | grep -oE 'v?[0-9]+\.[0-9]+\.[0-9]+' | tail -1)
YQ_MAJOR=$(echo "${YQ_VER_STR#v}" | cut -d. -f1)
YQ_MINOR=$(echo "${YQ_VER_STR#v}" | cut -d. -f2)
if [[ "$YQ_MAJOR" -lt 4 ]] || { [[ "$YQ_MAJOR" -eq 4 ]] && [[ "$YQ_MINOR" -lt 18 ]]; }; then
  echo "[FATAL] yq ${YQ_VER_STR} es muy antiguo. Los scripts 04a/04b requieren >= v4.18."
  echo "        Actualiza (sin sudo) con:"
  echo "          curl -fsSL https://github.com/mikefarah/yq/releases/download/v4.44.3/yq_linux_amd64 -o ~/.local/bin/yq"
  echo "          chmod +x ~/.local/bin/yq"
  exit 1
fi

qiime rescript --help >/dev/null 2>&1 || { echo "[FATAL] Plugin RESCRIPt no instalado en qiime2-2026"; exit 1; }

# ---- 04a  BOLD (SciLifeLab, sustituye a RESCRIPt get-bold-data) ------------
echo
echo ">>> [$(date +%T)] Iniciando 04a — BOLD (SciLifeLab v6 / BDP 15-May-2026)"
time bash workflow/04a_build_BOLD_db_SciLifeLab.sh
echo ">>> [$(date +%T)] 04a OK"

# ---- 04b  NCBI --------------------------------------------------------------
echo
echo ">>> [$(date +%T)] Iniciando 04b — NCBI (RESCRIPt)"
time bash workflow/04b_build_NCBI_db_RESCRIPt.sh
echo ">>> [$(date +%T)] 04b OK"

echo
echo "=========================================================="
echo " Pipeline 04 (BOLD + NCBI) completo — $(date -Iseconds)"
echo " Salidas BOLD: db/rescript/bold-coi-*.qza / *.fasta / *.tsv"
echo " Salidas NCBI: db/rescript/ncbi-coi-*.qza / *.fasta / *.tsv"
echo " Siguiente paso: workflow/04_taxonomy_BOLD_NCBI.R (hybrid assignment)"
echo "=========================================================="
