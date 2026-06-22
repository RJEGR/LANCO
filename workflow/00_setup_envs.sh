#!/usr/bin/env bash
# =============================================================================
# 00_setup_envs.sh
# Bootstrapping de entornos conda + paquetes R para el workflow LANCO/RUN29
# -----------------------------------------------------------------------------
# Crea (idempotentemente) dos entornos conda y verifica todas las dependencias:
#
#   1) lanco_coi         → pipeline DADA2 (R), cutadapt, vsearch, utilidades
#   2) qiime2-rescript   → QIIME2 + plugin RESCRIPt (solo pasos 04a/04b)
#
# Uso:
#   bash workflow/00_setup_envs.sh                  # crea + verifica todo
#   bash workflow/00_setup_envs.sh lanco            # solo lanco_coi
#   bash workflow/00_setup_envs.sh rescript         # solo qiime2-rescript
#   bash workflow/00_setup_envs.sh verify           # verifica sin reinstalar
#
# Requisitos previos del sistema:
#   - mamba o conda instalado
#   - bash >= 4
# =============================================================================

set -euo pipefail

# ---- Config ----------------------------------------------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT="$( cd "${SCRIPT_DIR}/.." && pwd )"
LOG_DIR="${ROOT}/logs/00_setup"
mkdir -p "$LOG_DIR"

ENV_LANCO="lanco_coi"
ENV_QIIME="qiime2-rescript"

# Versiones pineadas (ajustar si Ricardo quiere actualizar)
QIIME_RELEASE="2024.10"                      # QIIME2 Amplicon distribution
PY_LANCO="3.11"
PY_QIIME="3.10"                              # QIIME2 2024.10 usa 3.10

# Channels
CHANNELS_BIOCONDA="-c bioconda -c conda-forge"
CHANNELS_QIIME="-c qiime2 -c conda-forge -c bioconda -c defaults"

# Mamba o conda
if command -v mamba >/dev/null 2>&1; then
  CONDA=mamba
elif command -v conda >/dev/null 2>&1; then
  CONDA=conda
  echo "[WARN] mamba no encontrado, usando conda (más lento). Se recomienda instalar mamba."
else
  echo "[ERROR] Ni mamba ni conda están en PATH. Instala Miniforge desde:"
  echo "        https://github.com/conda-forge/miniforge"
  exit 1
fi

# Habilitar 'conda activate' dentro del script
CONDA_BASE=$(conda info --base 2>/dev/null || echo "$HOME/miniforge3")
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Algunos hooks de activación de conda (p.ej. activate-gfortran_osx-64.sh)
# referencian variables sin definir y rompen bajo 'set -u'; se desactiva
# nounset temporalmente alrededor de activate/deactivate.
activate_env () {
  set +u
  conda activate "$1"
  set -u
}

deactivate_env () {
  set +u
  conda deactivate
  set -u
}

# ---- Helpers ---------------------------------------------------------------
env_exists () {
  conda env list | awk '{print $1}' | grep -qx "$1"
}

log_run () {
  # log_run <label> <cmd...>
  local label="$1"; shift
  local logfile="${LOG_DIR}/${label}_$(date +%Y%m%d_%H%M%S).log"
  echo "[run]  ${label}  →  log: ${logfile}"
  "$@" 2>&1 | tee "$logfile"
}

# El "yq" de conda-forge/defaults es la herramienta de kislyuk (wrapper de jq,
# tope en 3.4.3); mikefarah/yq v4 -el que usan 02/04a/04b- nunca se publicó
# en conda, así que se instala el binario oficial directo en el bin del env activo.
install_yq_binary () {
  if command -v yq >/dev/null 2>&1 && yq --version 2>&1 | grep -q "version v4"; then
    echo "[skip] yq v4 ya está instalado ($(yq --version 2>&1))"
    return
  fi

  local os arch
  case "$(uname -s)" in
    Linux*)  os="linux"  ;;
    Darwin*) os="darwin" ;;
    *) echo "[ERROR] OS no soportado para yq: $(uname -s)"; exit 1 ;;
  esac
  case "$(uname -m)" in
    x86_64|amd64)  arch="amd64" ;;
    arm64|aarch64) arch="arm64" ;;
    *) echo "[ERROR] Arquitectura no soportada para yq: $(uname -m)"; exit 1 ;;
  esac

  local dest="${CONDA_PREFIX}/bin/yq"
  local url="https://github.com/mikefarah/yq/releases/latest/download/yq_${os}_${arch}"
  echo "[run]  Descargando yq v4 (${os}/${arch}) → ${dest}"
  curl -fsSL -o "$dest" "$url"
  chmod +x "$dest"
}

# ---- 1) Entorno lanco_coi --------------------------------------------------
setup_lanco () {
  echo
  echo "============================================================"
  echo " Entorno 1/2 :  ${ENV_LANCO}"
  echo "============================================================"

  if env_exists "$ENV_LANCO"; then
    echo "[skip] ${ENV_LANCO} ya existe (usa 'conda env remove -n ${ENV_LANCO}' para recrear)"
  else
    log_run "create_${ENV_LANCO}" $CONDA create -n "$ENV_LANCO" -y $CHANNELS_BIOCONDA \
      "python=${PY_LANCO}"   \
      "cutadapt=4.*"         \
      "vsearch=2.*"          \
      "blast=2.*"            \
      "biopython"            \
      "r-base=4.3.*"         \
      "r-yaml"               \
      "r-dplyr"              \
      "r-readr"              \
      "r-tibble"             \
      "r-tidyr"              \
      "r-ggplot2"            \
      "r-rprojroot"          \
      "r-biocmanager"
  fi

  activate_env "$ENV_LANCO"
  install_yq_binary

  # Instalar paquetes Bioconductor desde R (solo si no están)
  echo "[run]  Verificando/instalando paquetes R (Bioconductor)..."
  Rscript - <<'RSCRIPT' 2>&1 | tee "${LOG_DIR}/R_packages_$(date +%Y%m%d_%H%M%S).log"
needed <- c("dada2", "Biostrings", "ShortRead")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  cat("[R] Instalando vía BiocManager:", paste(missing, collapse=", "), "\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
  BiocManager::install(missing, ask = FALSE, update = FALSE)
} else {
  cat("[R] Todos los paquetes Bioconductor ya están instalados.\n")
}
# Reporte de versiones
cat("\n[R] Versiones:\n")
for (p in c("dada2","Biostrings","ShortRead","yaml","dplyr","readr","ggplot2")) {
  cat(sprintf("    %-12s %s\n", p, as.character(packageVersion(p))))
}
RSCRIPT
  deactivate_env
}

# ---- 2) Entorno qiime2-rescript --------------------------------------------
setup_qiime_rescript () {
  echo
  echo "============================================================"
  echo " Entorno 2/2 :  ${ENV_QIIME}  (QIIME2 ${QIIME_RELEASE} + RESCRIPt)"
  echo "============================================================"

  if env_exists "$ENV_QIIME"; then
    echo "[skip] ${ENV_QIIME} ya existe"
  else
    # QIIME2 se instala desde su YAML oficial para garantizar compatibilidad de versiones
    local OS_TAG
    case "$(uname -s)" in
      Linux*)  OS_TAG="linux-conda" ;;
      Darwin*) OS_TAG="osx-conda"   ;;
      *) echo "[ERROR] OS no soportado por QIIME2: $(uname -s)"; exit 1 ;;
    esac

    local QIIME_YML="qiime2-amplicon-${QIIME_RELEASE}-py310-${OS_TAG}.yml"
    local QIIME_URL="https://data.qiime2.org/distro/amplicon/${QIIME_YML}"

    echo "[run]  Descargando manifest QIIME2: ${QIIME_URL}"
    curl -fsSL -o "/tmp/${QIIME_YML}" "$QIIME_URL"

    # El canal de QIIME2 solo publica builds osx-64 (sin osx-arm64 nativo);
    # en Apple Silicon hay que forzar el subdir x86_64 y resolver vía Rosetta 2,
    # si no el solver reporta casi todos los paquetes qiime2/q2-*/bioconductor-* como "missing channel".
    if [[ "$(uname -s)" == "Darwin" && "$(uname -m)" == "arm64" ]]; then
      if ! arch -x86_64 /usr/bin/true >/dev/null 2>&1; then
        echo "[ERROR] Rosetta 2 no está instalado. Requerido para QIIME2 en Apple Silicon."
        echo "        Instálalo con: softwareupdate --install-rosetta"
        exit 1
      fi
      echo "[info] Apple Silicon detectado → forzando CONDA_SUBDIR=osx-64 (Rosetta 2) para ${ENV_QIIME}"
      log_run "create_${ENV_QIIME}" env CONDA_SUBDIR=osx-64 $CONDA env create -n "$ENV_QIIME" --file "/tmp/${QIIME_YML}"
      activate_env "$ENV_QIIME"
      conda config --env --set subdir osx-64
      deactivate_env
    else
      log_run "create_${ENV_QIIME}" $CONDA env create -n "$ENV_QIIME" --file "/tmp/${QIIME_YML}"
    fi
  fi

  # Instalar RESCRIPt vía pip (no está en el manifest oficial)
  activate_env "$ENV_QIIME"
  install_yq_binary
  if pip show rescript >/dev/null 2>&1; then
    echo "[skip] RESCRIPt ya instalado en ${ENV_QIIME}"
  else
    echo "[run]  pip install RESCRIPt..."
    log_run "rescript_install" pip install git+https://github.com/bokulich-lab/RESCRIPt.git
    # Refresh QIIME2 cache para que detecte el plugin
    qiime dev refresh-cache
  fi
  deactivate_env
}

# ---- 3) Verificación -------------------------------------------------------
verify () {
  echo
  echo "============================================================"
  echo " Verificación de dependencias"
  echo "============================================================"
  local errors=0

  # --- lanco_coi
  if env_exists "$ENV_LANCO"; then
    activate_env "$ENV_LANCO"
    echo
    echo "[$ENV_LANCO]"
    for tool in cutadapt vsearch blastn yq python R Rscript; do
      if command -v "$tool" >/dev/null 2>&1; then
        printf "  %-12s OK  (%s)\n" "$tool" "$($tool --version 2>&1 | head -1)"
      else
        printf "  %-12s MISSING\n" "$tool"; errors=$((errors+1))
      fi
    done
    echo "  R packages:"
    Rscript -e '
      pkgs <- c("dada2","Biostrings","ShortRead","yaml","dplyr","readr","ggplot2","rprojroot")
      for (p in pkgs) {
        ok <- suppressWarnings(requireNamespace(p, quietly=TRUE))
        cat(sprintf("    %-12s %s\n", p, if (ok) paste("OK",packageVersion(p)) else "MISSING"))
      }
    ' || errors=$((errors+1))
    deactivate_env
  else
    echo "[MISSING] entorno $ENV_LANCO no existe"
    errors=$((errors+1))
  fi

  # --- qiime2-rescript
  if env_exists "$ENV_QIIME"; then
    activate_env "$ENV_QIIME"
    echo
    echo "[$ENV_QIIME]"
    if command -v qiime >/dev/null 2>&1; then
      printf "  %-12s OK  (%s)\n" "qiime" "$(qiime --version 2>&1 | head -1)"
    else
      printf "  %-12s MISSING\n" "qiime"; errors=$((errors+1))
    fi
    if qiime rescript --help >/dev/null 2>&1; then
      printf "  %-12s OK\n" "rescript"
    else
      printf "  %-12s MISSING (pip install git+https://github.com/bokulich-lab/RESCRIPt.git)\n" "rescript"
      errors=$((errors+1))
    fi
    deactivate_env
  else
    echo "[MISSING] entorno $ENV_QIIME no existe"
    errors=$((errors+1))
  fi

  # --- NCBI API key (opcional pero recomendado)

  export NCBI_API_KEY=74e0e4bf2d8eaa3bb742c46316dbafe12909
  
  echo
  if [[ -n "${NCBI_API_KEY:-}" ]]; then
    echo "  NCBI_API_KEY  OK (~10x speedup en 04b)"
  else
    echo "  NCBI_API_KEY  no configurado (opcional)"
    echo "                Para acelerar 04b: export NCBI_API_KEY=tu_key"
    echo "                Obtener en: https://www.ncbi.nlm.nih.gov/account/settings/"
  fi

  echo
  if [[ $errors -eq 0 ]]; then
    echo "✓ Todas las dependencias OK. Listo para ejecutar el pipeline."
  else
    echo "✗ Se encontraron $errors problemas. Revisa los mensajes arriba."
    exit 1
  fi
}

# ---- Main -------------------------------------------------------------------
ACTION="${1:-all}"

case "$ACTION" in
  all)
    setup_lanco
    setup_qiime_rescript
    verify
    ;;
  lanco)
    setup_lanco
    ;;
  rescript)
    setup_qiime_rescript
    ;;
  verify)
    verify
    ;;
  *)
    echo "Uso: $0 [all|lanco|rescript|verify]"
    exit 1
    ;;
esac

echo
echo "============================================================"
echo " Próximo paso:"
echo "   conda activate ${ENV_LANCO}"
echo "   Rscript workflow/01_inspect_reads.R"
echo "============================================================"
