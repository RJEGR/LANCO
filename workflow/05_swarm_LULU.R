# =============================================================================
# 05_swarm_LULU.R   — colapso de ASVs en OTUs + curación de spurious clusters
# LANCO/RUN29 COI Leray-XT  (post-DADA2 SE)
# -----------------------------------------------------------------------------
# Implementa el pipeline post-DADA2 propuesto por Brandt et al. 2021
# (Mol Ecol Resour) y Antich et al. 2021 (BMC Bioinf):
#
#   1) Dereplicación + formato swarm (>ASV_NNNN_TotalAbundance)
#   2) swarm v3 clustering (d=13 para metazoo marino) → OTU centroids
#   3) Colapso de asv_table por cluster → otu_table_swarm
#   4) vsearch usearch_global all-vs-all → match_list (input LULU)
#   5) LULU curation → eliminación de daughter-OTUs espurios
#   6) Tracking ASVs → swarm-OTUs → LULU-OTUs
#
# Entrada:  results/03_dada2_SE/{asvs.fasta, asv_table.tsv}
# Salidas en results/05_swarm_lulu/:
#   - swarm/swarm_input.fasta              (dereplicated con abundance suffix)
#   - swarm/swarm_clusters.txt             (un cluster por línea, members)
#   - swarm/swarm_seeds.fasta              (centroids = OTUs representativos)
#   - swarm/swarm_struct.txt               (estructura interna por cluster)
#   - swarm/swarm_stats.txt                (theoretical_radius, etc.)
#   - otu_table_swarm.tsv                  (38 muestras × N_OTUs)
#   - lulu/match_list.txt                  (vsearch pairs ≥84% identity)
#   - lulu/curated_otu_table.tsv           (post-LULU, listo para taxonomía)
#   - lulu/curated_otus.fasta              (secuencias OTUs curados)
#   - lulu/lulu_summary.tsv                (decisión por OTU: kept|merged)
#   - track_clustering.tsv                 (resumen reducción ASVs→swarm→LULU)
# -----------------------------------------------------------------------------
# DEPENDENCIAS EXTERNAS (verificadas al inicio del script):
#   - swarm v3+   :  ver docs/INSTALL_dependencies.md
#   - vsearch v2+ :  ver docs/INSTALL_dependencies.md
#   - R-pkg lulu  :  devtools::install_github("tobiasgf/lulu")
#
# INSTALACIÓN validada 2026-06-23 (macOS osx-64):
#   conda create -n coi_pipeline -c conda-forge -c bioconda \
#          swarm=3.0 vsearch python-igraph -y
#   conda activate coi_pipeline
# NOTA: `conda install -c bioconda swarm` sin conda-forge falla por bug
# de metadata (dep falsa a python-igraph). El env aislado con conda-forge
# prioritario resuelve el conflicto.
# =============================================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(yaml)
  library(rprojroot)
})

ROOT   <- find_root(has_dir("RUN29"))
params <- yaml::read_yaml(file.path(ROOT, "config", "params.yml"))

# ---- 0) Verificación de dependencias externas ------------------------------
check_bin <- function(bin) {
  path <- Sys.which(bin)
  if (path == "") stop(sprintf(
    "Binario '%s' no encontrado en PATH.\n  → instalar con:\n    macOS: brew install %s\n    conda: conda install -c bioconda %s",
    bin, bin, bin))
  message(sprintf("[05] %-8s → %s", bin, path))
  path
}
SWARM   <- check_bin("swarm")
VSEARCH <- check_bin("vsearch")

if (!requireNamespace("lulu", quietly = TRUE)) {
  stop("Paquete R 'lulu' no instalado.\n",
       "  → instalar con: install.packages('devtools'); ",
       "devtools::install_github('tobiasgf/lulu')")
}
suppressPackageStartupMessages(library(lulu))

# ---- 1) Rutas y carga de inputs --------------------------------------------
in_dir  <- file.path(ROOT, params$paths$dada2_se_out)
out_dir <- file.path(ROOT, params$paths$swarm_lulu_out)
sw_dir  <- file.path(out_dir, "swarm")
lu_dir  <- file.path(out_dir, "lulu")
dir.create(sw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(lu_dir, showWarnings = FALSE, recursive = TRUE)

fa_in   <- file.path(in_dir, "asvs.fasta")
tab_in  <- file.path(in_dir, "asv_table.tsv")
stopifnot(file.exists(fa_in), file.exists(tab_in))

asv_seqs <- readDNAStringSet(fa_in)
asv_tab  <- read_tsv(tab_in, show_col_types = FALSE)
sample_ids <- asv_tab$sample_id
asv_mat <- as.matrix(asv_tab[, -1, drop = FALSE])
rownames(asv_mat) <- sample_ids   # filas=muestras, columnas=ASVs

message(sprintf("[05] Input: %d ASVs × %d muestras", ncol(asv_mat), nrow(asv_mat)))
stopifnot(ncol(asv_mat) == length(asv_seqs),
          all(colnames(asv_mat) == names(asv_seqs)))

# ---- 2) Preparar input swarm (header con abundance) ------------------------
# swarm requiere headers ">id_abundance" con la abundancia TOTAL del ASV.
total_abund <- colSums(asv_mat)
names(asv_seqs) <- sprintf("%s_%d", names(asv_seqs), total_abund)
swarm_input <- file.path(sw_dir, "swarm_input.fasta")
writeXStringSet(asv_seqs, swarm_input)
message("[05] swarm input escrito: ", swarm_input)

# ---- 3) Ejecutar swarm ------------------------------------------------------
sw          <- params$swarm_lulu$swarm
swarm_args  <- c("-d", sw$d, "-t", sw$threads,
                 "-o", file.path(sw_dir, "swarm_clusters.txt"),
                 "-s", file.path(sw_dir, "swarm_stats.txt"),
                 "-i", file.path(sw_dir, "swarm_struct.txt"),
                 "-w", file.path(sw_dir, "swarm_seeds.fasta"))
if (isTRUE(sw$fastidious) && sw$d == 1) swarm_args <- c(swarm_args, "-f")
swarm_args <- c(swarm_args, swarm_input)

message(sprintf("[05] swarm d=%d threads=%d fastidious=%s",
                sw$d, sw$threads, sw$fastidious))
res <- system2(SWARM, args = swarm_args, stdout = TRUE, stderr = TRUE)
writeLines(res, file.path(sw_dir, "swarm_log.txt"))
stopifnot(file.exists(file.path(sw_dir, "swarm_clusters.txt")))

# ---- 4) Parsear clusters → mapping ASV → swarm OTU -------------------------
# Cada línea = un cluster; miembros separados por espacios; primer miembro
# es el seed (representativo). Los IDs incluyen "_abundance" → strip.
strip_abund <- function(x) sub("_\\d+$", "", x)
cluster_lines <- readLines(file.path(sw_dir, "swarm_clusters.txt"))
cluster_list  <- strsplit(cluster_lines, " ", fixed = TRUE)
n_clusters    <- length(cluster_list)
otu_ids       <- sprintf("OTU_%04d", seq_len(n_clusters))

asv_to_otu <- tibble(
  asv = strip_abund(unlist(cluster_list)),
  otu = rep(otu_ids, lengths(cluster_list)),
  seed = rep(vapply(cluster_list, function(x) strip_abund(x[1]), character(1)),
             lengths(cluster_list))
)
write_tsv(asv_to_otu, file.path(out_dir, "asv_to_otu_map.tsv"))
message(sprintf("[05] swarm produjo %d OTUs desde %d ASVs (reducción %.1f%%)",
                n_clusters, ncol(asv_mat),
                100 * (1 - n_clusters / ncol(asv_mat))))

# ---- 5) Colapsar asv_table → otu_table_swarm -------------------------------
otu_mat <- matrix(0L, nrow = nrow(asv_mat), ncol = n_clusters,
                  dimnames = list(rownames(asv_mat), otu_ids))
for (i in seq_along(cluster_list)) {
  members <- strip_abund(cluster_list[[i]])
  members <- intersect(members, colnames(asv_mat))
  if (length(members) == 0) next
  otu_mat[, i] <- rowSums(asv_mat[, members, drop = FALSE])
}
write_tsv(as_tibble(otu_mat, rownames = "sample_id"),
          file.path(out_dir, "otu_table_swarm.tsv"))

# Renombrar seeds.fasta con IDs OTU_NNNN estables (swarm los guarda con _abund)
seeds_raw <- readDNAStringSet(file.path(sw_dir, "swarm_seeds.fasta"))
names(seeds_raw) <- strip_abund(names(seeds_raw))
# Reorden para que coincida con otu_ids (cluster order = línea en clusters.txt)
seed_per_cluster <- vapply(cluster_list, function(x) strip_abund(x[1]), character(1))
seeds_ord <- seeds_raw[seed_per_cluster]
names(seeds_ord) <- otu_ids
seeds_fa <- file.path(sw_dir, "swarm_seeds_renamed.fasta")
writeXStringSet(seeds_ord, seeds_fa)

# ---- 6) vsearch usearch_global all-vs-all → match_list para LULU -----------
vs          <- params$swarm_lulu$vsearch
match_list  <- file.path(lu_dir, "match_list.txt")
vs_args     <- c("--usearch_global", seeds_fa,
                 "--db",            seeds_fa,
                 "--self",
                 "--id",            sprintf("%.2f", vs$min_match_pct / 100),
                 "--iddef",         vs$iddef,
                 "--userout",       match_list,
                 "--userfields",    "query+target+id",
                 "--maxaccepts",    "0",
                 "--query_cov",     "0.9",
                 "--maxhits",       "10",
                 "--threads",       vs$threads)
message(sprintf("[05] vsearch all-vs-all (id>=%d%%, iddef=%d)",
                vs$min_match_pct, vs$iddef))
vs_log <- system2(VSEARCH, args = vs_args, stdout = TRUE, stderr = TRUE)
writeLines(vs_log, file.path(lu_dir, "vsearch_log.txt"))
stopifnot(file.exists(match_list))

# ---- 7) LULU curation -------------------------------------------------------
# lulu::lulu requiere otu_table con FILAS=OTUs, COLUMNAS=muestras (no como DADA2).
lu      <- params$swarm_lulu$lulu
ml_tab  <- read.table(match_list, header = FALSE, as.is = TRUE,
                      stringsAsFactors = FALSE)
# Algunas instalaciones vsearch escriben 0 filas si no hay matches → fallback.
if (nrow(ml_tab) == 0) {
  warning("[05] match_list vacío: ningún par OTU-OTU superó ",
          vs$min_match_pct, "% identidad. LULU no puede curar; copiando ",
          "otu_table_swarm como output final.")
  curated_table <- t(otu_mat)
  curated_otus  <- seeds_ord
  lulu_summary  <- tibble(otu = otu_ids, decision = "kept", parent = NA_character_)
} else {
  otu_table_for_lulu <- as.data.frame(t(otu_mat))   # filas=OTUs
  message(sprintf("[05] LULU: minimum_match=%d minimum_cooc=%.2f ratio_type=%s",
                  lu$minimum_match, lu$minimum_relative_cooccurence,
                  lu$minimum_ratio_type))
  curated <- lulu::lulu(
    otutable                       = otu_table_for_lulu,
    matchlist                      = ml_tab,
    minimum_ratio_type             = lu$minimum_ratio_type,
    minimum_ratio                  = lu$minimum_ratio,
    minimum_match                  = lu$minimum_match,
    minimum_relative_cooccurence   = lu$minimum_relative_cooccurence
  )
  message(sprintf("[05] LULU: %d OTUs retenidos, %d colapsados",
                  curated$curated_count, curated$discarded_count))

  curated_table <- as.matrix(curated$curated_table)   # filas=OTUs
  curated_otus  <- seeds_ord[rownames(curated_table)]
  lulu_summary  <- tibble(
    otu      = rownames(otu_table_for_lulu),
    decision = ifelse(rownames(otu_table_for_lulu) %in% rownames(curated_table),
                      "kept", "merged"),
    parent   = curated$otu_map$parent_id[match(rownames(otu_table_for_lulu),
                                               rownames(curated$otu_map))]
  )
}

# ---- 8) Escribir outputs finales -------------------------------------------
write_tsv(as_tibble(t(curated_table), rownames = "sample_id"),
          file.path(lu_dir, "curated_otu_table.tsv"))
writeXStringSet(curated_otus, file.path(lu_dir, "curated_otus.fasta"))
write_tsv(lulu_summary, file.path(lu_dir, "lulu_summary.tsv"))

# ---- 9) Tracking de reducción ----------------------------------------------
track <- tibble(
  stage = c("ASVs (DADA2 SE)", "OTUs (swarm d=", "OTUs curados (LULU)"),
  count = c(ncol(asv_mat), n_clusters, ncol(t(curated_table))),
  reads_total = c(sum(asv_mat), sum(otu_mat), sum(curated_table))
)
track$stage[2] <- sprintf("OTUs (swarm d=%d)", sw$d)
write_tsv(track, file.path(out_dir, "track_clustering.tsv"))

cat("\n=== Resumen 05_swarm_LULU ===\n")
print(as.data.frame(track), row.names = FALSE)
cat(sprintf("\nReducción total ASVs → OTUs curados: %d → %d (%.1f%%)\n",
            ncol(asv_mat), nrow(curated_table),
            100 * (1 - nrow(curated_table) / ncol(asv_mat))))
cat("\nOutputs en:", out_dir, "\n")
cat("\nAcciones recomendadas:\n")
cat(" - Inspeccionar lulu_summary.tsv: confirmar que merges tienen sentido.\n")
cat(" - Continuar con asignación taxonómica sobre:\n")
cat("     ", file.path(lu_dir, "curated_otus.fasta"), "\n")
cat(" - Comparar diversidad alfa (Shannon) entre asvs (sin curar) y OTUs curados.\n")
