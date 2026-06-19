# =============================================================================
# 04_taxonomy_BOLD_NCBI.R
# Asignación taxonómica de ASVs COI usando DBs RESCRIPt (BOLD + NCBI)
# -----------------------------------------------------------------------------
# Prerrequisitos (ejecutar una sola vez):
#   workflow/04a_build_BOLD_db_RESCRIPt.sh   → genera bold-coi-{seqs,tax}.{qza,fasta,tsv}
#   workflow/04b_build_NCBI_db_RESCRIPt.sh   → genera ncbi-coi-{seqs,tax}.{qza,fasta,tsv}
# -----------------------------------------------------------------------------
# Estrategia híbrida (params: taxonomy.strategy):
#   "vsearch_consensus" (DEFAULT, recomendado):
#       VSEARCH --usearch_global de ASVs contra BOLD → consenso top hits.
#       ASVs no asignados a Género o sin match → mismo procedimiento vs NCBI.
#       Es el método más reproducible y robusto a gaps en COI.
#   "naive_bayes":
#       Reinvoca QIIME2 (classify-sklearn) usando los .qza de RESCRIPt.
#       Requiere qiime2-rescript activo en el shell que llama a este script.
#   "dada2_R":
#       assignTaxonomy() de DADA2 sobre el FASTA exportado (headers DADA2-format).
#       Más simple pero más lento si la DB es grande.
# -----------------------------------------------------------------------------
# Salidas en results/04_taxonomy/:
#   - tax_BOLD.tsv          (asignación primaria)
#   - tax_NCBI.tsv          (rescate de huérfanos BOLD)
#   - tax_consensus.tsv     (tabla final consolidada)
#   - assignment_summary.tsv (% asignación por rango y DB)
# =============================================================================

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(yaml)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(rprojroot)
})

ROOT   <- find_root(has_dir("RUN29"))
params <- yaml::read_yaml(file.path(ROOT, "config", "params.yml"))

dada2_dir <- file.path(ROOT, params$paths$dada2_out)
tax_dir   <- file.path(ROOT, params$paths$taxonomy_out)
dir.create(tax_dir, showWarnings = FALSE, recursive = TRUE)

asv_fa <- file.path(dada2_dir, "asvs.fasta")
stopifnot(file.exists(asv_fa))
asvs <- readDNAStringSet(asv_fa)
message("[04] ASVs cargados: ", length(asvs))

STRATEGY <- params$taxonomy$strategy
message("[04] Estrategia: ", STRATEGY)

# ---- Helpers ----------------------------------------------------------------
RANK_NAMES <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

parse_dada2_taxon <- function(x) {
  # "k__Eukaryota;p__Arthropoda;c__Insecta;..." o "Eukaryota;Arthropoda;..."
  if (is.na(x) || x == "Unassigned") return(rep(NA_character_, 7))
  parts <- strsplit(x, ";", fixed = TRUE)[[1]]
  parts <- gsub("^[a-z]__", "", trimws(parts))
  length(parts) <- 7
  parts
}

# ============================================================================
# VSEARCH consensus (estrategia por defecto)
# ============================================================================
run_vsearch <- function(asv_fa, ref_fa, out_tsv,
                        pident = 0.85, maxaccepts = 10, threads = 8) {
  if (!file.exists(ref_fa)) {
    stop("Referencia no encontrada: ", ref_fa,
         "\nEjecuta primero workflow/04a o 04b para construir la DB.")
  }
  if (Sys.which("vsearch") == "") {
    stop("vsearch no encontrado. Instala: mamba install -c bioconda vsearch")
  }
  cmd <- paste(
    "vsearch --usearch_global", shQuote(asv_fa),
    "--db", shQuote(ref_fa),
    "--id", pident,
    "--top_hits_only --maxaccepts", maxaccepts,
    "--strand both --threads", threads,
    "--blast6out", shQuote(out_tsv)
  )
  message("[04] $ ", cmd)
  res <- system(cmd, intern = FALSE)
  if (res != 0) stop("vsearch falló con código ", res)
}

# Consenso a nivel de género/familia desde múltiples top hits
collapse_consensus <- function(blast6_tsv, tax_tsv, min_consensus = 0.51) {
  if (!file.exists(blast6_tsv) || file.info(blast6_tsv)$size == 0) {
    return(tibble())
  }
  hits <- read_tsv(blast6_tsv, show_col_types = FALSE,
                   col_names = c("qseqid","sseqid","pident","length",
                                 "mismatch","gapopen","qstart","qend",
                                 "sstart","send","evalue","bitscore"))
  taxmap <- read_tsv(tax_tsv, show_col_types = FALSE) |>
    rename(sseqid = `Feature ID`, Taxon = Taxon)

  hits |>
    inner_join(taxmap, by = "sseqid") |>
    group_by(qseqid) |>
    summarise(
      n_hits   = n(),
      best_pid = max(pident),
      Taxon = {
        ranks_split <- lapply(Taxon, parse_dada2_taxon)
        consensus_per_rank <- vapply(seq_along(RANK_NAMES), function(i) {
          vals <- vapply(ranks_split, `[`, character(1), i)
          vals <- vals[!is.na(vals)]
          if (length(vals) == 0) return(NA_character_)
          tb <- sort(table(vals), decreasing = TRUE)
          if (tb[1] / sum(tb) >= min_consensus) names(tb)[1] else NA_character_
        }, character(1))
        # Trunca al primer NA (consensus colapsa a nivel superior)
        first_na <- which(is.na(consensus_per_rank))[1]
        if (!is.na(first_na)) consensus_per_rank[first_na:7] <- NA_character_
        paste(consensus_per_rank, collapse = ";")
      },
      .groups = "drop"
    ) |>
    separate(Taxon, into = RANK_NAMES, sep = ";", fill = "right")
}

# ============================================================================
# Estrategia: vsearch_consensus
# ============================================================================
if (STRATEGY == "vsearch_consensus") {

  # --- 1) BOLD (primario) ---------------------------------------------------
  bold_fa  <- file.path(ROOT, params$taxonomy$bold$seqs_fasta)
  bold_tax <- file.path(ROOT, params$taxonomy$bold$tax_tsv)
  bold_out <- file.path(tax_dir, "vsearch_BOLD.b6")

  message("[04] VSEARCH ASVs vs BOLD...")
  run_vsearch(asv_fa, bold_fa, bold_out,
              pident     = params$taxonomy$classification$perc_identity,
              maxaccepts = params$taxonomy$classification$maxaccepts,
              threads    = params$taxonomy$classification$threads)

  bold_tax_df <- collapse_consensus(bold_out, bold_tax,
                                    min_consensus = params$taxonomy$classification$min_consensus)
  if (nrow(bold_tax_df) > 0) {
    write_tsv(bold_tax_df, file.path(tax_dir, "tax_BOLD.tsv"))
  } else {
    write_tsv(tibble(qseqid = character()), file.path(tax_dir, "tax_BOLD.tsv"))
  }
  message("[04] BOLD asignó: ", nrow(bold_tax_df), "/", length(asvs), " ASVs")

  # --- 2) NCBI (rescate de huérfanos) --------------------------------------
  asv_ids_all <- names(asvs)
  resolved_by_bold_genus <- bold_tax_df |>
    filter(!is.na(Genus)) |>
    pull(qseqid)
  orphans <- setdiff(asv_ids_all, resolved_by_bold_genus)
  message("[04] ASVs huérfanos de BOLD (sin Género) → NCBI: ", length(orphans))

  if (length(orphans) > 0) {
    orphan_fa <- file.path(tax_dir, "asvs_orphans.fasta")
    writeXStringSet(asvs[orphans], orphan_fa)

    ncbi_fa  <- file.path(ROOT, params$taxonomy$ncbi$seqs_fasta)
    ncbi_tax <- file.path(ROOT, params$taxonomy$ncbi$tax_tsv)
    ncbi_out <- file.path(tax_dir, "vsearch_NCBI.b6")

    message("[04] VSEARCH huérfanos vs NCBI...")
    run_vsearch(orphan_fa, ncbi_fa, ncbi_out,
                pident     = params$taxonomy$classification$perc_identity,
                maxaccepts = params$taxonomy$classification$maxaccepts,
                threads    = params$taxonomy$classification$threads)
    ncbi_tax_df <- collapse_consensus(ncbi_out, ncbi_tax,
                                      min_consensus = params$taxonomy$classification$min_consensus)
    if (nrow(ncbi_tax_df) > 0) {
      write_tsv(ncbi_tax_df, file.path(tax_dir, "tax_NCBI.tsv"))
    } else {
      ncbi_tax_df <- tibble(qseqid = character())
      write_tsv(ncbi_tax_df, file.path(tax_dir, "tax_NCBI.tsv"))
    }
    message("[04] NCBI rescató: ", nrow(ncbi_tax_df), "/", length(orphans))
  } else {
    ncbi_tax_df <- tibble(qseqid = character())
  }

  # --- 3) Consenso final ---------------------------------------------------
  # BOLD tiene prioridad cuando ambas DBs asignan; NCBI complementa
  bold_assigned <- bold_tax_df |>
    select(qseqid, all_of(RANK_NAMES), best_pid_bold = best_pid, n_hits_bold = n_hits)
  ncbi_assigned <- ncbi_tax_df |>
    rename(best_pid_ncbi = best_pid, n_hits_ncbi = n_hits) |>
    select(qseqid, all_of(RANK_NAMES), best_pid_ncbi, n_hits_ncbi)

  consensus <- tibble(qseqid = names(asvs)) |>
    left_join(bold_assigned, by = "qseqid", suffix = c("", ".bold")) |>
    left_join(ncbi_assigned, by = "qseqid", suffix = c(".bold", ".ncbi"))

  # Combinar: para cada rango, BOLD primero; si NA usar NCBI
  for (r in RANK_NAMES) {
    bold_col <- paste0(r, ".bold"); ncbi_col <- paste0(r, ".ncbi")
    if (bold_col %in% names(consensus) && ncbi_col %in% names(consensus)) {
      consensus[[r]] <- coalesce(consensus[[bold_col]], consensus[[ncbi_col]])
    }
  }
  consensus <- consensus |>
    mutate(source = case_when(
      !is.na(Kingdom) & !is.na(Genus.bold) ~ "BOLD",
      !is.na(Kingdom) & !is.na(Genus.ncbi) ~ "NCBI",
      !is.na(Kingdom)                       ~ "BOLD_partial",
      TRUE                                   ~ "Unassigned"
    )) |>
    select(qseqid, all_of(RANK_NAMES),
           best_pid_bold, n_hits_bold, best_pid_ncbi, n_hits_ncbi, source)

  write_tsv(consensus, file.path(tax_dir, "tax_consensus.tsv"))

  # --- 4) Resumen ----------------------------------------------------------
  summary_tbl <- tibble(
    rank = RANK_NAMES,
    assigned = vapply(RANK_NAMES, function(r) sum(!is.na(consensus[[r]])), integer(1)),
    pct      = round(100 * vapply(RANK_NAMES, function(r) sum(!is.na(consensus[[r]])), integer(1)) /
                       nrow(consensus), 1)
  )
  write_tsv(summary_tbl, file.path(tax_dir, "assignment_summary.tsv"))

  cat("\n=== Resumen 04_taxonomy (vsearch_consensus) ===\n")
  print(summary_tbl)
  cat("\nFuente de asignación:\n")
  print(table(consensus$source))
}

# ============================================================================
# Estrategia: naive_bayes (delega a QIIME2 vía system())
# ============================================================================
if (STRATEGY == "naive_bayes") {
  cat("\n=== Estrategia: naive_bayes (QIIME2 classify-sklearn) ===\n")
  cat("Pasos en shell (ejecutar fuera de R con qiime2-rescript activo):\n\n")
  cat("# 1) Importar ASVs como QZA\n")
  cat("qiime tools import \\\n")
  cat("  --input-path", asv_fa, "\\\n")
  cat("  --output-path", file.path(tax_dir, "asvs.qza"),
      "\\\n  --type 'FeatureData[Sequence]'\n\n")
  cat("# 2) Clasificar con BOLD\n")
  cat("qiime feature-classifier classify-sklearn \\\n")
  cat("  --i-classifier", file.path(ROOT, params$taxonomy$bold$classifier_qza), "\\\n")
  cat("  --i-reads",      file.path(tax_dir, "asvs.qza"), "\\\n")
  cat("  --o-classification", file.path(tax_dir, "tax_BOLD_nb.qza"), "\\\n")
  cat("  --p-n-jobs", params$taxonomy$classification$threads, "\n\n")
  cat("# 3) Exportar a TSV: qiime tools export ...\n")
}

# ============================================================================
# Estrategia: dada2_R (assignTaxonomy() puro)
# ============================================================================
if (STRATEGY == "dada2_R") {
  bold_fa <- file.path(ROOT, params$taxonomy$bold$seqs_fasta)
  message("[04] assignTaxonomy() vs BOLD FASTA (DADA2-format)...")
  tax_bold <- assignTaxonomy(
    seqs        = as.character(asvs),
    refFasta    = bold_fa,
    minBoot     = params$taxonomy$classification$dada2_minBoot,
    tryRC       = params$taxonomy$classification$dada2_tryRC,
    multithread = TRUE,
    verbose     = TRUE
  )
  out <- as_tibble(tax_bold, rownames = "seq") |>
    mutate(qseqid = names(asvs)[match(seq, as.character(asvs))]) |>
    select(qseqid, everything(), -seq)
  write_tsv(out, file.path(tax_dir, "tax_BOLD_dada2.tsv"))
  message("Salida: ", file.path(tax_dir, "tax_BOLD_dada2.tsv"))
}

cat("\nAcciones recomendadas:\n")
cat(" - Revisar tax_consensus.tsv y assignment_summary.tsv.\n")
cat(" - Asignar a phyloseq: tax_consensus + asv_table + samples_RUN29_COI.tsv.\n")
cat(" - Comparar pct asignación BOLD vs NCBI (material de discusión del workshop).\n")
