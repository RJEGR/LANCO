# =============================================================================
# 06b_classify_OTUs_NCBI.R
# Parser + consensus taxonomy sobre VSEARCH vs NCBI COI,
# aplicado SOLO a OTUs no asignados en 06a (BOLD).
# Merge final con las asignaciones BOLD → tax_consensus.tsv
# LANCO/RUN29 COI Leray-XT — post 06a.
# -----------------------------------------------------------------------------
# Input:
#   - results/06_taxonomy/ncbi/otus_vs_ncbi.b6         (VSEARCH blast6out)
#   - results/06_taxonomy/bold/otus_bold_assigned.tsv  (asignadas fase 1)
#   - results/06_taxonomy/bold/otus_unassigned_BOLD.fasta
#   - db/rescript/ncbi-coi-tax.tsv
# Outputs en results/06_taxonomy/:
#   - ncbi/otus_ncbi_assigned.tsv
#   - ncbi/otus_ncbi_unassigned.tsv
#   - ncbi/ncbi_classification_stats.tsv
#   - tax_consensus.tsv          (BOLD + NCBI + unassigned, source-tagged)
#   - final_taxonomy_stats.tsv   (resumen global)
# =============================================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(yaml)
  library(rprojroot)
  library(purrr)
})

ROOT   <- find_root(has_dir("."))
params <- yaml::read_yaml(file.path(ROOT, "config", "params.yml"))
cls    <- params$taxonomy$classification

# ---- Rutas ------------------------------------------------------------------
tax_dir  <- file.path(ROOT, params$paths$taxonomy_out)
bold_dir <- file.path(tax_dir, "bold")
ncbi_dir <- file.path(tax_dir, "ncbi")
dir.create(ncbi_dir, recursive = TRUE, showWarnings = FALSE)

blast6_in <- file.path(ncbi_dir, "otus_vs_ncbi.b6")
bold_ass  <- file.path(bold_dir, "otus_bold_assigned.tsv")
ncbi_tax  <- file.path(ROOT, params$taxonomy$ncbi$tax_tsv)
stopifnot(file.exists(blast6_in), file.exists(bold_ass), file.exists(ncbi_tax))

# ---- 1) Cargar hits NCBI + tax de referencia -------------------------------
b6_cols <- c("query","subject","pct_id","aln_len","mismatch","gapopen",
             "qstart","qend","sstart","send","evalue","bitscore")
hits <- read_tsv(blast6_in, col_names = b6_cols, show_col_types = FALSE)
message(sprintf("[06b] VSEARCH NCBI hits: %d (%d queries únicos)",
                nrow(hits), length(unique(hits$query))))

ncbi_tx <- read_tsv(ncbi_tax, col_names = c("subject","Taxon"),
                    skip = 1, show_col_types = FALSE)

parse_taxon <- function(tx_string, ranks) {
  parts <- strsplit(tx_string, ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  prefix_to_rank <- c(k = "superkingdom", p = "phylum", c = "class",
                      o = "order", f = "family", g = "genus", s = "species")
  out <- setNames(rep(NA_character_, length(ranks)), ranks)
  for (p in parts) {
    m <- regmatches(p, regexec("^([kpcofgs])__(.*)$", p))[[1]]
    if (length(m) < 3) next
    rk <- prefix_to_rank[m[2]]
    val <- m[3]
    if (rk %in% ranks && nzchar(val)) out[rk] <- val
  }
  out
}

ncbi_tx_parsed <- ncbi_tx |>
  mutate(tx_list = map(Taxon, parse_taxon, ranks = cls$consensus$ranks)) |>
  unnest_wider(tx_list) |>
  select(-Taxon)

# ---- 2) Merge + filtros perc_id por rango ----------------------------------
hits_tx <- hits |>
  inner_join(ncbi_tx_parsed, by = "subject") |>
  select(query, subject, pct_id, all_of(cls$consensus$ranks))

rank_thresholds <- c(
  species = cls$consensus$min_perc_id_species,
  genus   = cls$consensus$min_perc_id_genus,
  family  = cls$consensus$min_perc_id_family,
  order   = cls$consensus$min_perc_id_order
)
for (rk in names(rank_thresholds)) {
  if (rk %in% names(hits_tx)) {
    hits_tx[[rk]] <- ifelse(hits_tx$pct_id < rank_thresholds[rk], NA_character_,
                             hits_tx[[rk]])
  }
}

# ---- 3) Consensus por OTU (idéntico algo que 06a) --------------------------
consensus_row <- function(sub_df, ranks, min_cons) {
  n_hits <- nrow(sub_df)
  out <- tibble(otu_id = sub_df$query[1],
                n_hits = n_hits,
                max_pct_id = max(sub_df$pct_id, na.rm = TRUE))
  passed_upper <- TRUE
  for (rk in ranks) {
    if (!passed_upper) { out[[rk]] <- NA_character_; next }
    vals <- sub_df[[rk]]
    non_na <- vals[!is.na(vals)]
    if (length(non_na) == 0) {
      out[[rk]] <- NA_character_; passed_upper <- FALSE; next
    }
    tab <- sort(table(non_na), decreasing = TRUE)
    if (tab[1] / n_hits >= min_cons) {
      out[[rk]] <- names(tab)[1]
    } else {
      out[[rk]] <- NA_character_; passed_upper <- FALSE
    }
  }
  out
}

message("[06b] Computando consensus NCBI...")
# Skeleton vacío garantiza columnas correctas si hits_tx tiene 0 filas
ncbi_consensus_skel <- tibble(
  otu_id = character(), n_hits = integer(), max_pct_id = double()
)
for (rk in cls$consensus$ranks) ncbi_consensus_skel[[rk]] <- character()

if (nrow(hits_tx) == 0) {
  ncbi_consensus <- ncbi_consensus_skel
} else {
  ncbi_consensus <- hits_tx |>
    group_by(query) |>
    group_split() |>
    map_dfr(consensus_row,
            ranks    = cls$consensus$ranks,
            min_cons = cls$consensus$min_consensus)
  if (!"otu_id" %in% names(ncbi_consensus)) ncbi_consensus <- ncbi_consensus_skel
}

# OTUs de entrada a NCBI que no obtuvieron hits
unassigned_bold_ids <- names(readDNAStringSet(
  file.path(bold_dir, "otus_unassigned_BOLD.fasta")))
otus_no_hits_ncbi <- setdiff(unassigned_bold_ids,
                             ncbi_consensus$otu_id %||% character())
if (length(otus_no_hits_ncbi) > 0) {
  no_hit_rows <- tibble(otu_id = otus_no_hits_ncbi, n_hits = 0L,
                        max_pct_id = NA_real_)
  for (rk in cls$consensus$ranks) no_hit_rows[[rk]] <- NA_character_
  ncbi_consensus <- bind_rows(ncbi_consensus, no_hit_rows)
}

# ---- 4) Split assigned / unassigned NCBI -----------------------------------
uc <- cls$unassigned_criteria
rank_order <- cls$consensus$ranks
depth <- vapply(seq_len(nrow(ncbi_consensus)),
                function(i) sum(!is.na(ncbi_consensus[i, rank_order])),
                integer(1))
max_depth_ok <- match(uc$max_rank_assigned, rank_order)
is_assigned <- (
  !is.na(ncbi_consensus$max_pct_id) &
  ncbi_consensus$max_pct_id >= uc$min_perc_id &
  depth > max_depth_ok
)

ncbi_assigned <- ncbi_consensus[is_assigned, ] |>
  mutate(source = "NCBI") |>
  select(otu_id, source, max_pct_id, n_hits, all_of(cls$consensus$ranks))

ncbi_unassigned <- ncbi_consensus[!is_assigned, ] |>
  mutate(source = "unassigned",
         reason = case_when(
           is.na(max_pct_id)            ~ "no_hits_either_db",
           max_pct_id < uc$min_perc_id  ~ sprintf("low_id_%.1f", max_pct_id),
           TRUE                         ~ sprintf("shallow_rank_depth_%d",
                                                  depth[!is_assigned])
         )) |>
  select(otu_id, source, reason, max_pct_id, n_hits,
         all_of(cls$consensus$ranks))

write_tsv(ncbi_assigned,   file.path(ncbi_dir, "otus_ncbi_assigned.tsv"),
          na = "")
write_tsv(ncbi_unassigned, file.path(ncbi_dir, "otus_ncbi_unassigned.tsv"),
          na = "")

# ---- 5) Merge final: BOLD + NCBI + unassigned → tax_consensus.tsv ---------
# Forzar tipos: sin col_types explícito, read_tsv adivina character para
# columnas numéricas con NAs, rompiendo bind_rows con la tabla NCBI (double).
bold_col_types <- cols(
  otu_id     = col_character(),
  source     = col_character(),
  max_pct_id = col_double(),
  n_hits     = col_integer(),
  .default   = col_character()
)
bold_assigned <- read_tsv(bold_ass, col_types = bold_col_types,
                          na = c("", "NA"))

# Normalizar columnas (unassigned no tiene "source"/"reason" en el mismo layout)
core_cols <- c("otu_id","source","max_pct_id","n_hits", cls$consensus$ranks)
tax_final <- bind_rows(
  bold_assigned |> select(all_of(core_cols)),
  ncbi_assigned |> select(all_of(core_cols)),
  ncbi_unassigned |>
    mutate(source = "unassigned") |>
    select(all_of(core_cols))
) |>
  arrange(otu_id)

# Añadir columna confidence: rango más profundo asignado
depth_map <- function(row) {
  d <- sum(!is.na(row[cls$consensus$ranks]))
  if (d == 0) return("none")
  cls$consensus$ranks[d]
}
tax_final$deepest_rank <- vapply(seq_len(nrow(tax_final)),
                                 function(i) depth_map(tax_final[i, ]),
                                 character(1))

write_tsv(tax_final, file.path(tax_dir, "tax_consensus.tsv"), na = "")

# ---- 6) Stats globales -----------------------------------------------------
stats_global <- tibble(
  total_otus              = nrow(tax_final),
  assigned_bold           = nrow(bold_assigned),
  assigned_ncbi           = nrow(ncbi_assigned),
  unassigned              = nrow(ncbi_unassigned),
  reached_species         = sum(tax_final$deepest_rank == "species"),
  reached_genus_or_deeper = sum(tax_final$deepest_rank %in%
                                  c("genus","species")),
  reached_family_or_deeper = sum(tax_final$deepest_rank %in%
                                  c("family","genus","species"))
)
write_tsv(stats_global, file.path(tax_dir, "final_taxonomy_stats.tsv"))

cat("\n=== Resumen 06b_classify_OTUs_NCBI + Merge Final ===\n")
print(as.data.frame(stats_global), row.names = FALSE)
cat(sprintf("\nCobertura BOLD:   %d OTUs (%.1f%%)\n",
            nrow(bold_assigned),
            100 * nrow(bold_assigned) / nrow(tax_final)))
cat(sprintf("Rescate NCBI:     %d OTUs (%.1f%%)\n",
            nrow(ncbi_assigned),
            100 * nrow(ncbi_assigned) / nrow(tax_final)))
cat(sprintf("Sin asignar:      %d OTUs (%.1f%%)\n",
            nrow(ncbi_unassigned),
            100 * nrow(ncbi_unassigned) / nrow(tax_final)))
cat("\nOutput consolidado:", file.path(tax_dir, "tax_consensus.tsv"), "\n")
