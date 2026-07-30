# =============================================================================
# 06a_classify_OTUs_BOLD.R
# Parser + consensus taxonomy sobre resultados VSEARCH usearch_global vs BOLD.
# LANCO/RUN29 COI Leray-XT — post swarm+LULU (paso 05).
# -----------------------------------------------------------------------------
# Se ejecuta DESPUÉS del vsearch (invocado por 06a_classify_OTUs_BOLD.slurm).
# Input:
#   - results/06_taxonomy/bold/otus_vs_bold.b6         (VSEARCH blast6out)
#   - results/05_swarm_lulu/lulu/curated_otus.fasta    (secuencias OTUs)
#   - db/rescript/bold-coi-tax.tsv                     (mapping seqid → taxonomy)
#   - config/params.yml
# Outputs en results/06_taxonomy/bold/:
#   - otus_bold_assigned.tsv        (otu_id, source=BOLD, ranks..., perc_id)
#   - otus_bold_unassigned.tsv      (otu_id + criterio de descarte)
#   - otus_unassigned_BOLD.fasta    (input para 06b_..._NCBI)
#   - bold_classification_stats.tsv (resumen métrico)
# -----------------------------------------------------------------------------
# Consensus algo (por OTU, para cada rango en config.consensus.ranks):
#   1) Toma los top-N hits (post --top_hits_only) del OTU.
#   2) Aplica umbral min_perc_id_<rank> para descartar hits confiables solo
#      hasta rangos altos (ej: hit al 88% no puede resolver a species).
#   3) Calcula la fracción del hit-más-frecuente por rango.
#   4) Si fracción >= min_consensus (0.51), asigna ese taxon; si no,
#      unassigned para ese rango y todos los descendientes.
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
dir.create(bold_dir, recursive = TRUE, showWarnings = FALSE)

blast6_in <- file.path(bold_dir, "otus_vs_bold.b6")
otus_fa   <- file.path(ROOT, params$paths$swarm_lulu_out, "lulu", "curated_otus.fasta")
bold_tax  <- file.path(ROOT, params$taxonomy$bold$tax_tsv)
stopifnot(file.exists(blast6_in), file.exists(otus_fa), file.exists(bold_tax))

# ---- 1) Cargar hits VSEARCH + secuencias OTU + tax de referencia -----------
b6_cols <- c("query","subject","pct_id","aln_len","mismatch","gapopen",
             "qstart","qend","sstart","send","evalue","bitscore")
hits <- read_tsv(blast6_in, col_names = b6_cols, show_col_types = FALSE)
message(sprintf("[06a] VSEARCH hits leídos: %d (%d queries únicos)",
                nrow(hits), length(unique(hits$query))))

otu_seqs <- readDNAStringSet(otus_fa)
all_otus <- names(otu_seqs)

# BOLD tax.tsv formato QIIME2/RESCRIPt: Feature ID <TAB> Taxon
# Taxon = "k__X;p__Y;c__Z;o__A;f__B;g__C;s__D"
bold_tx <- read_tsv(bold_tax, col_names = c("subject","Taxon"),
                    skip = 1, show_col_types = FALSE)

parse_taxon <- function(tx_string, ranks) {
  # Split "k__X;p__Y;..." → named vector por rank
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

bold_tx_parsed <- bold_tx |>
  mutate(tx_list = map(Taxon, parse_taxon, ranks = cls$consensus$ranks)) |>
  unnest_wider(tx_list) |>
  select(-Taxon)

# ---- 2) Merge hits ↔ taxonomy + filtros por umbral perc_id por rango -------
hits_tx <- hits |>
  inner_join(bold_tx_parsed, by = "subject") |>
  select(query, subject, pct_id, all_of(cls$consensus$ranks))

# Colocar NA en el rango cuando el pct_id no alcanza el umbral de confianza
rank_thresholds <- c(
  species     = cls$consensus$min_perc_id_species,
  genus       = cls$consensus$min_perc_id_genus,
  family      = cls$consensus$min_perc_id_family,
  order       = cls$consensus$min_perc_id_order
)
for (rk in names(rank_thresholds)) {
  if (rk %in% names(hits_tx)) {
    hits_tx[[rk]] <- ifelse(hits_tx$pct_id < rank_thresholds[rk], NA_character_,
                             hits_tx[[rk]])
  }
}

# ---- 3) Consensus por OTU ---------------------------------------------------
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
    top_frac <- tab[1] / n_hits
    if (top_frac >= min_cons) {
      out[[rk]] <- names(tab)[1]
    } else {
      out[[rk]] <- NA_character_
      passed_upper <- FALSE   # sin consenso aquí → hijos también NA
    }
  }
  out
}

message("[06a] Computando consensus por OTU...")
# Skeleton garantiza columnas correctas incluso si hits_tx = 0 filas
consensus_skel <- tibble(
  otu_id = character(), n_hits = integer(), max_pct_id = double()
)
for (rk in cls$consensus$ranks) consensus_skel[[rk]] <- character()

if (nrow(hits_tx) == 0) {
  consensus_tbl <- consensus_skel
} else {
  consensus_tbl <- hits_tx |>
    group_by(query) |>
    group_split() |>
    map_dfr(consensus_row,
            ranks    = cls$consensus$ranks,
            min_cons = cls$consensus$min_consensus)
  if (!"otu_id" %in% names(consensus_tbl)) consensus_tbl <- consensus_skel
}

# OTUs sin hit alguno en BOLD
otus_no_hits <- setdiff(all_otus, consensus_tbl$otu_id %||% character())
if (length(otus_no_hits) > 0) {
  no_hit_rows <- tibble(
    otu_id = otus_no_hits, n_hits = 0L, max_pct_id = NA_real_
  )
  for (rk in cls$consensus$ranks) no_hit_rows[[rk]] <- NA_character_
  consensus_tbl <- bind_rows(consensus_tbl, no_hit_rows)
}

# ---- 4) Clasificar asignados vs unassigned ---------------------------------
uc <- cls$unassigned_criteria
rank_order <- cls$consensus$ranks
rank_depth <- function(row) {
  # profundidad = rango más específico con asignación (0 = ninguna, 7 = species)
  sum(!is.na(row[rank_order]))
}
depth <- vapply(seq_len(nrow(consensus_tbl)),
                function(i) rank_depth(consensus_tbl[i, ]), integer(1))
max_depth_ok <- match(uc$max_rank_assigned, rank_order)

is_assigned <- (
  !is.na(consensus_tbl$max_pct_id) &
  consensus_tbl$max_pct_id >= uc$min_perc_id &
  depth > max_depth_ok           # más profundo que "class"
)

assigned <- consensus_tbl[is_assigned, ] |>
  mutate(source = "BOLD") |>
  select(otu_id, source, max_pct_id, n_hits, all_of(cls$consensus$ranks))

unassigned <- consensus_tbl[!is_assigned, ] |>
  mutate(reason = case_when(
    is.na(max_pct_id)                   ~ "no_hits",
    max_pct_id < uc$min_perc_id         ~ sprintf("low_id_%.1f", max_pct_id),
    TRUE                                ~ sprintf("shallow_rank_depth_%d", depth[!is_assigned])
  )) |>
  select(otu_id, reason, max_pct_id, n_hits, all_of(cls$consensus$ranks))

# ---- 5) Escribir outputs ---------------------------------------------------
write_tsv(assigned,   file.path(bold_dir, "otus_bold_assigned.tsv"),
          na = "")
write_tsv(unassigned, file.path(bold_dir, "otus_bold_unassigned.tsv"),
          na = "")

# FASTA con los OTUs unassigned → input de 06b
unassigned_seqs <- otu_seqs[unassigned$otu_id]
writeXStringSet(unassigned_seqs,
                file.path(bold_dir, "otus_unassigned_BOLD.fasta"))

# Stats
stats <- tibble(
  total_otus         = length(all_otus),
  assigned_bold      = nrow(assigned),
  unassigned_bold    = nrow(unassigned),
  no_hits            = sum(unassigned$reason == "no_hits"),
  low_identity       = sum(startsWith(unassigned$reason, "low_id_")),
  shallow_taxonomy   = sum(startsWith(unassigned$reason, "shallow_rank_"))
)
write_tsv(stats, file.path(bold_dir, "bold_classification_stats.tsv"))

cat("\n=== Resumen 06a_classify_OTUs_BOLD ===\n")
print(as.data.frame(stats), row.names = FALSE)
cat(sprintf("\nAsignados: %d / %d (%.1f%%)\n",
            nrow(assigned), length(all_otus),
            100 * nrow(assigned) / length(all_otus)))
cat("\nOutputs en:", bold_dir, "\n")
cat(" → Continuar con: workflow/06b_classify_OTUs_NCBI.slurm\n")
cat(" → Input NCBI:", file.path(bold_dir, "otus_unassigned_BOLD.fasta"), "\n")
