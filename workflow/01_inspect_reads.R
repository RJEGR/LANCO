# =============================================================================
# 01_inspect_reads.R
# Inspección de calidad pre-cutadapt (LANCO/RUN29, COI Leray-XT)
# -----------------------------------------------------------------------------
# Objetivo (workshop):
#   1) Cargar muestras desde metadata/samples_RUN29.tsv
#   2) Conteo de reads por muestra
#   3) Perfiles de calidad agregados (forward y reverse)
#   4) Detección/cuantificación de primers en los reads — confirma si cutadapt
#      es necesario y con qué tasa de éxito esperar
# -----------------------------------------------------------------------------
# Salidas en results/01_quality/:
#   - read_counts_raw.tsv
#   - quality_profile_R1.pdf, quality_profile_R2.pdf
#   - primer_hits.tsv (fracción de reads con primer detectado por muestra)
# =============================================================================

suppressPackageStartupMessages({
  library(dada2)        # plotQualityProfile, getN
  library(ShortRead)    # FastqStreamer, readFastq
  library(Biostrings)   # DNAString, vmatchPattern
  library(yaml)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# ---- Config ------------------------------------------------------------------
ROOT   <- rprojroot::find_root(rprojroot::has_dir("RUN29"))
params <- yaml::read_yaml(file.path(ROOT, "config", "params.yml"))

samples <- read_tsv(file.path(ROOT, params$paths$metadata), show_col_types = FALSE)
fnFs    <- file.path(ROOT, samples$R1)
fnRs    <- file.path(ROOT, samples$R2)
names(fnFs) <- samples$sample_id
names(fnRs) <- samples$sample_id

out_dir <- file.path(ROOT, params$paths$results, "01_quality")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1) Conteo de reads por muestra -----------------------------------------
message("[01] Contando reads...")
read_counts <- tibble(
  sample_id = samples$sample_id,
  group     = samples$group,
  reads_R1  = vapply(fnFs, function(f) as.integer(ShortRead::countFastq(f)[1, "records"]),
                     integer(1)),
  reads_R2  = vapply(fnRs, function(f) as.integer(ShortRead::countFastq(f)[1, "records"]),
                     integer(1))
)
write_tsv(read_counts, file.path(out_dir, "read_counts_raw.tsv"))

# ---- 2) Perfiles de calidad agregados ---------------------------------------
# Muestreamos hasta 12 muestras por grupo para legibilidad del PDF.
set.seed(42)
pick <- samples |>
  group_by(group) |>
  slice_sample(n = min(12, n()), replace = FALSE) |>
  pull(sample_id)

message("[01] Generando quality profiles (n=", length(pick), " muestras)...")
pQ_R1 <- plotQualityProfile(fnFs[pick]) +
  ggtitle("Forward (R1) — pre-cutadapt")
pQ_R2 <- plotQualityProfile(fnRs[pick]) +
  ggtitle("Reverse (R2) — pre-cutadapt")

ggsave(file.path(out_dir, "quality_profile_R1.pdf"), pQ_R1, width = 12, height = 8)
ggsave(file.path(out_dir, "quality_profile_R2.pdf"), pQ_R2, width = 12, height = 8)

# ---- 3) Detección de primers Leray-XT ---------------------------------------
# Convierte IUPAC degenerada -> regex y cuenta hits en los primeros N reads.
iupac_to_regex <- function(s) {
  map <- c(A="A", C="C", G="G", T="T",
           R="[AG]", Y="[CT]", S="[GC]", W="[AT]", K="[GT]", M="[AC]",
           B="[CGT]", D="[AGT]", H="[ACT]", V="[ACG]", N="[ACGT]")
  paste0(map[strsplit(s, "")[[1]]], collapse = "")
}

fwd_rx <- paste0("^", iupac_to_regex(params$primers$fwd_seq))
rev_rx <- paste0("^", iupac_to_regex(params$primers$rev_seq))

count_primer_hits <- function(fastq_gz, regex, n = 5000) {
  fq <- ShortRead::readFastq(fastq_gz)
  fq <- fq[seq_len(min(n, length(fq)))]
  s  <- as.character(ShortRead::sread(fq))
  sum(grepl(regex, s))
}

message("[01] Detectando primers en primeros 5000 reads/muestra...")
primer_hits <- tibble(
  sample_id      = samples$sample_id,
  group          = samples$group,
  hits_fwd_R1    = vapply(fnFs, count_primer_hits, integer(1), regex = fwd_rx),
  hits_rev_R2    = vapply(fnRs, count_primer_hits, integer(1), regex = rev_rx)
) |>
  mutate(
    pct_fwd_R1 = round(100 * hits_fwd_R1 / 5000, 1),
    pct_rev_R2 = round(100 * hits_rev_R2 / 5000, 1)
  )
write_tsv(primer_hits, file.path(out_dir, "primer_hits.tsv"))

# ---- Resumen consola --------------------------------------------------------
cat("\n=== Resumen 01_inspect_reads ===\n")
cat("Muestras totales:           ", nrow(samples), "\n")
cat("Reads R1 (mediana):         ", median(read_counts$reads_R1), "\n")
cat("Primer fwd presente (mediana %): ", median(primer_hits$pct_fwd_R1), "%\n")
cat("Primer rev presente (mediana %): ", median(primer_hits$pct_rev_R2), "%\n")
cat("Salidas en:                 ", out_dir, "\n")
cat("\nAcciones recomendadas:\n")
cat(" - Revisar quality_profile_R{1,2}.pdf y ajustar truncLen en params.yml.\n")
cat(" - Si pct_fwd_R1 o pct_rev_R2 < 70 %% en alguna muestra → revisar primers/contaminación.\n")
cat(" - Proceder a workflow/02_cutadapt_primers.sh.\n")
