# =============================================================================
# 03_dada2_pipeline_SE.R   (Solución 2 — single-end forward-only)
# DADA2 para COI Leray-XT (LANCO/RUN29) — solo R1, 38 muestras L
# -----------------------------------------------------------------------------
# Estrategia SINGLE-END:
#   - R2 MiSeq 2×301 colapsa Q<20 después del ciclo ~180; PE+justConcatenate
#     retuvo ~30 % (pipeline _PE.R archivado).
#   - R1 mantiene calidad alta hasta pos ~220 → usamos solo R1.
#   - ASV resultante (~220 nt) cubre los 4 dominios variables 5' del amplicón
#     Leray, suficiente para asignación a género/familia.
#   - Endorse: Callahan, DADA2 GitHub #176; Antich et al. 2021 (BMC Bioinf.).
#
# Trade-off frente al PE+concat anterior:
#   + Retención esperada 70–90 % (vs 30 %).
#   + ASV biológicamente continuo → NUMT-screen por traducción REACTIVO.
#   - Pierde ~100 nt 3' del amplicón → menor poder resolutivo intra-especie.
# -----------------------------------------------------------------------------
# Entrada:  results/02_cutadapt/<sample>_R1.fastq.gz  (R2 IGNORADO)
# Salidas en results/03_dada2_SE/:
#   - filtered/<sample>_R1.fastq.gz
#   - quality_postfilter_R1.pdf       (diagnóstico post-Q20)
#   - errors_R1.rds, error_model_R1.pdf
#   - dada_R1.rds
#   - seqtab.rds, seqtab_nochim.rds
#   - filter_trim_stats.tsv
#   - track_reads.tsv                 (pérdidas por etapa)
#   - asvs.fasta, asv_table.tsv       (IDs estables ASV_0001 ...)
# =============================================================================

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(yaml)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(rprojroot)
})

ROOT   <- find_root(has_dir("RUN29"))
params <- yaml::read_yaml(file.path(ROOT, "config", "params.yml"))

# Coerción defensiva: yaml lee 1e-40 como character en algunas versiones.
params$dada2_inference$OMEGA_A <- as.numeric(params$dada2_inference$OMEGA_A)
stopifnot(is.finite(params$dada2_inference$OMEGA_A),
          params$dada2_inference$OMEGA_A > 0,
          params$dada2_inference$OMEGA_A < 1)

samples <- read_tsv(file.path(ROOT, params$paths$metadata), show_col_types = FALSE)
stopifnot(all(samples$target_detected == "COI_LerayXT"))
message("[03-SE] Procesando ", nrow(samples), " muestras COI Leray-XT (single-end R1)")

sids     <- samples$sample_id
cut_dir  <- file.path(ROOT, params$paths$cutadapt_out)

# Salidas SEPARADAS del pipeline PE para no pisar resultados previos.
out_dir  <- file.path(ROOT, "results", "03_dada2_SE")
filt_dir <- file.path(out_dir, "filtered")
dir.create(filt_dir, showWarnings = FALSE, recursive = TRUE)

fnFs   <- file.path(cut_dir,  paste0(sids, "_R1.fastq.gz"))
filtFs <- file.path(filt_dir, paste0(sids, "_R1.fastq.gz"))
names(fnFs) <- names(filtFs) <- sids

missing <- !file.exists(fnFs)
if (any(missing)) {
  stop("Faltan archivos cutadapt R1 para: ", paste(sids[missing], collapse = ", "),
       "\nEjecuta primero workflow/02_cutadapt_primers.sh")
}

# ---- 1) filterAndTrim single-end --------------------------------------------
ft <- params$filter_and_trim
message(sprintf("[03-SE] filterAndTrim R1-only  (truncQ=%d, maxEE=%d, truncLen=%d)",
                ft$truncQ, ft$maxEE_R1, ft$truncLen_R1))

out_ft <- filterAndTrim(
  fwd          = fnFs,  filt = filtFs,
  truncLen     = ft$truncLen_R1,
  maxN         = ft$maxN,
  maxEE        = ft$maxEE_R1,
  truncQ       = ft$truncQ,
  rm.phix      = ft$rm_phix,
  compress     = TRUE,
  multithread  = ft$multithread,
  verbose      = TRUE
)
saveRDS(out_ft, file.path(out_dir, "filter_trim_stats.rds"))

ft_stats <- as_tibble(out_ft, rownames = "file") |>
  mutate(sample_id = sub("_R1\\.fastq\\.gz$", "", file),
         pct_kept  = round(100 * reads.out / reads.in, 1))
write_tsv(ft_stats, file.path(out_dir, "filter_trim_stats.tsv"))

flagged <- ft_stats |> filter(pct_kept < 60)
if (nrow(flagged) > 0) {
  warning("Muestras con <60 % retenido tras filter SE: ",
          paste(flagged$sample_id, collapse = ", "),
          "\n→ Acción recomendada: inspeccionar muestras individuales o subir truncLen.")
}

# Descartar muestras vacías post-filtro
keep   <- file.exists(filtFs) & file.info(filtFs)$size > 50
filtFs <- filtFs[keep]; sids <- sids[keep]
message("[03-SE] Muestras tras filter (no vacías): ", length(sids))
stopifnot(length(sids) > 0)

# ---- 2) Quality profile POST-filter ----------------------------------------
message("[03-SE] Quality profile post-filterAndTrim R1...")
set.seed(42)
pick <- sample(sids, min(12, length(sids)))
qp_R1 <- plotQualityProfile(filtFs[pick]) +
  ggtitle("Forward (R1) — post-filterAndTrim single-end")
ggsave(file.path(out_dir, "quality_postfilter_R1.pdf"), qp_R1, width = 12, height = 8)

# ---- 3) learnErrors --------------------------------------------------------
message("[03-SE] learnErrors R1...")
errF <- learnErrors(filtFs,
                    multithread = params$dada2_inference$multithread,
                    verbose = TRUE)
saveRDS(errF, file.path(out_dir, "errors_R1.rds"))
ggsave(file.path(out_dir, "error_model_R1.pdf"),
       plotErrors(errF, nominalQ = TRUE), width = 10, height = 8)

# ---- 4) dada() inferencia --------------------------------------------------
message("[03-SE] dada() R1 (pool=", params$dada2_inference$pool, ")...")
dadaFs <- dada(filtFs, err = errF,
               pool        = params$dada2_inference$pool,
               OMEGA_A     = params$dada2_inference$OMEGA_A,
               multithread = params$dada2_inference$multithread)
saveRDS(dadaFs, file.path(out_dir, "dada_R1.rds"))

# ---- 5) Sequence table (sin mergePairs en single-end) ----------------------
seqtab <- makeSequenceTable(dadaFs)
message("[03-SE] seqtab raw:        ", nrow(seqtab), " muestras x ", ncol(seqtab), " ASVs")
saveRDS(seqtab, file.path(out_dir, "seqtab.rds"))

seqtab_nc <- removeBimeraDenovo(
  seqtab, method = params$chimera$method,
  multithread = TRUE, verbose = TRUE
)
message(sprintf("[03-SE] seqtab post-chimera: %d ASVs (%.1f%% reads conservados)",
                ncol(seqtab_nc), 100 * sum(seqtab_nc) / sum(seqtab)))

# ---- 6) Filtro de longitud COI Leray (R1 only) -----------------------------
lens   <- nchar(colnames(seqtab_nc))
keep_l <- lens >= params$coi_filtering$min_length & lens <= params$coi_filtering$max_length
message(sprintf("[03-SE] Filtro longitud (%d-%d nt): %d/%d ASVs conservados",
                params$coi_filtering$min_length, params$coi_filtering$max_length,
                sum(keep_l), length(lens)))
seqtab_nc <- seqtab_nc[, keep_l, drop = FALSE]

# ---- 7) NUMT screen por traducción (ACTIVO en single-end) ------------------
# Justificación: en SE el ASV es R1 biológicamente continuo, sin linker N que
# rompa el marco de lectura → translate() es válido. Código genético 5
# (NCBI table 5, Invertebrate Mitochondrial). NUMTs típicamente generan stop
# codons internos → se descartan.
if (isTRUE(params$coi_filtering$apply_translation_check)) {
  message("[03-SE] NUMT screen (código genético ",
          params$coi_filtering$genetic_code, ")...")
  gc <- getGeneticCode(as.character(params$coi_filtering$genetic_code))
  asv_dna <- DNAStringSet(colnames(seqtab_nc))
  has_orf <- vapply(asv_dna, function(s) {
    any(vapply(0:2, function(off) {
      n_trim <- ((nchar(s) - off) %/% 3) * 3
      if (n_trim < 60) return(FALSE)
      aa <- translate(subseq(s, start = 1 + off, end = off + n_trim),
                      genetic.code = gc, if.fuzzy.codon = "solve")
      !grepl("\\*", as.character(aa))
    }, logical(1)))
  }, logical(1))
  message("[03-SE] ASVs con ORF abierto (no NUMT-like): ",
          sum(has_orf), "/", length(has_orf))
  seqtab_nc <- seqtab_nc[, has_orf, drop = FALSE]
}

saveRDS(seqtab_nc, file.path(out_dir, "seqtab_nochim.rds"))

# ---- 8) Exportar ASVs.fasta con IDs estables -------------------------------
asv_seqs <- colnames(seqtab_nc)
asv_ids  <- sprintf("ASV_%04d", seq_along(asv_seqs))
names(asv_seqs) <- asv_ids
writeXStringSet(DNAStringSet(asv_seqs), file.path(out_dir, "asvs.fasta"))

asv_tab <- seqtab_nc
colnames(asv_tab) <- asv_ids
write_tsv(as_tibble(asv_tab, rownames = "sample_id"),
          file.path(out_dir, "asv_table.tsv"))

# ---- 9) Tabla de seguimiento de reads --------------------------------------
getN <- function(x) sum(getUniques(x))
basenames_ft <- sub("_R1\\.fastq\\.gz$", "", rownames(out_ft))

track <- tibble(
  sample_id    = sids,
  input        = out_ft[match(sids, basenames_ft), "reads.in"],
  filtered_R1  = out_ft[match(sids, basenames_ft), "reads.out"],
  denoised_R1  = vapply(dadaFs[sids], getN, integer(1)),
  nonchim_filt = rowSums(seqtab_nc[sids, , drop = FALSE])
) |>
  mutate(
    pct_filtered = round(100 * filtered_R1  / input, 1),
    pct_retained = round(100 * nonchim_filt / input, 1)
  )

write_tsv(track, file.path(out_dir, "track_reads.tsv"))

cat("\n=== Resumen 03_dada2_pipeline_SE ===\n")
cat("Estrategia:                    SINGLE-END R1 only (truncLen=",
    ft$truncLen_R1, ")\n", sep="")
cat("Muestras procesadas:           ", nrow(track), "\n")
cat("ASVs finales (post-NUMT):      ", ncol(seqtab_nc), "\n")
cat("Reads totales conservados:     ", sum(seqtab_nc), "\n")
cat("Mediana % tras filterAndTrim:  ", median(track$pct_filtered),  "%\n")
cat("Mediana % retenido end-to-end: ", median(track$pct_retained),  "%\n")
cat("Outputs en:                    ", out_dir, "\n")
cat("\nAcciones recomendadas:\n")
cat(" - Revisar quality_postfilter_R1.pdf: confirmar Q20+ sostenido hasta truncLen.\n")
cat(" - Si pct_filtered < 60 %% en >5 muestras: revisar muestras individuales.\n")
cat(" - Comparar ncol(asvs) y diversidad vs pipeline _PE archivado.\n")
cat(" - Continuar:  workflow/04_taxonomy_BOLD_NCBI.R sobre", out_dir, "/asvs.fasta\n")
