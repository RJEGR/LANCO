# =============================================================================
# 03_dada2_pipeline_PE.R  *** ARCHIVADO 2026-06-23 ***
# Versión PAIRED-END con mergePairs(justConcatenate=TRUE).
# Resultado real (RUN29, 38 muestras L): pct_retained mediana ~30 %,
# insuficiente para análisis aguas abajo.
# Conservado como snapshot reproducible — NO ejecutar en producción.
# Pipeline activo: workflow/03_dada2_pipeline_SE.R (single-end R1 only).
# =============================================================================
# (Contenido original sigue intacto desde aquí.)
# -----------------------------------------------------------------------------
# DADA2 para COI Leray-XT (LANCO/RUN29) — solo muestras L (n=38)
# -----------------------------------------------------------------------------
# Estrategia de quality trimming ESTRICTA:
#   - truncQ = 20  → cada read se trunca en la PRIMERA base con Q<20
#   - maxEE = (1, 2)  → tolerancia de error esperado endurecida
#   - truncLen ajustado al colapso real de calidad observado en muestras L:
#       R1 mantiene Q>=20 hasta ~pos 175 (raw); post-cutadapt: truncLen 220
#       R2 colapsa después de ~pos 180 (raw);     post-cutadapt: truncLen 150
#   - Sum truncLen = 370 ≥ 313 + 12 (overlap mínimo para mergePairs) ✓
# -----------------------------------------------------------------------------
# Entrada:  results/02_cutadapt/<sample>_R{1,2}.fastq.gz  (solo L1..L38)
# Salidas en results/03_dada2/:
#   - filtered/<sample>_R{1,2}.fastq.gz
#   - quality_postfilter_R{1,2}.pdf   (diagnóstico post-Q20)
#   - errors_R{1,2}.rds
#   - dada_R{1,2}.rds
#   - mergers.rds
#   - seqtab.rds, seqtab_nochim.rds
#   - track_reads.tsv                 (pérdidas por etapa)
#   - asvs.fasta, asv_table.tsv       (IDs estables ASV_0001 ...)
# -----------------------------------------------------------------------------
# COI-specific post-processing:
#   - Filtro de longitud 299–320 nt centrado en amplicón Leray (313 nt)
#   - Filtro por traducción (NUMT screen), código genético 5 (Invertebrate Mt)
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
params <- yaml::read_yaml(file.path(ROOT, "config", "params_PE.yml"))

# Coerción defensiva: yaml::read_yaml lee notación científica (1e-40) como
# character en algunas versiones, lo que rompe dada(OMEGA_A=). Forzamos numeric.
params$dada2_inference$OMEGA_A <- as.numeric(params$dada2_inference$OMEGA_A)
stopifnot(is.finite(params$dada2_inference$OMEGA_A),
          params$dada2_inference$OMEGA_A > 0,
          params$dada2_inference$OMEGA_A < 1)

samples <- read_tsv(file.path(ROOT, params$paths$metadata), show_col_types = FALSE)
stopifnot(all(samples$target_detected == "COI_LerayXT"))
message("[03] Procesando ", nrow(samples), " muestras COI Leray-XT (grupo L)")

sids     <- samples$sample_id
cut_dir  <- file.path(ROOT, params$paths$cutadapt_out)
filt_dir <- file.path(ROOT, params$paths$filtered_out)
out_dir  <- file.path(ROOT, params$paths$dada2_out)
dir.create(filt_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)

fnFs   <- file.path(cut_dir,  paste0(sids, "_R1.fastq.gz"))
fnRs   <- file.path(cut_dir,  paste0(sids, "_R2.fastq.gz"))
filtFs <- file.path(filt_dir, paste0(sids, "_R1.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(sids, "_R2.fastq.gz"))
names(fnFs) <- names(fnRs) <- names(filtFs) <- names(filtRs) <- sids

missing <- !file.exists(fnFs) | !file.exists(fnRs)
if (any(missing)) {
  stop("Faltan archivos cutadapt para: ", paste(sids[missing], collapse = ", "),
       "\nEjecuta primero workflow/02_cutadapt_primers.sh")
}

# ---- 1) filterAndTrim Q20 estricto ------------------------------------------
ft <- params$filter_and_trim
message(sprintf("[03] filterAndTrim Q20-strict  (truncQ=%d, maxEE=(%d,%d), truncLen=(%d,%d))",
                ft$truncQ, ft$maxEE_R1, ft$maxEE_R2, ft$truncLen_R1, ft$truncLen_R2))

out_ft <- filterAndTrim(
  fwd          = fnFs,  filt    = filtFs,
  rev          = fnRs,  filt.rev = filtRs,
  truncLen     = c(ft$truncLen_R1, ft$truncLen_R2),
  maxN         = ft$maxN,
  maxEE        = c(ft$maxEE_R1, ft$maxEE_R2),
  truncQ       = ft$truncQ,           # ← Q20 estricto
  rm.phix      = ft$rm_phix,
  compress     = TRUE,
  multithread  = ft$multithread,
  verbose      = TRUE
)
saveRDS(out_ft, file.path(out_dir, "filter_trim_stats.rds"))

# Diagnóstico: porcentaje retenido por muestra tras Q20
ft_stats <- as_tibble(out_ft, rownames = "file") |>
  mutate(sample_id = sub("_R1\\.fastq\\.gz$", "", file),
         pct_kept  = round(100 * reads.out / reads.in, 1))
write_tsv(ft_stats, file.path(out_dir, "filter_trim_stats.tsv"))

flagged <- ft_stats |> filter(pct_kept < 50)
if (nrow(flagged) > 0) {
  warning("Muestras con <50 % retenido tras Q20-strict: ",
          paste(flagged$sample_id, collapse = ", "),
          "\n→ Acción recomendada: revisar quality_postfilter_R*.pdf y considerar relajar truncQ.")
}

# Descartar muestras vacías post-filtro
keep   <- file.exists(filtFs) & file.info(filtFs)$size > 50
filtFs <- filtFs[keep]; filtRs <- filtRs[keep]; sids <- sids[keep]
message("[03] Muestras tras filter (no vacías): ", length(sids))

# ---- 2) Quality profile POST-filter (verifica que Q20 trimming funcionó) ----
message("[03] Quality profiles post-filterAndTrim...")
set.seed(42)
pick <- sample(sids, min(12, length(sids)))
qp_R1 <- plotQualityProfile(filtFs[pick]) +
  ggtitle("Forward (R1) — post-filterAndTrim Q20-strict")
qp_R2 <- plotQualityProfile(filtRs[pick]) +
  ggtitle("Reverse (R2) — post-filterAndTrim Q20-strict")
ggsave(file.path(out_dir, "quality_postfilter_R1.pdf"), qp_R1, width = 12, height = 8)
ggsave(file.path(out_dir, "quality_postfilter_R2.pdf"), qp_R2, width = 12, height = 8)

# ---- 3) learnErrors --------------------------------------------------------
message("[03] learnErrors R1...")
errF <- learnErrors(filtFs, multithread = params$dada2_inference$multithread, verbose = TRUE)
saveRDS(errF, file.path(out_dir, "errors_R1.rds"))

message("[03] learnErrors R2...")
errR <- learnErrors(filtRs, multithread = params$dada2_inference$multithread, verbose = TRUE)
saveRDS(errR, file.path(out_dir, "errors_R2.rds"))

ggsave(file.path(out_dir, "error_model_R1.pdf"), plotErrors(errF, nominalQ = TRUE), width = 10, height = 8)
ggsave(file.path(out_dir, "error_model_R2.pdf"), plotErrors(errR, nominalQ = TRUE), width = 10, height = 8)

# ---- 4) dada() inferencia --------------------------------------------------
message("[03] dada() forward (pool=", params$dada2_inference$pool, ")...")
dadaFs <- dada(filtFs, err = errF,
               pool        = params$dada2_inference$pool,
               OMEGA_A     = params$dada2_inference$OMEGA_A,
               multithread = params$dada2_inference$multithread)
saveRDS(dadaFs, file.path(out_dir, "dada_R1.rds"))

message("[03] dada() reverse...")
dadaRs <- dada(filtRs, err = errR,
               pool        = params$dada2_inference$pool,
               OMEGA_A     = params$dada2_inference$OMEGA_A,
               multithread = params$dada2_inference$multithread)
saveRDS(dadaRs, file.path(out_dir, "dada_R2.rds"))

# ---- 5) mergePairs ---------------------------------------------------------
# Con truncLen Q20-reality (140+160=300 nt) < Leray amplicón (313 nt), NO hay
# overlap suficiente. justConcatenate=TRUE concatena R1+10×N+revComp(R2) → 310 nt.
just_concat <- isTRUE(params$merge_pairs$justConcatenate)
message("[03] mergePairs (justConcatenate=", just_concat, ")...")
mergers <- mergePairs(
  dadaFs, filtFs, dadaRs, filtRs,
  minOverlap      = params$merge_pairs$minOverlap,
  maxMismatch     = params$merge_pairs$maxMismatch,
  justConcatenate = just_concat,
  verbose         = TRUE
)
saveRDS(mergers, file.path(out_dir, "mergers.rds"))

# ---- 6) Sequence table + chimera removal ----------------------------------
seqtab <- makeSequenceTable(mergers)
message("[03] seqtab raw:        ", nrow(seqtab), " muestras x ", ncol(seqtab), " ASVs")
saveRDS(seqtab, file.path(out_dir, "seqtab.rds"))

seqtab_nc <- removeBimeraDenovo(
  seqtab, method = params$chimera$method,
  multithread = TRUE, verbose = TRUE
)
message(sprintf("[03] seqtab post-chimera: %d ASVs (%.1f%% reads conservados)",
                ncol(seqtab_nc), 100 * sum(seqtab_nc) / sum(seqtab)))

# ---- 7) Filtro de longitud COI Leray ---------------------------------------
lens   <- nchar(colnames(seqtab_nc))
keep_l <- lens >= params$coi_filtering$min_length & lens <= params$coi_filtering$max_length
message(sprintf("[03] Filtro longitud (%d-%d nt): %d/%d ASVs conservados",
                params$coi_filtering$min_length, params$coi_filtering$max_length,
                sum(keep_l), length(lens)))
seqtab_nc <- seqtab_nc[, keep_l, drop = FALSE]

# ---- 8) NUMT screen por traducción -----------------------------------------
# Con justConcatenate=TRUE el ASV es R1 + 10×N + revComp(R2), discontinuo en
# coordenadas biológicas. La traducción a través del linker NNNNNNNNNN rompe
# el marco de lectura artificialmente, generando stop codons falsos. Por eso
# el check se DESACTIVA cuando hay concatenación. Si se quisiera, habría que
# traducir cada mitad por separado y luego unirlas — fuera de scope del workshop.
if (isTRUE(params$coi_filtering$apply_translation_check) && !just_concat) {
  message("[03] Check de traducción (código genético ",
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
  message("[03] ASVs con ORF abierto (no NUMT-like): ",
          sum(has_orf), "/", length(has_orf))
  seqtab_nc <- seqtab_nc[, has_orf, drop = FALSE]
} else if (just_concat) {
  message("[03] NUMT screen DESACTIVADO (incompatible con justConcatenate=TRUE).")
}

saveRDS(seqtab_nc, file.path(out_dir, "seqtab_nochim.rds"))

# ---- 9) Exportar ASVs.fasta con IDs estables -------------------------------
asv_seqs <- colnames(seqtab_nc)
asv_ids  <- sprintf("ASV_%04d", seq_along(asv_seqs))
names(asv_seqs) <- asv_ids
writeXStringSet(DNAStringSet(asv_seqs), file.path(out_dir, "asvs.fasta"))

asv_tab <- seqtab_nc
colnames(asv_tab) <- asv_ids
write_tsv(as_tibble(asv_tab, rownames = "sample_id"),
          file.path(out_dir, "asv_table.tsv"))

# ---- 10) Tabla de seguimiento de reads -------------------------------------
getN <- function(x) sum(getUniques(x))
basenames_ft <- sub("_R1\\.fastq\\.gz$", "", rownames(out_ft))

track <- tibble(
  sample_id    = sids,
  input        = out_ft[match(sids, basenames_ft), "reads.in"],
  q20_filtered = out_ft[match(sids, basenames_ft), "reads.out"],
  denoised_F   = vapply(dadaFs[sids], getN, integer(1)),
  denoised_R   = vapply(dadaRs[sids], getN, integer(1)),
  merged       = vapply(mergers[sids], getN, integer(1)),
  nonchim_filt = rowSums(seqtab_nc[sids, , drop = FALSE])
) |>
  mutate(
    pct_q20      = round(100 * q20_filtered / input,    1),
    pct_retained = round(100 * nonchim_filt / input,    1)
  )

write_tsv(track, file.path(out_dir, "track_reads.tsv"))

cat("\n=== Resumen 03_dada2_pipeline ===\n")
cat("Muestras procesadas:           ", nrow(track), "\n")
cat("ASVs finales (post-NUMT):      ", ncol(seqtab_nc), "\n")
cat("Reads totales conservados:     ", sum(seqtab_nc), "\n")
cat("Mediana % retenido tras Q20:   ", median(track$pct_q20),       "%\n")
cat("Mediana % retenido end-to-end: ", median(track$pct_retained),  "%\n")
cat("Outputs en:                    ", out_dir, "\n")
cat("\nAcciones recomendadas:\n")
cat(" - Revisar quality_postfilter_R{1,2}.pdf: confirmar Q20+ sostenido.\n")
cat(" - Si pct_q20 < 40 %% en >5 muestras: relajar truncQ a 15 en params.yml.\n")
cat(" - Construir DBs RESCRIPt: workflow/04a_build_BOLD_db_RESCRIPt.sh\n")
cat("                            workflow/04b_build_NCBI_db_RESCRIPt.sh\n")
cat(" - Asignar taxonomía:        workflow/04_taxonomy_BOLD_NCBI.R\n")
