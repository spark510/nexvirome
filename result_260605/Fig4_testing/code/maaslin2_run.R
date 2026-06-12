#!/usr/bin/env Rscript
# Generic MaAsLin2 Control-vs-Asthma differential-abundance run for one TPM
# relative-abundance matrix (sample x feature) + metadata (sample x group).
#
# MaAsLin2 default for relative-abundance/TPM input: normalization=TSS,
# transform=LOG, analysis_method=LM. Reference = Control (so a positive coef =
# higher in Asthma). BH-FDR q-values; in shotgun work q<0.25 is the conventional
# MaAsLin discovery threshold (we also report q<0.05/0.10).
#
# Usage:
#   Rscript maaslin2_run.R <tpm.csv> <meta.csv> <outdir> <label>
#     tpm.csv   sample x feature relative abundance (row.names = sample)
#     meta.csv  sample x metadata (must contain a 'group' column: Control/Asthma)
#     outdir    MaAsLin2 output directory
#     label     printed tag
suppressMessages(library(Maaslin2))

args <- commandArgs(trailingOnly = TRUE)
tpm_fn <- args[1]; meta_fn <- args[2]; outdir <- args[3]
label  <- if (length(args) >= 4) args[4] else basename(tpm_fn)

feat <- read.csv(tpm_fn,  row.names = 1, check.names = FALSE)   # sample x feature
meta <- read.csv(meta_fn, row.names = 1, check.names = FALSE)   # sample x group

# align
common <- intersect(rownames(feat), rownames(meta))
feat <- feat[common, , drop = FALSE]
meta <- meta[common, , drop = FALSE]
meta$group <- factor(meta$group, levels = c("Control", "Asthma"))  # ref = Control

nc <- sum(meta$group == "Control"); na <- sum(meta$group == "Asthma")
cat(sprintf("[%s] %d features x %d samples (Control %d / Asthma %d)\n",
            label, ncol(feat), nrow(feat), nc, na))

if (ncol(feat) < 1 || nc < 2 || na < 2) {
  cat("  too few features/samples; skipping\n"); quit(save = "no")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
fit <- Maaslin2(
  input_data       = feat,
  input_metadata   = meta,
  output           = outdir,
  fixed_effects    = c("group"),
  reference        = c("group,Control"),
  normalization    = "TSS",
  transform        = "LOG",
  analysis_method  = "LM",
  min_prevalence   = 0.1,
  max_significance = 0.25,        # MaAsLin's own q threshold for "significant"
  plot_heatmap     = FALSE,
  plot_scatter     = FALSE,
  cores            = 1
)

# summarise all_results
res <- fit$results
res <- res[order(res$qval), ]
cat(sprintf("  results: %d feature tests | q<0.05: %d | q<0.10: %d | q<0.25: %d\n",
            nrow(res),
            sum(res$qval < 0.05, na.rm = TRUE),
            sum(res$qval < 0.10, na.rm = TRUE),
            sum(res$qval < 0.25, na.rm = TRUE)))
hits <- res[!is.na(res$qval) & res$qval < 0.25, ]
if (nrow(hits)) {
  cat("  q<0.25 hits (coef>0 = higher in Asthma):\n")
  print(format(hits[, c("feature","value","coef","stderr","pval","qval","N.not.0")],
               digits = 3), row.names = FALSE)
} else cat("  no q<0.25 hits\n")
