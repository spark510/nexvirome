#!/usr/bin/env Rscript
# Generic MaAsLin3 Control-vs-Asthma differential-abundance run for one TPM
# relative-abundance matrix (sample x feature) + metadata (sample x group).
#
# MaAsLin3 fits BOTH an abundance model (CLR-style, among detected) AND a
# prevalence model (presence/absence), each FDR-corrected (BH qval). Reference =
# Control, so a positive coefficient = higher / more prevalent in Asthma.
# Input is relative abundance (TPM); MaAsLin3 default normalization=TSS,
# transform=LOG.
#
# Usage:
#   Rscript maaslin3_run.R <tpm.csv> <meta.csv> <outdir> <label>
#     tpm.csv   sample x feature relative abundance (row.names = sample)
#     meta.csv  sample x metadata (must contain 'group': Control/Asthma)
#     outdir    MaAsLin3 output directory
#     label     printed tag
suppressMessages(library(maaslin3))

args <- commandArgs(trailingOnly = TRUE)
tpm_fn <- args[1]; meta_fn <- args[2]; outdir <- args[3]
label  <- if (length(args) >= 4) args[4] else basename(tpm_fn)

feat <- read.csv(tpm_fn,  row.names = 1, check.names = FALSE)   # sample x feature
meta <- read.csv(meta_fn, row.names = 1, check.names = FALSE)   # sample x group

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
fit <- tryCatch(
  maaslin3(
    input_data       = feat,
    input_metadata   = meta,
    output           = outdir,
    formula          = ~ group,
    reference        = c("group,Control"),
    normalization    = "TSS",
    transform        = "LOG",
    min_prevalence   = 0.10,
    max_significance = 0.25,            # MaAsLin discovery q threshold
    plot_summary_plot = FALSE,
    plot_associations = FALSE,
    cores            = 1),
  error = function(e) { cat("  MaAsLin3 error:", conditionMessage(e), "\n"); NULL })

if (is.null(fit)) quit(save = "no")

# all_results.tsv holds both abundance and prevalence associations
res_fn <- file.path(outdir, "all_results.tsv")
if (!file.exists(res_fn)) { cat("  no all_results.tsv produced\n"); quit(save="no") }
res <- read.delim(res_fn, stringsAsFactors = FALSE)
res <- res[order(res$qval_individual), ]

cat(sprintf("  %d associations | q<0.05: %d | q<0.10: %d | q<0.25: %d\n",
            nrow(res),
            sum(res$qval_individual < 0.05, na.rm = TRUE),
            sum(res$qval_individual < 0.10, na.rm = TRUE),
            sum(res$qval_individual < 0.25, na.rm = TRUE)))
hits <- res[!is.na(res$qval_individual) & res$qval_individual < 0.25, ]
if (nrow(hits)) {
  cat("  q<0.25 hits (coef>0 = higher/more prevalent in Asthma):\n")
  show <- intersect(c("feature","metadata","value","association","coef",
                      "pval_individual","qval_individual","N","N_not_zero"),
                    colnames(hits))
  print(format(hits[, show], digits = 3), row.names = FALSE)
} else cat("  no q<0.25 hits\n")
