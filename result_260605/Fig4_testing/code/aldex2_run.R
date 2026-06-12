#!/usr/bin/env Rscript
# Generic ALDEx2 Control-vs-Asthma differential-abundance run for one count
# matrix. ALDEx2 takes RAW INTEGER COUNTS (it does CLR + Dirichlet Monte-Carlo
# internally). Welch-t (we.*) and Wilcoxon (wi.*) tests + CLR effect size, with
# BH-FDR (*.eBH). Reference group = "Control" (so positive diff.btw / effect =
# higher in Asthma).
#
# Usage:
#   Rscript aldex2_run.R <counts.csv> <sample_group.csv> <out.csv> <label>
#     counts.csv        feature x sample integer matrix (row.names = feature)
#     sample_group.csv  columns: sample, group  (Control / Asthma)
#     out.csv           output table
#     label             printed tag (e.g. "RNA non-phage genus")
#
# Output columns: feature, prevalence, n_Control, n_Asthma, we.ep, we.eBH,
#   wi.ep, wi.eBH, effect, diff.btw, diff.win, rab.all, rab.Control, rab.Asthma
suppressMessages(library(ALDEx2))

args <- commandArgs(trailingOnly = TRUE)
cts_fn <- args[1]; grp_fn <- args[2]; out_fn <- args[3]
label  <- if (length(args) >= 4) args[4] else basename(cts_fn)

cts <- as.matrix(read.csv(cts_fn, row.names = 1, check.names = FALSE))
grp <- read.csv(grp_fn, stringsAsFactors = FALSE)

cond <- grp$group[match(colnames(cts), grp$sample)]
keep <- !is.na(cond)
cts  <- cts[, keep, drop = FALSE]; cond <- cond[keep]
# reference = Control: make it the first factor level
cond <- factor(cond, levels = c("Control", "Asthma"))

prev <- rowSums(cts > 0)
cts  <- cts[rowSums(cts) > 0, , drop = FALSE]
n_ctrl <- sum(cond == "Control"); n_asth <- sum(cond == "Asthma")
cat(sprintf("[%s] %d features x %d samples (Control %d / Asthma %d)\n",
            label, nrow(cts), ncol(cts), n_ctrl, n_asth))

if (nrow(cts) < 1 || n_ctrl < 2 || n_asth < 2) {
  cat("  too few features/samples for ALDEx2; writing empty result\n")
  write.csv(data.frame(), out_fn, row.names = FALSE); quit(save = "no")
}

set.seed(1)
x  <- aldex.clr(cts, conds = as.character(cond), mc.samples = 128,
                denom = "all", verbose = FALSE)
tt <- aldex.ttest(x)                 # we.ep/we.eBH/wi.ep/wi.eBH
ef <- aldex.effect(x)                # effect/diff.btw/diff.win/rab.*

res <- data.frame(feature = rownames(tt),
                  prevalence = prev[rownames(tt)],
                  n_Control = n_ctrl, n_Asthma = n_asth,
                  tt, ef, row.names = NULL, check.names = FALSE)
res <- res[order(res$wi.ep), ]
write.csv(res, out_fn, row.names = FALSE)

# console summary
cat(sprintf("  -> %s\n", out_fn))
cat(sprintf("  wi.ep<0.05: %d | wi.eBH<0.05: %d | wi.eBH<0.10: %d | |effect|>1: %d\n",
            sum(res$wi.ep  < 0.05, na.rm = TRUE),
            sum(res$wi.eBH < 0.05, na.rm = TRUE),
            sum(res$wi.eBH < 0.10, na.rm = TRUE),
            sum(abs(res$effect) > 1, na.rm = TRUE)))
sig <- res[!is.na(res$wi.eBH) & res$wi.eBH < 0.10, ]
if (nrow(sig)) {
  cat("  FDR<0.10 hits:\n")
  print(format(sig[, c("feature","prevalence","effect","diff.btw","wi.ep","wi.eBH")],
               digits = 3), row.names = FALSE)
} else cat("  no FDR<0.10 hits\n")
