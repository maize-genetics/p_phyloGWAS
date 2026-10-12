## Overlay PR / ROC curves of the zero-shot deleterious-vs-neutral benchmark
## for several datasets on one figure.
##
## Usage:
##   Rscript src/compare_datasets.R <outPNG> <label1>:<fasta1>:<scoreTSV1> [<label2>:<fasta2>:<scoreTSV2> ...]
##
## Each dataset is processed exactly as in visualization.R: gene-level z-scaling
## of scaledZeroShotScore, WT-relative contrast (-Diff_zeroShot), deleterious
## (_mut_) as foreground vs neutral (_neutralMut_) as background. Gene key is
## strsplit2(header, "_| ")[,1].

library(Biostrings)
library(PRROC)

args = commandArgs(trailingOnly = TRUE)
outPNG = args[1]
specs  = args[-1]

contrast = function(fastaPath, scoreTSV, groupRegex = "_| ") {
  fa = readAAStringSet(fastaPath)
  fastaID = sub("\\s.*$", "", names(fa))
  st = read.table(scoreTSV, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  zs = setNames(st$scaledZeroShotScore, st$seqID)[fastaID]
  lab = rep(0, length(fa)); lab[grepl("_mut_", names(fa))] = 2; lab[grepl("_neutralMut_", names(fa))] = 1
  gene = limma::strsplit2(names(fa), groupRegex)[,1]
  zsg  = unlist(tapply(zs, gene, scale))
  labr = unlist(tapply(lab, gene, function(x) x))
  gr   = unlist(tapply(gene, gene, function(x) x))
  wt   = setNames(zsg[labr == 0], sort(unique(gene)))
  d    = -(zsg - wt[gr])
  list(fg = na.omit(d[labr == 2]), bg = na.omit(d[labr == 1]))
}

cols = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
res = list()
for (i in seq_along(specs)) {
  p = strsplit(specs[i], ":")[[1]]
  ct = contrast(p[2], p[3])
  res[[p[1]]] = list(
    pr  = pr.curve(scores.class0 = ct$fg, scores.class1 = ct$bg, curve = TRUE, rand.compute = TRUE),
    roc = roc.curve(scores.class0 = ct$fg, scores.class1 = ct$bg, curve = TRUE, rand.compute = TRUE),
    n_fg = length(ct$fg), n_bg = length(ct$bg))
}

png(outPNG, height = 8.7, width = 8.7*2, units = "cm", pointsize = 8, res = 600)
par(mfrow = c(1,2))
plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "Recall", ylab = "Precision", main = "PR curve")
for (i in seq_along(res)) {
  lines(res[[i]]$pr$curve[,1:2], col = cols[i], lwd = 2)
  abline(h = res[[i]]$n_fg / (res[[i]]$n_fg + res[[i]]$n_bg), col = cols[i], lty = 3)
}
legend("bottomleft", bty = "n", col = cols[seq_along(res)], lwd = 2,
       legend = sprintf("%s (auc = %.3f)", names(res), sapply(res, function(r) r$pr$auc.integral)))

plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "FPR", ylab = "TPR", main = "ROC curve")
abline(0, 1, col = "grey", lty = 2)
for (i in seq_along(res)) lines(res[[i]]$roc$curve[,1:2], col = cols[i], lwd = 2)
legend("bottomright", bty = "n", col = cols[seq_along(res)], lwd = 2,
       legend = sprintf("%s (auc = %.3f)", names(res), sapply(res, function(r) r$roc$auc)))
dev.off()

for (n in names(res)) cat(sprintf("%-8s n_del=%d n_neut=%d PR=%.3f ROC=%.3f\n", n, res[[n]]$n_fg, res[[n]]$n_bg,
                                   res[[n]]$pr$auc.integral, res[[n]]$roc$auc))
cat("wrote", outPNG, "\n")
