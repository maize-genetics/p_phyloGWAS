## Zero-shot score benchmark: does the LLM zero-shot score separate
## deleterious mutants (label 2, "_mut_") from neutral mutants (label 1,
## "_neutralMut_"), relative to their WT (label 0)?
##
## The score is unsupervised (no fitting happens anywhere in this pipeline),
## so there is nothing for a held-out gene split to protect against leakage on
## - all genes are used together for the reported curves.
##
## Usage:
##   Rscript src/visualization.R <fasta> <zeroShotScoreTSV> <outputDir> [outPrefix] [groupRegex]
##
## <fasta>            same FASTA used for steps 4 and 10
## <zeroShotScoreTSV> output of 10_logit2zeroShot.R (seqID, zeroShotScore, scaledZeroShotScore)
## <outputDir>        where to write figures + summary tables
## <outPrefix>        optional; prefix for output files (default "At_ESM2_650M")
## <groupRegex>       optional; regex passed to strsplit2 to recover the gene ID
##                     from a seqID (default "_| ", matching "<gene>_mut_N" /
##                     "<gene>_neutralMut_N" / "<gene> <rest of header>")

library(Biostrings)
library(PRROC)

args = commandArgs(trailingOnly = TRUE)
fastaPath   = args[1]
scoreTSV    = args[2]
outDir      = args[3]
outPrefix   = ifelse(length(args) >= 4, args[4], "At_ESM2_650M")
groupRegex  = ifelse(length(args) >= 5, args[5], "_| ")

dir.create(file.path(outDir, "figure"), recursive = TRUE, showWarnings = FALSE)

# ---- read inputs ----

Atfasta = readAAStringSet(fastaPath)
fastaID = sub("\\s.*$", "", names(Atfasta)) # first header token, matches npz/score keys

scoreTab = read.table(scoreTSV, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
zeroshotAt = setNames(scoreTab$scaledZeroShotScore, scoreTab$seqID)[fastaID] # align to FASTA order by name
names(zeroshotAt) = names(Atfasta)
if (any(is.na(zeroshotAt))) {
  warning(sum(is.na(zeroshotAt)), " FASTA records had no matching zero-shot score")
}

lab = rep(0, length(Atfasta))
lab[grepl("_mut_", names(Atfasta))] = 2
lab[grepl("_neutralMut_", names(Atfasta))] = 1

gene = limma::strsplit2(names(Atfasta), groupRegex)[,1]
allGeneSorted = sort(unique(gene))

# ---- per-gene mean score by label (WT / neutral / deleterious) ----

avgScore = t(sapply(split(seq_along(gene), gene)[allGeneSorted], function(idx) {
  tapply(zeroshotAt[idx], lab[idx], mean)[as.character(0:2)]
}))
colnames(avgScore) = c("WT", "neutralMut", "mut")
avgScoreTab = data.frame(Gene = allGeneSorted, avgScore, row.names = NULL)
write.table(avgScoreTab, file.path(outDir, paste0(outPrefix, "_perGeneMeanByLabel.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- gene-level z-scaling, then WT-relative contrast per mutant ----
# (unlist(tapply(...)) reorders into sorted-gene order; lab/seq are re-derived
# through the same tapply so everything stays aligned to that order)

zeroshotAt_scaleByGene = unlist(tapply(zeroshotAt, gene, scale))
lab_reorder  = unlist(tapply(lab, gene, function(x) x))
seqReorder   = unlist(tapply(as.character(Atfasta), gene, function(x) x))
geneOfScaled = unlist(tapply(gene, gene, function(x) x)) # gene of each entry, in the same tapply order

WT_zs  = setNames(zeroshotAt_scaleByGene[lab_reorder == 0], allGeneSorted)
WT_seq = setNames(seqReorder[lab_reorder == 0], allGeneSorted)

diff_zs     = zeroshotAt_scaleByGene - WT_zs[geneOfScaled]
WT_seq_all  = WT_seq[geneOfScaled]

outTab = data.frame(Gene = geneOfScaled,
                     WT_seq = WT_seq_all,
                     mut_seq = seqReorder,
                     Diff_zeroShot = diff_zs,
                     label = as.factor(lab_reorder))
outTab_filtered = outTab[outTab$label != 0,]

# ---- PR / ROC curves on all genes (unsupervised score - no split needed) ----

fg = na.omit(-outTab_filtered$Diff_zeroShot[outTab_filtered$label == 2]) # deleterious
bg = na.omit(-outTab_filtered$Diff_zeroShot[outTab_filtered$label == 1]) # neutral
pr  = pr.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE, rand.compute = TRUE)
roc = roc.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE, rand.compute = TRUE)

summaryTab = data.frame(set = "all genes",
                         n_deleterious = length(fg), n_neutral = length(bg),
                         PR_AUC = pr$auc.integral, ROC_AUC = roc$auc,
                         median_deleterious = median(fg), median_neutral = median(bg))
write.table(summaryTab, file.path(outDir, paste0(outPrefix, "_summary.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summaryTab)

# ---- Fig 3b: PR + ROC curves, all genes ----

png(file.path(outDir, "figure", paste0(outPrefix, "_PR_ROC.png")), height = 8.7, width = 8.7*4/3, units = "cm", pointsize = 8, res = 600)
par(mfrow = c(1,2))
plot(pr$curve, type = "l", ylim = c(0,1), col = "blue", lwd = 2,
     xlab = "Recall", ylab = "Precision",
     main = sprintf("PR curve (auc = %.3f)", pr$auc.integral))
points(pr$rand$curve, type = 'l', lty = 2, col = "grey")

plot(roc$curve, type = "l", ylim = c(0,1), col = "blue", lwd = 2,
     xlab = "FPR", ylab = "TPR",
     main = sprintf("ROC curve (auc = %.3f)", roc$auc))
points(roc$rand$curve, type = 'l', lty = 2, col = "grey")
dev.off()

# ---- supplementary: score distributions by label ----

png(file.path(outDir, "figure", paste0(outPrefix, "_scoreDist.png")), height = 8.7, width = 8.7*4/3, units = "cm", pointsize = 8, res = 600)
par(mfrow = c(1,2), cex.main = 0.85)
contrastByLabel = split(-outTab_filtered$Diff_zeroShot, outTab_filtered$label)
boxplot(contrastByLabel[c("1","2")], names = c("neutral", "deleterious"),
        ylab = "-Diff zero-shot score (WT-relative)", col = c("grey80","tomato"),
        main = "Score by mutation class")

plot(density(bg), col = "grey40", lwd = 2, xlab = "-Diff zero-shot score",
     main = "Score distributions (all genes)",
     xlim = range(c(fg, bg)))
lines(density(fg), col = "tomato", lwd = 2)
legend("topleft", legend = c("neutral", "deleterious"), col = c("grey40","tomato"), lwd = 2, bty = "n")
dev.off()

cat("\nWrote:\n",
    file.path(outDir, paste0(outPrefix, "_perGeneMeanByLabel.txt")), "\n",
    file.path(outDir, paste0(outPrefix, "_summary.txt")), "\n",
    file.path(outDir, "figure", paste0(outPrefix, "_PR_ROC.png")), "\n",
    file.path(outDir, "figure", paste0(outPrefix, "_scoreDist.png")), "\n")
