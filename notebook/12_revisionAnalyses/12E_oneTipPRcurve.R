# Analysis #4: precision-recall of the one-tip-per-species association models against
# the full set. Positives = OGs with raw p < 0.001 in the full-set model; score =
# -log10 p in the one-tip-per-species model; OGs fitted in both runs only. The null
# PR-AUC equals the positive rate (a random ranking), drawn as a dashed line per envPC.
#
# Usage: Rscript notebook/12_revisionAnalyses/12E_oneTipPRcurve.R

suppressMessages({library(data.table); library(ggplot2); library(PRROC)})

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")
outDir <- file.path(PHYLOGWAS_ROOT, "output/revisionAnalyses")

curves <- list(); stats <- list()
for (k in 1:3) {
  full <- fread(file.path(PHYLOGWAS_ROOT, sprintf("output/finalModels_20251002/envPC_%d/ASREML_res.txt", k)),
                select = c("OG", "p"))
  sub  <- fread(file.path(outDir, sprintf("perSubsetModels/oneTipPerSpecies/envPC_%d/ASREML_res.txt", k)),
                select = c("OG", "p"))
  d <- merge(full, sub, by = "OG", suffixes = c("_full", "_sub"))[!is.na(p_full) & !is.na(p_sub)]
  pos <- d$p_full < 1e-3; score <- -log10(d$p_sub)
  pr <- pr.curve(scores.class0 = score[pos], scores.class1 = score[!pos], curve = TRUE)
  curves[[k]] <- data.table(envPC = paste0("envPC", k), recall = pr$curve[, 1], precision = pr$curve[, 2])
  stats[[k]] <- data.table(envPC = paste0("envPC", k), n = nrow(d), n_pos = sum(pos),
                           null_auc = mean(pos), pr_auc = pr$auc.davis.goadrich)
}
cv <- rbindlist(curves); st <- rbindlist(stats)
print(st)

pal <- c(envPC1 = "#1b6ca8", envPC2 = "#d0892b", envPC3 = "#3a9e6f")
st[, lab := sprintf("%s: AUC %.2f (null %.3f)", envPC, pr_auc, null_auc)]
p <- ggplot(cv, aes(recall, precision, colour = envPC)) +
  geom_hline(data = st, aes(yintercept = null_auc, colour = envPC), linetype = 2, linewidth = .4) +
  geom_line(linewidth = .7) +
  scale_colour_manual(values = pal, labels = setNames(st$lab, st$envPC), name = NULL) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Recall", y = "Precision",
       title = "One tip per species vs full set",
       subtitle = "Positives: full-set p < 0.001; score: one-tip -log10 p\nDashed lines: null AUC (positive rate)") +
  theme_bw(9) + theme(panel.grid.minor = element_blank(),
                      legend.position = c(.02, .1), legend.justification = c(0, 0),
                      legend.background = element_blank())
out <- file.path(outDir, "prcurve_oneTipPerSpecies_vs_full.png")
ggsave(out, p, width = 11, height = 10, units = "cm", dpi = 600)
cat("wrote", out, "\n")
