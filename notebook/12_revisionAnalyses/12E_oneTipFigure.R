# Analysis #4 summary figure: one tip per species (569) vs full set (707), envPC1-3.
#   row 1 -- per-branch shift z (ASR), full vs one-tip, on all 1,136 shared branches;
#            coloured by transition-clade call (harmonised comparison, see
#            12E_branchShiftConcordance.R)
#   row 2 -- association-model significance, -log10 p full vs one-tip, per OG
#   row 3 -- precision-recall: positives = full-set p < 0.001, score = one-tip -log10 p;
#            dashed = null AUC (positive rate)
#
# Inputs: output/revisionAnalyses/oneTipPerSpeciesASR/branchShift_full_vs_oneTip.csv
# (from 12E_branchShiftConcordance.R) and the full / one-tip ASREML_res.txt tables.
# Output: output/revisionAnalyses/one_tip_per_species.{png,pdf}, 17.4 cm wide.
#
# Usage: Rscript notebook/12_revisionAnalyses/12E_oneTipFigure.R

suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(PRROC)})

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")
outDir <- file.path(PHYLOGWAS_ROOT, "output/revisionAnalyses")
FS <- 7                                   # base font size, pt
P_POS <- 1e-3

th <- theme_bw(FS) + theme(panel.grid.minor = element_blank(), plot.title = element_text(size = FS + 1, hjust = .5),
                           legend.position = "none", plot.background = element_blank())
# Panel text in a fixed corner, a few pt in from the panel edge, lines aligned to that side.
corner <- function(txt, where = "tl") {
  left <- grepl("l", where); top <- grepl("t", where)
  annotation_custom(grid::textGrob(txt,
    x = if (left) unit(4, "pt") else unit(1, "npc") - unit(4, "pt"),
    y = if (top) unit(1, "npc") - unit(4, "pt") else unit(.12, "npc"),   # bottom: clear of null line
    hjust = if (left) 0 else 1, vjust = if (top) 1 else 0,
    gp = grid::gpar(fontsize = FS - .5, lineheight = .95)))
}

# Scatters use the same range on both axes (1:1 aspect), shared across the panels of a row.
equal_axes <- function(x, y) { r <- range(c(x, y)); coord_fixed(xlim = r, ylim = r) }

# ---- row 1: ASR branch shift z ----
br <- fread(file.path(outDir, "oneTipPerSpeciesASR/branchShift_full_vs_oneTip.csv"))
pal_tr <- c(neither = "grey80", both = "#2471a3", `full only` = "#c0392b", `one-tip only` = "#e69f00")
row1 <- lapply(1:3, function(k) {
  d <- br[envPC == k][order(transition != "neither")]
  rho <- cor(d$z_full, d$z_onetip, method = "spearman")
  ggplot(d, aes(z_full, z_onetip, colour = transition)) +
    geom_abline(linetype = 2, colour = "grey50", linewidth = .3) +
    geom_point(size = .7, alpha = .8) +
    scale_colour_manual(values = pal_tr, name = "Transition clade",
                        breaks = c("both", "full only", "one-tip only", "neither")) +
    corner(sprintf("ρ = %.2f\nn = %s branches", rho, format(nrow(d), big.mark = ","))) +
    equal_axes(br$z_full, br$z_onetip) +
    labs(title = paste0("envPC", k), x = "Branch shift z, full set",
         y = if (k == 1) "Branch shift z, one tip" else NULL) + th +
    # Transition-clade legend inside the envPC3 panel only, in its empty lower-right corner.
    (if (k == 3) theme(legend.position = c(.98, .02), legend.justification = c(1, 0),
                       legend.title = element_text(size = FS - .5), legend.text = element_text(size = FS - 1),
                       legend.key.size = unit(FS, "pt"), legend.background = element_blank())
     else NULL)
})
row1[[3]] <- row1[[3]] + guides(colour = guide_legend(override.aes = list(size = 1.4, alpha = 1)))

# ---- rows 2-3: association model ----
sig <- lapply(1:3, function(k) {
  full <- fread(file.path(PHYLOGWAS_ROOT, sprintf("output/finalModels_20251002/envPC_%d/ASREML_res.txt", k)),
                select = c("OG", "p"))
  sub  <- fread(file.path(outDir, sprintf("perSubsetModels/oneTipPerSpecies/envPC_%d/ASREML_res.txt", k)),
                select = c("OG", "p"))
  d <- merge(full, sub, by = "OG", suffixes = c("_full", "_sub"))[!is.na(p_full) & !is.na(p_sub)]
  d[, `:=`(lf = -log10(p_full), ls = -log10(p_sub), pos = p_full < P_POS)]
  d
})
pal_sig <- c(`TRUE` = "#2471a3", `FALSE` = "grey75")
CUT_FULL <- "#c0392b"   # full-set p = 0.001 (defines positives)
CUT_SUB  <- "#e69f00"   # one-tip p = 0.05 (recall threshold)
row2 <- lapply(1:3, function(k) {
  d <- sig[[k]][order(pos)]
  rho <- cor(d$lf, d$ls, method = "spearman")
  rec <- mean(d$ls[d$pos] > -log10(.05))
  ggplot(d, aes(lf, ls, colour = pos)) +
    geom_abline(linetype = 2, colour = "grey50", linewidth = .3) +
    geom_vline(xintercept = -log10(P_POS), linetype = 2, colour = CUT_FULL, linewidth = .4) +
    geom_hline(yintercept = -log10(.05), linetype = 2, colour = CUT_SUB, linewidth = .4) +
    geom_point(size = .55, alpha = .6) +
    scale_colour_manual(values = pal_sig) +
    corner(sprintf("ρ = %.2f\nrecall (p < 0.05) = %.1f%%\nn = %s OGs", rho, 100 * rec,
                   format(nrow(d), big.mark = ",")), "tl") +
    equal_axes(unlist(lapply(sig, `[[`, "lf")), unlist(lapply(sig, `[[`, "ls"))) +
    labs(x = expression(-log[10] * p ~ "full set"),
         y = if (k == 1) expression(-log[10] * p ~ "one tip") else NULL) + th
})
row3 <- lapply(1:3, function(k) {
  d <- sig[[k]]
  pr <- pr.curve(scores.class0 = d$ls[d$pos], scores.class1 = d$ls[!d$pos], curve = TRUE)
  cv <- data.table(recall = pr$curve[, 1], precision = pr$curve[, 2])
  ggplot(cv, aes(recall, precision)) +
    geom_hline(yintercept = mean(d$pos), linetype = 2, colour = "grey45", linewidth = .35) +
    geom_line(colour = "#2471a3", linewidth = .6) +
    corner(sprintf("PR-AUC = %.2f\nnull AUC = %.3f\n%d positives", pr$auc.davis.goadrich,
                   mean(d$pos), sum(d$pos)), "bl") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Recall", y = if (k == 1) "Precision" else NULL) + th
})

# Line and colour key goes in the manuscript caption: rows 1-2 dashed grey = 1:1; rows 2-3
# blue = OGs with full-set p < 0.001; row 2 red dashed = full-set p = 0.001, orange
# dashed = one-tip p = 0.05.
fig <- wrap_plots(c(row1, row2, row3), ncol = 3)
for (ext in c("png", "pdf")) {
  f <- file.path(outDir, paste0("one_tip_per_species.", ext))
  if (ext == "png") ggsave(f, fig, width = 17.4, height = 18, units = "cm", dpi = 600)
  else ggsave(f, fig, width = 17.4, height = 18, units = "cm", device = cairo_pdf)
}
cat("wrote", file.path(outDir, "one_tip_per_species.{png,pdf}"), "\n")
