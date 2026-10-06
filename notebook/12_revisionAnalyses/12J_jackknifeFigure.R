# Figure for analysis #5 (leave-one-subtribe-out jackknife, envPC1): dropped subtribes
# ordered by phylogeny, with the jackknife's rank-correlation stability (rho) and the
# envPC1 spread of the taxa each drop removes.
#
#   left   -- subtribe-level cladogram: the species tree pruned to one tip per subtribe,
#             labels coloured by tribe
#   middle -- rho of -log10(p) between the jackknife and full-set runs
#   right  -- envPC1 values of the dropped subtribe's taxa (box + points), labelled with
#             the whisker length (1.5 x IQR rule, as drawn) and the full-set percentiles
#             of the two whisker ends; dashed line = full-set median
#
# Subtribes are placed on the cladogram by their largest pure clade, because several are
# not strictly monophyletic (Anthistiriinae, Ischaeminae). No caption is drawn; it is
# written in the manuscript.
#
# Inputs are the 12J outputs, the species tree, taxonomy and the envPC trait table.
# envPC1 values are on the 707-taxon trait basis, as in 12J.
#
# Usage: Rscript notebook/12_revisionAnalyses/12J_jackknifeFigure.R

suppressMessages({
  library(ape); library(ggtree); library(ggplot2); library(patchwork); library(data.table)
})

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")
p <- function(...) file.path(PHYLOGWAS_ROOT, ...)

jk  <- fread(p("output/revisionAnalyses/jackknife_transitionLoss_fixedASR_envPC1.csv"))
env <- read.table(p("output/envData_707Poaceae_20250804.txt"), header = TRUE)
env <- data.table(assemblyID = rownames(env), envPC1 = env$envPC_1)
tax <- fread(p("data/Poaceae_taxonomic_groups_subtribeAdded_SKH20250402.csv"))[assemblyID %in% env$assemblyID]

tre <- read.tree(p("output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk"))
tre <- keep.tip(tre, intersect(tre$tip.label, env$assemblyID))
tre <- root(tre, "ASM1935983v1", resolve.root = TRUE)

# ---- one representative tip per subtribe ----
# Several subtribes are not strictly monophyletic (Anthistiriinae is paraphyletic
# around Andropogoneae's core; Ischaeminae is scattered). Each subtribe is placed by
# a tip from its largest pure clade, i.e. the largest clade whose tips all belong to
# that subtribe; this marks where the subtribe's core sits.
subs <- jk$subtribe
lab_of <- setNames(tax$subtribe, tax$assemblyID)
nt <- Ntip(tre)
clades <- phangorn::Descendants(tre, (nt + 1):(nt + tre$Nnode), type = "tips")
core <- sapply(subs, function(s) {
  pure <- Filter(function(d) all(lab_of[tre$tip.label[d]] %in% s), clades)
  best <- if (length(pure)) pure[[which.max(lengths(pure))]] else which(lab_of[tre$tip.label] %in% s)[1]
  n_s <- sum(lab_of[tre$tip.label] %in% s)
  message(sprintf("%-15s largest pure clade %d of %d tips", s, length(best), n_s))
  tre$tip.label[best[1]]
})
rep_tip <- core
st <- keep.tip(tre, unname(rep_tip))
st$tip.label <- names(rep_tip)[match(st$tip.label, rep_tip)]
st <- ladderize(st, right = FALSE)
st$edge.length <- NULL                      # cladogram: topology only

tribe_of <- tax[subtribe %in% subs, .(tribe = tribe[1]), by = subtribe]
jk <- merge(jk, tribe_of, by = "subtribe")

pt <- ggtree(st, ladderize = FALSE, branch.length = "none")
ord <- pt$data[pt$data$isTip, c("label", "y")]
jk[, y := ord$y[match(subtribe, ord$label)]]
n <- nrow(jk)
ylims <- c(.5, n + .5)

tribe_cols <- c(Andropogoneae = "#a0522d", Poeae = "#2e7d6b", Triticeae = "#4a6fa5", Oryzeae = "#7a5c99")
lab <- jk[, .(label = subtribe, tribe, tiplab = sprintf("%s (%d)", subtribe, taxa_dropped))]
pt <- pt %<+% lab + geom_tiplab(aes(label = tiplab, colour = tribe), size = 3.6, offset = .15,
                                fontface = "bold", show.legend = FALSE) +
  geom_point(data = function(d) d[d$isTip, ], aes(colour = tribe), size = 2.2) +
  scale_colour_manual(values = tribe_cols, name = "Tribe") +
  xlim(0, max(pt$data$x) + 6.5) + scale_y_continuous(limits = ylims, expand = c(0, 0)) +
  labs(title = "Dropped subtribe (taxa)") +
  theme(legend.position = "none", plot.title = element_text(size = 11, face = "bold"))

# ---- rho ----
prho <- ggplot(jk, aes(rho, y)) +
  geom_segment(aes(x = 1, xend = rho, yend = y), colour = "grey70") +
  geom_point(size = 3, colour = "#2b5d8a") +
  geom_text(aes(label = sprintf("%.2f", rho)), nudge_y = .32, size = 3) +
  scale_y_continuous(limits = ylims, expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(limits = c(.75, 1), breaks = seq(.75, 1, .05)) +
  labs(x = expression(rho ~ "(-log"[10] * "p, jackknife vs full)"), y = NULL,
       title = "Ranking stability") +
  theme_bw(11) + theme(plot.title = element_text(size = 11, face = "bold"),
                       panel.grid.minor = element_blank())

# ---- envPC1 spread ----
pts <- merge(tax[subtribe %in% subs, .(assemblyID, subtribe)], env, by = "assemblyID")
pts[, y := jk$y[match(subtribe, jk$subtribe)]]
# Where each subtribe's min and max fall in the full-set envPC1 distribution.
# Whisker extent as ggplot draws it: the most extreme points within 1.5 x IQR of the
# hinges (type-7 quantiles). Labelled with its length and the full-set percentiles of
# the two whisker ends.
whiskers <- function(x) {
  q <- quantile(x, c(.25, .75)); w <- 1.5 * diff(q)
  c(lo = min(x[x >= q[1] - w]), hi = max(x[x <= q[2] + w]))
}
pc_ecdf <- ecdf(env$envPC1)
spread <- pts[, {wk <- whiskers(envPC1)
                 .(wlo = wk[["lo"]], whi = wk[["hi"]],
                   plo = 100 * pc_ecdf(wk[["lo"]]), phi = 100 * pc_ecdf(wk[["hi"]]))}, by = y]
jk <- merge(jk, spread, by = "y")
jk[, spread_lab := sprintf("whisker %.1f (%.0f%%\u2013%.0f%%)", whi - wlo, plo, phi)]
print(jk[order(y, decreasing = TRUE), .(subtribe, rho, pc1_range, whisker = round(whi - wlo, 2),
                                         plo = round(plo), phi = round(phi))])
penv <- ggplot(pts, aes(envPC1, y)) +
  geom_vline(xintercept = median(env$envPC1), linetype = 2, colour = "grey50") +
  geom_boxplot(aes(group = y), width = .55, outlier.shape = NA, fill = "#e8b04a", alpha = .45) +
  geom_jitter(height = .12, width = 0, size = .8, alpha = .55, colour = "#8a5a10") +
  geom_text(data = jk, aes(x = min(env$envPC1), y = y, label = spread_lab),
            hjust = 0, nudge_y = .38, size = 2.9, colour = "grey30") +
  scale_y_continuous(limits = ylims, expand = c(0, 0), breaks = NULL) +
  coord_cartesian(xlim = range(env$envPC1)) +
  labs(x = "envPC1 (temperate \u2194 tropical)", y = NULL,
       title = "envPC1 spread of dropped taxa") +
  theme_bw(11) + theme(plot.title = element_text(size = 11, face = "bold"),
                       panel.grid.minor = element_blank())

fig <- pt + prho + penv + plot_layout(widths = c(1.2, 1, 1.5))
out <- p("output/revisionAnalyses/jackknife_phylo_rho_envPC1.png")
ggsave(out, fig, width = 13, height = 5.8, dpi = 300)
cat("wrote", out, "\n")
