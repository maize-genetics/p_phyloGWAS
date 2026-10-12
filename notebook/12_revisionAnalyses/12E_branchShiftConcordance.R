# Analysis #4 follow-up: per-branch state change, full set (707) vs one-tip-per-species
# (569), on every branch present in both trees.
#
# The transition call (06B, 12E) needs (a) non-overlapping 95% CIs of the parent and
# child ancestral states and (b) crossing the 30th/70th percentile cutoff. Part (a) is
# turned into a test here:
#   z = (child - parent) / sqrt(var_parent + var_child)
# from fastAnc's state estimates and variances (tips have variance 0). The covariance
# between parent and child is ignored, so the test is conservative. Part (b) has no
# p-value; it is shown only as which branches are called transitions, at clade level
# (nested calls merged with keep_most_ancestral, as in the transition counts).
#
# Retained branches: a one-tip branch is defined by its child clade C and parent clade P
# (tip sets on the 569 taxa). It is retained if C and P are also clades of the full tree
# once restricted to the 569 taxa; its full-tree counterpart runs from MRCA(P) to
# MRCA(C). Where dropped duplicate tips hang off that path, the full-tree "branch" is a
# short chain of edges, compared as one change.
#
# Usage: Rscript notebook/12_revisionAnalyses/12E_branchShiftConcordance.R

suppressMessages({
  library(ape); library(phangorn); library(phytools); library(data.table); library(ggplot2)
})

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")
outDir <- file.path(PHYLOGWAS_ROOT, "output/revisionAnalyses/oneTipPerSpeciesASR")

# --- same inputs and rooting as 12E / 06B ---
metadata <- read.delim(file.path(PHYLOGWAS_ROOT, "data/Poaceae_metadata_filtered_2025.08.28.tsv"))
spTre <- read.tree(file.path(PHYLOGWAS_ROOT, "output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk"))
spTre.rooted <- root(spTre, "ASM1935983v1", resolve.root = TRUE)
spTre.rooted$edge.length[is.na(spTre.rooted$edge.length)] <- 0.1
ePCs <- readRDS(file.path(PHYLOGWAS_ROOT, "output/ePC_20250804.rds"))
traits <- ePCs$synthetic.environmental.traits

fullTips <- Reduce(intersect, list(spTre.rooted$tip.label, rownames(ePCs$environmental.features),
                                   metadata$assemblyID, rownames(traits)))
oneTips <- intersect(fullTips, rownames(read.table(
  file.path(PHYLOGWAS_ROOT, "output/revisionAnalyses/traitTables/envData_oneTipPerSpecies.txt"), header = TRUE)))
treeFull <- multi2di(keep.tip(spTre.rooted, fullTips))
treeOne  <- multi2di(keep.tip(spTre.rooted, oneTips))
cat(sprintf("full tree %d tips, one-tip tree %d tips\n", Ntip(treeFull), Ntip(treeOne)))

# keep_most_ancestral(): as in 12E / 06B cell 45 -- of nested calls, keep the oldest.
keep_most_ancestral <- function(tree, nodes) {
  nodes <- unique(as.integer(nodes))
  if (length(nodes) <= 1) return(nodes)
  depth <- node.depth.edgelength(tree)[nodes]
  tipset <- lapply(nodes, function(nd) sort(Descendants(tree, nd, type = "tips")[[1]]))
  redundant <- logical(length(nodes))
  for (i in seq_along(nodes)) for (j in seq_along(nodes)) {
    if (i == j) next
    A <- tipset[[i]]; B <- tipset[[j]]
    if (length(A) <= length(B) && all(A %in% B)) {
      if (length(A) != length(B) || depth[i] > depth[j]) { redundant[i] <- TRUE; break }
    }
  }
  nodes[!redundant]
}
# Transition clades of a tree under the 06B rule, after merging nested calls.
transition_clades <- function(tree, s, cut) {
  E <- tree$edge; a <- E[, 1]; b <- E[, 2]
  toL <- s$lo[a] > s$up[b] & s$ace[b] < cut[1] & s$ace[a] > cut[1]
  toH <- s$lo[b] > s$up[a] & s$ace[b] > cut[2] & s$ace[a] < cut[2]
  list(toH = keep_most_ancestral(tree, b[toH]), toL = keep_most_ancestral(tree, b[toL]))
}

node_states <- function(tree, trait) {
  fit <- fastAnc(tree, trait[tree$tip.label], var = TRUE, CI = TRUE)
  n <- Ntip(tree)
  data.table(node = 1:(n + tree$Nnode),
             ace = c(trait[tree$tip.label], fit$ace), var = c(rep(0, n), fit$var),
             lo = c(trait[tree$tip.label], fit$CI95[, 1]), up = c(trait[tree$tip.label], fit$CI95[, 2]))
}
key <- function(tips) paste(sort(tips), collapse = "|")
clade_tips <- function(tree, nd) tree$tip.label[Descendants(tree, nd, type = "tips")[[1]]]

# Full-tree node for each one-tip clade (restricted to one-tip taxa), if it is a clade there.
nFull <- Ntip(treeFull) + treeFull$Nnode
fullKey <- vapply(seq_len(nFull), function(nd) key(intersect(clade_tips(treeFull, nd), oneTips)), "")
# Several full nodes can share a restricted tip set (chains above dropped tips); MRCA = the
# deepest one, i.e. the one with the most full-tree tips.
fullSize <- vapply(seq_len(nFull), function(nd) length(clade_tips(treeFull, nd)), 1L)
ord <- order(fullSize)                       # later (larger) entries overwrite smaller ones...
mrcaOf <- setNames(seq_len(nFull)[ord], fullKey[ord])
mrcaOf <- mrcaOf[!duplicated(names(mrcaOf), fromLast = TRUE)]   # ...keep the deepest
# ...but the child end of a branch should be the shallowest node with that set.
childOf <- setNames(seq_len(nFull)[ord], fullKey[ord])
childOf <- childOf[!duplicated(names(childOf))]

rows <- list(); cladeCounts <- list()
for (i in 1:3) {
  trait <- setNames(traits[, i], rownames(traits))
  sF <- node_states(treeFull, trait); sO <- node_states(treeOne, trait)
  cutFull <- quantile(trait[fullTips], c(.3, .7))
  # Clade-level calls (nested calls merged), keyed by descendant taxa on the one-tip set.
  # A full-tree clade made only of dropped duplicate tips has an empty key ("gone").
  cF <- transition_clades(treeFull, sF, cutFull); cO <- transition_clades(treeOne, sO, cutFull)
  cladeKeys <- function(tree, cl, restrict) unlist(lapply(names(cl), function(dir)
    setNames(rep(dir, length(cl[[dir]])),
             vapply(cl[[dir]], function(nd) key(intersect(clade_tips(tree, nd), restrict)), ""))))
  kF <- cladeKeys(treeFull, cF, oneTips); kO <- cladeKeys(treeOne, cO, oneTips)
  gone <- names(kF) == ""
  bothDir <- intersect(paste(names(kF), kF), paste(names(kO), kO))
  cladeCounts[[i]] <- data.table(envPC = i,
    full_clades = length(kF), onetip_clades = length(kO),
    both = length(bothDir),
    full_only = sum(!gone) - length(bothDir), gone_with_dropped_tips = sum(gone),
    onetip_only = length(kO) - length(bothDir))
  # Harmonised full-set calls: the 06B rule applied to each one-tip branch's whole stretch
  # of the full tree (MRCA(P) -> MRCA(C)), so a shift split over 2-3 full-tree edges by a
  # dropped duplicate tip counts once; nested calls are then merged on the one-tip tree,
  # the same way as the one-tip calls.
  rule <- function(s, a, b) {
    if (s$lo[a] > s$up[b] && s$ace[b] < cutFull[1] && s$ace[a] > cutFull[1]) "toL"
    else if (s$lo[b] > s$up[a] && s$ace[b] > cutFull[2] && s$ace[a] < cutFull[2]) "toH" else "none"
  }
  pathCall <- setNames(rep("none", Ntip(treeOne) + treeOne$Nnode), seq_len(Ntip(treeOne) + treeOne$Nnode))
  E <- treeOne$edge
  for (e in seq_len(nrow(E))) {
    p <- E[e, 1]; ch <- E[e, 2]
    kP <- key(clade_tips(treeOne, p)); kC <- key(clade_tips(treeOne, ch))
    if (!(kP %in% names(mrcaOf)) || !(kC %in% names(childOf))) next    # not a clade in the full tree
    fp <- mrcaOf[[kP]]; fc <- childOf[[kC]]
    pathCall[ch] <- rule(sF, fp, fc)
    z <- function(s, a, b) (s$ace[b] - s$ace[a]) / sqrt(s$var[a] + s$var[b])
    rows[[length(rows) + 1]] <- data.table(envPC = i, child = kC, terminal = ch <= Ntip(treeOne),
      z_full = z(sF, fp, fc), z_onetip = z(sO, p, ch),
      call_full = if (kC %in% names(kF)) kF[[kC]] else "none",        # merged clade calls
      call_onetip = if (kC %in% names(kO)) kO[[kC]] else "none", node_onetip = ch)
  }
  hF <- list(toH = keep_most_ancestral(treeOne, which(pathCall == "toH")),
             toL = keep_most_ancestral(treeOne, which(pathCall == "toL")))
  kH <- cladeKeys(treeOne, hF, oneTips)
  bothH <- intersect(paste(names(kH), kH), paste(names(kO), kO))
  cladeCounts[[i]][, `:=`(harmonised_full_clades = length(kH), harmonised_both = length(bothH),
                          harmonised_full_only = length(kH) - length(bothH),
                          harmonised_onetip_only = length(kO) - length(bothH))]
  for (r in seq_along(rows)) if (rows[[r]]$envPC == i)
    rows[[r]][, call_full_harmonised := if (child %in% names(kH)) kH[[child]] else "none"]
}
br <- rbindlist(rows)
br[, `:=`(p_full = 2 * pnorm(-abs(z_full)), p_onetip = 2 * pnorm(-abs(z_onetip)))]
# Colour by the harmonised comparison (full-set calls on whole stretches, merged on the
# one-tip tree); the edge-by-edge full-tree calls stay in call_full.
br[, transition := fifelse(call_full_harmonised != "none" & call_full_harmonised == call_onetip, "both",
                   fifelse(call_full_harmonised != "none", "full only",
                   fifelse(call_onetip != "none", "one-tip only", "neither")))]
fwrite(br[, !c("child", "node_onetip")], file.path(outDir, "branchShift_full_vs_oneTip.csv"))

nOneEdges <- nrow(treeOne$edge)
summ <- br[, .(retained = .N, of_onetip_branches = nOneEdges,
               rho_z = round(cor(z_full, z_onetip, method = "spearman"), 3),
               rho_logp = round(cor(-log10(p_full), -log10(p_onetip), method = "spearman"), 3),
               rho_z_internal = round(cor(z_full[!terminal], z_onetip[!terminal], method = "spearman"), 3),
               trans_both = sum(transition == "both"), trans_full_only = sum(transition == "full only"),
               trans_onetip_only = sum(transition == "one-tip only")), by = envPC]
print(summ)
cc <- rbindlist(cladeCounts)
cat("\nTransition clades (nested calls merged; full-set cutoffs).\n",
    "harmonised_*: full-set rule on whole stretches, merged on the one-tip tree;\n",
    "gone_with_dropped_tips applies to both comparisons.\n"); print(cc)
fwrite(cc, file.path(outDir, "branchShift_transitionClades_full_vs_oneTip.csv"))
fwrite(summ, file.path(outDir, "branchShift_full_vs_oneTip_summary.csv"))

# --- scatter: signed z, full vs one-tip ---
lab <- summ[, .(envPC, txt = sprintf("ρ = %.2f (all)\nρ = %.2f (internal)\nn = %d", rho_z, rho_z_internal, retained))]
br[, envPC_lab := factor(paste0("envPC", envPC))]; lab[, envPC_lab := factor(paste0("envPC", envPC))]
br <- br[order(transition != "neither")]          # transitions drawn on top
pal <- c(neither = "grey80", both = "#2471a3", `full only` = "#c0392b", `one-tip only` = "#e69f00")
p <- ggplot(br, aes(z_full, z_onetip, colour = transition)) +
  geom_abline(linetype = 2, colour = "grey50", linewidth = .3) +
  geom_point(size = .6, alpha = .7) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE,
            hjust = -.08, vjust = 1.15, size = 2.6, lineheight = .9) +
  facet_wrap(~ envPC_lab, nrow = 1) +
  scale_colour_manual(values = pal, name = "Transition clade") +
  labs(x = "Branch shift z, full set (707)", y = "Branch shift z, one tip per species (569)") +
  theme_bw(8) + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                      strip.background = element_blank())
ggsave(file.path(outDir, "branchShift_full_vs_oneTip.png"), p, width = 17.4, height = 7.5, units = "cm", dpi = 600)
cat("wrote", file.path(outDir, "branchShift_full_vs_oneTip.{csv,png}"), "\n")
