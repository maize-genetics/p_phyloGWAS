# Analysis #5 follow-up (revision-analyses.md): how many ASR transition events does
# each leave-one-subtribe-out jackknife actually remove?
#
# METHOD NOTE -- why the ASR is computed ONCE and then held fixed:
# The obvious implementation is to re-run the ASR on each jackknifed taxon set and
# compare transition counts. That conflates two different things: (a) transition-
# bearing lineages genuinely removed by the drop, and (b) re-estimation effects --
# dropping taxa shifts the quantile(trait,.3)/quantile(trait,.7) cutoffs and the
# fastAnc CIs, so edges that previously missed the threshold get newly called. Run
# that way, three subtribes (Oryzinae +3, Triticinae/Aveninae +2) appear to *gain*
# transitions despite removing zero transition-carrying tips -- a pure artifact.
#
# So: the full-707 ASR is computed once, its transition clades are held FIXED, and
# each jackknife is scored purely structurally -- what fraction of each fixed clade's
# tips survive the drop. No re-detection, hence no "lost vs retained" call to make;
# a clade is either wiped (no tips left), partially eroded, or untouched.
#
# Empirically the partial category is nearly empty: transition clades and subtribes
# are both monophyletic, so a clade is almost always wholly inside or wholly outside
# a dropped subtribe. "Wiped" is therefore the natural unit of transition loss.
#
# envPC1 only, matching the envPC1-only jackknife runs. Lightweight: one fastAnc on
# 707 tips plus set arithmetic -- seconds.
#
# Usage: Rscript notebook/12_revisionAnalyses/12J_jackknifeTransitionLoss.R

suppressMessages({
  library(ape); library(phangorn); library(phytools); library(ggtree)
})

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")
ENVPC_IDX      <- as.integer(Sys.getenv("ENVPC_IDX", unset = "1"))
MIN_SUBTRIBE_N <- as.integer(Sys.getenv("MIN_SUBTRIBE_N", unset = "21"))  # >20 taxa

# --- keep_most_ancestral(): copied verbatim from 06B_visualizationEnvAdapt.ipynb cell 45 ---
keep_most_ancestral <- function(tree, nodes) {
  stopifnot(inherits(tree, "phylo"))
  nodes <- unique(as.integer(nodes))
  max_id <- Ntip(tree) + tree$Nnode
  nodes <- nodes[nodes >= 1 & nodes <= max_id]
  if (length(nodes) <= 1) return(nodes)
  depth  <- node.depth.edgelength(tree)[nodes]
  tipset <- lapply(nodes, function(nd) sort(Descendants(tree, nd, type = "tips")[[1]]))
  is_subset <- function(A, B) length(A) <= length(B) && all(A %in% B)
  redundant <- logical(length(nodes))
  for (i in seq_along(nodes)) {
    for (j in seq_along(nodes)) {
      if (i == j) next
      if (is_subset(tipset[[i]], tipset[[j]])) {
        if (length(tipset[[i]]) != length(tipset[[j]])) { redundant[i] <- TRUE; break }
        else if (depth[i] > depth[j])                   { redundant[i] <- TRUE; break }
      }
    }
  }
  nodes[!redundant]
}

# --- base taxon set, matching 06B / 12D / 12E exactly ---
metadata <- read.delim(file.path(PHYLOGWAS_ROOT, "data/Poaceae_metadata_filtered_2025.08.28.tsv"), header = TRUE)
colnames(metadata)[3] <- "latest_name"
rownames(metadata) <- metadata$assemblyID

spTre <- read.tree(file.path(PHYLOGWAS_ROOT, "output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk"))
spTre.rooted <- root(spTre, "ASM1935983v1", resolve.root = TRUE)
spTre.rooted$edge.length[is.na(spTre.rooted$edge.length)] <- 0.1

ePCs <- readRDS(file.path(PHYLOGWAS_ROOT, "output/ePC_20250804.rds"))
baseCommon <- intersect(spTre.rooted$tip.label, rownames(ePCs$environmental.features))
baseCommon <- intersect(baseCommon, metadata$assemblyID)

# --- full-set ASR, computed ONCE; its transition clades are the fixed reference ---
tree.sub <- keep.tip(spTre.rooted, baseCommon)
tree.sub <- keep.tip(tree.sub, intersect(baseCommon, rownames(ePCs$synthetic.environmental.traits)))
tree.sub <- multi2di(tree.sub)

trait <- na.omit(ePCs$synthetic.environmental.traits[, ENVPC_IDX])
names(trait) <- rownames(ePCs$synthetic.environmental.traits)
trait <- trait[tree.sub$tip.label]

fit <- fastAnc(tree.sub, trait, var = TRUE, CI = TRUE)
d <- rbind(
  data.frame(node = 1:Ntip(tree.sub), ace = trait, var = 0, upper = trait, lower = trait),
  data.frame(node = (Ntip(tree.sub) + 1):(Ntip(tree.sub) + Nnode(tree.sub)),
             ace = fit$ace, var = fit$var, upper = fit$CI95[, 2], lower = fit$CI95[, 1]))
d$node <- as.numeric(d$node)
plotTree <- dplyr::full_join(tree.sub, d, by = "node")

edgeState      <- apply(plotTree@phylo$edge, c(1, 2), function(x) d[d$node == x, 2])
edgeStateUpper <- apply(plotTree@phylo$edge, c(1, 2), function(x) d[d$node == x, 4])
edgeStateLower <- apply(plotTree@phylo$edge, c(1, 2), function(x) d[d$node == x, 5])

toL <- edgeStateLower[, 1] > edgeStateUpper[, 2] & edgeState[, 2] < quantile(trait, .3) & edgeState[, 1] > quantile(trait, .3)
toH <- edgeStateLower[, 2] > edgeStateUpper[, 1] & edgeState[, 2] > quantile(trait, .7) & edgeState[, 1] < quantile(trait, .7)

nodetoH <- keep_most_ancestral(plotTree@phylo, plotTree@phylo$edge[toH, 2])
nodetoL <- keep_most_ancestral(plotTree@phylo, plotTree@phylo$edge[toL, 2])

tipsOf <- function(nds) lapply(nds, function(x) sort(tree.sub$tip.label[Descendants(tree.sub, x, type = "tips")[[1]]]))
clades    <- c(tipsOf(nodetoH), tipsOf(nodetoL))
direction <- c(rep("toH", length(nodetoH)), rep("toL", length(nodetoL)))

cat(sprintf("Fixed reference (envPC_%d, full set): %d tips, %d transition clades (%d toH + %d toL), %d transition-carrying tips\n\n",
            ENVPC_IDX, Ntip(tree.sub), length(clades), length(nodetoH), length(nodetoL),
            length(unique(unlist(clades)))))

# --- per-subtribe structural loss ---
tax <- read.csv(file.path(PHYLOGWAS_ROOT, "data/Poaceae_taxonomic_groups_subtribeAdded_SKH20250402.csv"),
                stringsAsFactors = FALSE)
tax <- tax[!duplicated(tax$assemblyID), ]
rownames(tax) <- tax$assemblyID

subtribeSizes <- table(tax$subtribe[tax$assemblyID %in% baseCommon & !is.na(tax$subtribe)])
subtribes <- names(sort(subtribeSizes[subtribeSizes >= MIN_SUBTRIBE_N], decreasing = TRUE))
cat(sprintf("Scoring %d subtribes with >= %d taxa in the base set\n\n", length(subtribes), MIN_SUBTRIBE_N))

res <- data.frame()
for (st in subtribes) {
  dropped <- intersect(rownames(tax)[!is.na(tax$subtribe) & tax$subtribe == st], baseCommon)
  keep    <- setdiff(baseCommon, dropped)

  before <- sapply(clades, length)
  after  <- sapply(clades, function(cl) length(intersect(cl, keep)))

  res <- rbind(res, data.frame(
    subtribe         = st,
    taxa_dropped     = length(dropped),
    clades_wiped     = sum(after == 0),                     # transition removed outright
    clades_partial   = sum(after > 0 & after < before),     # eroded but still present
    clades_untouched = sum(after == before),
    wiped_toH        = sum(after == 0 & direction == "toH"),
    wiped_toL        = sum(after == 0 & direction == "toL"),
    tips_lost        = sum(before - after),
    pct_tips_lost    = round(100 * sum(before - after) / sum(before), 1)))
}

# --- join the jackknife stability metrics, if that run has completed ---
jkPath <- file.path(PHYLOGWAS_ROOT, "output/revisionAnalyses/jackknife_n20_envPC1_summary.csv")
if (file.exists(jkPath) && ENVPC_IDX == 1) {
  jk  <- read.csv(jkPath)
  res <- merge(res, jk[, c("subtribe", "rho", "recall_p001", "recall_p05")], by = "subtribe", all.x = TRUE)
  res <- res[order(-res$clades_wiped, -res$tips_lost), ]
  cat("--- transition loss vs jackknife stability ---\n")
  print(res, row.names = FALSE)
  ok <- !is.na(res$rho)
  cat("\n--- correlation with rho ---\n")
  for (v in c("clades_wiped", "tips_lost", "pct_tips_lost", "taxa_dropped")) {
    cat(sprintf("  %-14s Pearson %+.3f   Spearman %+.3f\n", v,
                cor(res[[v]][ok], res$rho[ok]),
                cor(res[[v]][ok], res$rho[ok], method = "spearman")))
  }
} else {
  res <- res[order(-res$clades_wiped, -res$tips_lost), ]
  print(res, row.names = FALSE)
}

outPath <- file.path(PHYLOGWAS_ROOT,
                     sprintf("output/revisionAnalyses/jackknife_transitionLoss_fixedASR_envPC%d.csv", ENVPC_IDX))
write.csv(res, outPath, row.names = FALSE)
cat(sprintf("\nWrote %s\n", outPath))
