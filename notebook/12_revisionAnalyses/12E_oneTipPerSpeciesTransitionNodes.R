# Analysis #4 (revision-analyses.md): one-tip-per-species ASR transition-node
# stability check. Faithfully reproduces notebook/06_envirotyping/06B_visualizationEnvAdapt.ipynb's
# ASR + transition-identification logic (cells 0-3, 45-46, 49) -- same rooting, same
# CI-based toL/toH edge classification, same keep_most_ancestral() collapsing --
# restricted to the one-tip-per-species taxon set (output/revisionAnalyses/traitTables/
# envData_oneTipPerSpecies.txt's 569 taxa), for all 3 envPCs (unlike
# 12D_andropogoneaeTransitionNodes.R, which only needed envPC1 for the power
# simulation -- this script is purely for the stability check, no downstream
# power-sim consumer).
#
# Lightweight: ~569 tips, three traits' fastAnc() calls -- seconds, not hours.

suppressMessages({
  library(ape)
  library(phangorn)
  library(phytools)
  library(ggtree)   # for full_join.phylo, matching the original notebook exactly
  library(rjson)
})

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")

# --- keep_most_ancestral(): copied verbatim from 06B_visualizationEnvAdapt.ipynb cell 45 ---
keep_most_ancestral <- function(tree, nodes) {
  stopifnot(inherits(tree, "phylo"))
  nodes <- unique(as.integer(nodes))

  max_id <- Ntip(tree) + tree$Nnode
  nodes <- nodes[nodes >= 1 & nodes <= max_id]
  if (length(nodes) <= 1) return(nodes)

  depth <- node.depth.edgelength(tree)
  depth <- depth[nodes]
  names(depth) <- as.character(nodes)

  tipset <- lapply(nodes, function(nd) sort(Descendants(tree, nd, type = "tips")[[1]]))
  names(tipset) <- as.character(nodes)

  is_subset <- function(A, B) length(A) <= length(B) && all(A %in% B)

  redundant <- logical(length(nodes))
  for (i in seq_along(nodes)) {
    Ai <- tipset[[i]]; di <- depth[i]
    for (j in seq_along(nodes)) {
      if (i == j) next
      Bj <- tipset[[j]]; dj <- depth[j]
      if (is_subset(Ai, Bj)) {
        identical_sets <- length(Ai) == length(Bj)
        if (!identical_sets) { redundant[i] <- TRUE; break }
        else if (di > dj) { redundant[i] <- TRUE; break }
      }
    }
  }
  nodes[!redundant]
}

# --- Steps 0-3, matching 06B cells 0-3 exactly ---
metadata <- read.delim(file.path(PHYLOGWAS_ROOT, "data/Poaceae_metadata_filtered_2025.08.28.tsv"), header = TRUE)
metadata$spTaxa <- paste(metadata[, 1], metadata[, 3], sep = ":")
colnames(metadata)[3] <- "latest_name"
rownames(metadata) <- metadata$assemblyID

spTre <- read.tree(file.path(PHYLOGWAS_ROOT, "output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk"))
spTre.rooted <- root(spTre, "ASM1935983v1", resolve.root = TRUE)
spTre.rooted$edge.length[is.na(spTre.rooted$edge.length)] <- 0.1

ePCs <- readRDS(file.path(PHYLOGWAS_ROOT, "output/ePC_20250804.rds"))

commonID <- intersect(spTre.rooted$tip.label, rownames(ePCs$environmental.features))
commonID <- intersect(commonID, metadata$assemblyID)
spTre.rooted.filtered <- keep.tip(spTre.rooted, commonID)

# --- one-tip-per-species restriction (the only deviation from the original notebook) ---
# Rooting already happened on the full tree above, before any pruning, so (unlike
# 12_runPermulation_perOGModel.R) no outgroup needs to survive this subsetting -- but
# the taxon list itself already force-includes the outgroup via 12A's writeSubset().
fullID <- commonID   # full taxon set, for the fixed (full-set) percentile cutoffs below
oneTipTaxa <- rownames(read.table(file.path(PHYLOGWAS_ROOT, "output/revisionAnalyses/traitTables/envData_oneTipPerSpecies.txt"), header = TRUE))
commonID <- intersect(commonID, oneTipTaxa)
cat(sprintf("One-tip-per-species taxa with tree + ePC + metadata coverage: %d\n", length(commonID)))

rooted_tree <- keep.tip(spTre.rooted.filtered, commonID)

outDir <- file.path(PHYLOGWAS_ROOT, "output/revisionAnalyses/oneTipPerSpeciesASR")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Transition calls use the 30th/70th percentiles of the trait as cutoffs. Two versions
# for the one-tip set:
#   own   -- percentiles of the 569-taxon subset itself (the original behaviour)
#   fixed -- percentiles of the full taxon set, as in the full-set ASR, so that a change
#            in counts reflects the dropped tips rather than shifted cutoffs.
# The ASR (fastAnc) is the same for both; only the cutoffs differ. The full-set ASR is
# recomputed here too, so the clade-matching below needs no external inputs.
asr_calls <- function(tips, i, cuts) {
  tree.sub <- multi2di(keep.tip(spTre.rooted.filtered, tips))
  trait <- ePCs$synthetic.environmental.traits[tree.sub$tip.label, i]
  names(trait) <- tree.sub$tip.label
  fit <- fastAnc(tree.sub, trait, var = TRUE, CI = TRUE)
  d <- rbind(data.frame(node = 1:Ntip(tree.sub), ace = trait, var = 0, upper = trait, lower = trait),
             data.frame(node = (Ntip(tree.sub) + 1):(Ntip(tree.sub) + Nnode(tree.sub)),
                        ace = fit$ace, var = fit$var, upper = fit$CI95[, 2], lower = fit$CI95[, 1]))
  E <- tree.sub$edge
  st <- matrix(d$ace[E], ncol = 2); up <- matrix(d$upper[E], ncol = 2); lo <- matrix(d$lower[E], ncol = 2)
  out <- list(tree = tree.sub, trait = trait)
  for (cn in names(cuts)) {
    cl <- cuts[[cn]]
    toL <- lo[, 1] > up[, 2] & st[, 2] < cl[1] & st[, 1] > cl[1]
    toH <- lo[, 2] > up[, 1] & st[, 2] > cl[2] & st[, 1] < cl[2]
    nodes <- list(H = keep_most_ancestral(tree.sub, E[toH, 2]), L = keep_most_ancestral(tree.sub, E[toL, 2]))
    out[[cn]] <- lapply(nodes, function(ns) list(nodes = ns,
      tipsets = lapply(ns, function(nd) sort(tree.sub$tip.label[Descendants(tree.sub, nd, type = "tips")[[1]]]))))
  }
  out
}

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
best_j <- function(sets, against) vapply(sets, function(a)
  if (length(against)) max(vapply(against, jaccard, numeric(1), a = a)) else 0, numeric(1))

fullTips <- intersect(fullID, rownames(ePCs$synthetic.environmental.traits))
oneTips  <- intersect(commonID, rownames(ePCs$synthetic.environmental.traits))
results <- list(); counts <- list(); match <- list()
for (i in 1:3) {
  cat(sprintf("\n--- envPC_%d ---\n", i))
  cutFull <- quantile(ePCs$synthetic.environmental.traits[fullTips, i], c(.3, .7))
  cutOwn  <- quantile(ePCs$synthetic.environmental.traits[oneTips, i], c(.3, .7))
  full <- asr_calls(fullTips, i, list(full = cutFull))
  one  <- asr_calls(oneTips, i, list(own = cutOwn, fixed = cutFull))
  cat(sprintf("Full (%d tips): toH=%d toL=%d | one-tip (%d tips) own: toH=%d toL=%d, fixed: toH=%d toL=%d\n",
              Ntip(full$tree), length(full$full$H$nodes), length(full$full$L$nodes), Ntip(one$tree),
              length(one$own$H$nodes), length(one$own$L$nodes), length(one$fixed$H$nodes), length(one$fixed$L$nodes)))
  counts[[i]] <- data.frame(envPC = i, cut30_own = cutOwn[1], cut70_own = cutOwn[2],
    cut30_full = cutFull[1], cut70_full = cutFull[2],
    toH_full = length(full$full$H$nodes), toL_full = length(full$full$L$nodes),
    toH_own = length(one$own$H$nodes), toL_own = length(one$own$L$nodes),
    toH_fixed = length(one$fixed$H$nodes), toL_fixed = length(one$fixed$L$nodes))

  # Clade matching by shared tips, on the one-tip taxa. Recall direction: is each full-set
  # transition clade still called? A full clade whose tips are all duplicate-species
  # assemblies has no tips left after downsampling ("removed").
  for (dir in c("H", "L")) {
    fullSets <- lapply(full$full[[dir]]$tipsets, intersect, oneTips)
    removed <- lengths(fullSets) == 0
    for (cn in c("own", "fixed")) {
      oneSets <- one[[cn]][[dir]]$tipsets
      jr <- best_j(fullSets[!removed], oneSets)          # full -> one-tip (recall)
      jp <- best_j(oneSets, fullSets[!removed])          # one-tip -> full (precision)
      match[[length(match) + 1]] <- data.frame(envPC = i, direction = paste0("to", dir), cutoff = cn,
        full_clades = length(fullSets), removed_by_downsampling = sum(removed),
        recovered_J50 = sum(jr >= .5), recovered_any = sum(jr > 0), lost = sum(jr == 0),
        mean_J_full_to_onetip = round(mean(jr), 3),
        onetip_clades = length(oneSets), onetip_matched_J50 = sum(jp >= .5),
        mean_J_onetip_to_full = round(mean(jp), 3))
    }
  }
  results[[i]] <- list(full = full$full, own = one$own, fixed = one$fixed)
}

cc <- do.call(rbind, counts); mm <- do.call(rbind, match)
saveRDS(results, file.path(outDir, "oneTipPerSpecies_transitionNodes.rds"))
write.csv(cc, file.path(outDir, "oneTipPerSpecies_transitionCounts_ownVsFixedCutoff.csv"), row.names = FALSE)
write.csv(mm, file.path(outDir, "oneTipPerSpecies_cladeMatching.csv"), row.names = FALSE)
cat("\nTransition counts:\n"); print(cc, row.names = FALSE)
cat("\nClade matching (Jaccard on one-tip taxa):\n"); print(mm, row.names = FALSE)
write.tree(rooted_tree, file.path(outDir, "oneTipPerSpecies_tree.nwk"))
write(rjson::toJSON(lapply(results, function(x) list(x$own$H$nodes, x$own$L$nodes))),
      file = file.path(outDir, "oneTipPerSpecies_realTransition.json"))

cat(sprintf("\nWrote outputs to %s\n", outDir))
