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
oneTipTaxa <- rownames(read.table(file.path(PHYLOGWAS_ROOT, "output/revisionAnalyses/traitTables/envData_oneTipPerSpecies.txt"), header = TRUE))
commonID <- intersect(commonID, oneTipTaxa)
cat(sprintf("One-tip-per-species taxa with tree + ePC + metadata coverage: %d\n", length(commonID)))

rooted_tree <- keep.tip(spTre.rooted.filtered, commonID)

outDir <- file.path(PHYLOGWAS_ROOT, "output/revisionAnalyses/oneTipPerSpeciesASR")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

results <- list()
for (i in 1:3) {
  cat(sprintf("\n--- envPC_%d ---\n", i))
  tree.sub <- keep.tip(rooted_tree, intersect(commonID, rownames(ePCs$synthetic.environmental.traits)))
  tree.sub <- multi2di(tree.sub)
  cat(sprintf("Final ASR tree: %d tips\n", Ntip(tree.sub)))

  trait <- na.omit(ePCs$synthetic.environmental.traits[, i])
  names(trait) <- rownames(ePCs$synthetic.environmental.traits)
  trait <- trait[tree.sub$tip.label]

  fit <- fastAnc(tree.sub, trait, var = TRUE, CI = TRUE)

  td <- data.frame(node = 1:Ntip(tree.sub), ace = trait, var = 0, upper = trait, lower = trait)
  nd <- data.frame(node = (Ntip(tree.sub) + 1):(Ntip(tree.sub) + Nnode(tree.sub)),
                    ace = fit$ace, var = fit$var, upper = fit$CI95[, 2], lower = fit$CI95[, 1])
  d <- rbind(td, nd)
  d$node <- as.numeric(d$node)
  plotTree <- dplyr::full_join(tree.sub, d, by = "node")

  edgeState <- apply(plotTree@phylo$edge, c(1, 2), function(x) d[d$node == x, 2])
  edgeStateUpper <- apply(plotTree@phylo$edge, c(1, 2), function(x) d[d$node == x, 4])
  edgeStateLower <- apply(plotTree@phylo$edge, c(1, 2), function(x) d[d$node == x, 5])

  toL <- edgeStateLower[, 1] > edgeStateUpper[, 2] & edgeState[, 2] < quantile(trait, .3) & edgeState[, 1] > quantile(trait, .3)
  toH <- edgeStateLower[, 2] > edgeStateUpper[, 1] & edgeState[, 2] > quantile(trait, .7) & edgeState[, 1] < quantile(trait, .7)
  cat(sprintf("Raw transition edges: toH=%d toL=%d\n", sum(toH, na.rm = TRUE), sum(toL, na.rm = TRUE)))

  nodetoH <- keep_most_ancestral(plotTree@phylo, plotTree@phylo$edge[toH, 2])
  nodetoL <- keep_most_ancestral(plotTree@phylo, plotTree@phylo$edge[toL, 2])
  cat(sprintf("After keep_most_ancestral: toH=%d toL=%d\n", length(nodetoH), length(nodetoL)))

  # Tip-membership sets for each transition node, restricted to the common one-tip-per-species
  # taxon set -- used downstream to compare against the full-707-taxon calls by shared-tip
  # overlap rather than raw node number (node IDs aren't comparable across differently-pruned
  # trees).
  tipsetH <- lapply(nodetoH, function(nd) sort(tree.sub$tip.label[Descendants(tree.sub, nd, type = "tips")[[1]]]))
  tipsetL <- lapply(nodetoL, function(nd) sort(tree.sub$tip.label[Descendants(tree.sub, nd, type = "tips")[[1]]]))

  results[[i]] <- list(nodetoH = nodetoH, nodetoL = nodetoL, tipsetH = tipsetH, tipsetL = tipsetL)
}

saveRDS(results, file.path(outDir, "oneTipPerSpecies_transitionNodes.rds"))
write.tree(rooted_tree, file.path(outDir, "oneTipPerSpecies_tree.nwk"))
write(rjson::toJSON(lapply(results, function(x) list(x$nodetoH, x$nodetoL))),
      file = file.path(outDir, "oneTipPerSpecies_realTransition.json"))

cat(sprintf("\nWrote outputs to %s\n", outDir))
