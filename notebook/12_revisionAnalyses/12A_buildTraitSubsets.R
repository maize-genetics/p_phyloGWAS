# Builds the taxon-subset trait tables needed for the revision analyses (see
# revision-analyses.md): Andropogoneae-only, one-tip-per-species, and 24
# leave-one-subtribe-out subsets. Each output file is read.table(header=T)-
# compatible with 12_runPermulation_perOGModel.R's `dir.trait` argument --
# built via read.table/write.table on the original file, not hand-rolled, so
# the format matches exactly.
#
# Data-wrangling only, no model fitting -- safe to run regardless of system load.

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")

envData <- read.table(file.path(PHYLOGWAS_ROOT, "output/envData_707Poaceae_20250804.txt"), header = TRUE)
stopifnot(nrow(envData) == 707)

tax <- read.csv(file.path(PHYLOGWAS_ROOT, "data/Poaceae_taxonomic_groups_subtribeAdded_SKH20250402.csv"),
                 stringsAsFactors = FALSE)
# One exact duplicate row exists (Zm-B73-REFERENCE-NAM-5.0, identical content both times) --
# drop defensively so rownames<- doesn't error; verified this doesn't change any counts
# since a byte-identical duplicate carries no new information either way.
nBeforeDedup <- nrow(tax)
tax <- tax[!duplicated(tax$assemblyID), ]
if (nrow(tax) < nBeforeDedup) cat(sprintf("Dropped %d exact-duplicate taxonomy row(s)\n", nBeforeDedup - nrow(tax)))
rownames(tax) <- tax$assemblyID

qc <- read.delim(file.path(PHYLOGWAS_ROOT, "data/Poaceae_metadata_filtered_2025.08.28.tsv"),
                  stringsAsFactors = FALSE)
qc <- qc[!duplicated(qc$assemblyID), ]
rownames(qc) <- qc$assemblyID

allTaxa <- rownames(envData)
stopifnot(all(allTaxa %in% rownames(qc)))  # QC metrics must cover every model taxon

# AN21TNTL0178 (Sorghum arundinaceum) has no taxonomy-CSV row -- excluded by
# default from every taxonomy-driven subset below (see revision-analyses.md).
unclassified <- setdiff(allTaxa, rownames(tax))
cat("Unclassified taxa (excluded from taxonomy-driven subsets):", paste(unclassified, collapse = ", "), "\n")
classified <- setdiff(allTaxa, unclassified)

outDir <- file.path(PHYLOGWAS_ROOT, "output/revisionAnalyses/traitTables")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# 12_runPermulation_perOGModel.R prunes the tree to the trait table's taxa
# *before* rooting (keep.tip() then root(spTre, "ASM1935983v1")) -- so this one
# taxon (Pharus latifolius, the actual designated outgroup, clade="outgroup")
# must be present in every subset trait table or the script errors out at
# rooting, before any model ever fits. Confirmed via a real failed run, not
# assumed. It's `subtribe`-NA and (as far as this taxonomy table goes) the sole
# assembly for its species, so one-tip-per-species and all 24 subtribe-drops
# already keep it naturally -- only the Andropogoneae-only subset needs it
# force-included, adding a single non-Andropogoneae row (1/484) to that
# analysis's population. This mirrors the original full-707 analysis, which
# already includes this same outgroup taxon as part of its own population.
ROOT_OUTGROUP <- "ASM1935983v1"
stopifnot(ROOT_OUTGROUP %in% allTaxa)

writeSubset <- function(taxa, filename) {
  taxa <- union(taxa, ROOT_OUTGROUP)  # safety net: force-include on every subset, not just #1
  sub <- envData[taxa, , drop = FALSE]
  write.table(sub, file.path(outDir, filename), quote = FALSE)
  cat(sprintf("  %-45s %d taxa\n", filename, nrow(sub)))
}

manifest <- list()

## 1. Andropogoneae-only ------------------------------------------------------
andro <- classified[tax[classified, "tribe"] == "Andropogoneae"]
cat("\nAndropogoneae-only subset (+1 forced outgroup taxon for tree rooting):\n")
writeSubset(andro, "envData_andropogoneae.txt")
manifest[["andropogoneae"]] <- data.frame(subset = "andropogoneae", nTaxa = length(union(andro, ROOT_OUTGROUP)), droppedGroup = "non-Andropogoneae")

## 2. One-tip-per-species ------------------------------------------------------
# Best QC per duplicated species: highest BUSCO `complete`, then lowest `duplicated`.
species <- tax[classified, "Species"]
ord <- order(species, -qc[classified, "complete"], qc[classified, "duplicated"])
classifiedOrdered <- classified[ord]
speciesOrdered <- species[ord]
keepOneTip <- classifiedOrdered[!duplicated(speciesOrdered)]
cat("\nOne-tip-per-species subset:\n")
writeSubset(keepOneTip, "envData_oneTipPerSpecies.txt")
manifest[["oneTipPerSpecies"]] <- data.frame(subset = "oneTipPerSpecies", nTaxa = length(keepOneTip), droppedGroup = "duplicate-species tips")

## 3. Leave-one-subtribe-out (subtribes with >=5 taxa within the 707) --------
subtribeCounts <- table(tax[classified, "subtribe"])
subtribeCounts <- subtribeCounts[names(subtribeCounts) != "" & !is.na(names(subtribeCounts))]
dropSubtribes <- names(subtribeCounts[subtribeCounts >= 5])
cat(sprintf("\nLeave-one-subtribe-out: %d subtribes with >=5 taxa\n", length(dropSubtribes)))

for (st in dropSubtribes) {
  keep <- allTaxa[!(allTaxa %in% classified[tax[classified, "subtribe"] == st])]
  fname <- sprintf("envData_dropSubtribe_%s.txt", gsub("[^A-Za-z0-9]", "", st))
  writeSubset(keep, fname)
  manifest[[paste0("drop_", st)]] <- data.frame(subset = paste0("drop_", st), nTaxa = length(keep), droppedGroup = st)
}

manifestDf <- do.call(rbind, manifest)
write.csv(manifestDf, file.path(outDir, "manifest.csv"), row.names = FALSE)
cat(sprintf("\nWrote %d trait-table subsets + manifest.csv to %s\n", nrow(manifestDf), outDir))

# Safety net: every output file must contain ROOT_OUTGROUP or the model script
# will fail at tree-rooting before fitting a single OG.
for (f in list.files(outDir, pattern = "^envData_.*\\.txt$", full.names = TRUE)) {
  taxaInFile <- read.table(f, header = TRUE)
  stopifnot(ROOT_OUTGROUP %in% rownames(taxaInFile))
}
cat(sprintf("Verified: %s present in every output file.\n", ROOT_OUTGROUP))
