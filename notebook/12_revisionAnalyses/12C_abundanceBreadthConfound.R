# Analysis #3 (revision-analyses.md): does candidate status depend on protein
# abundance level + tissue-expression breadth (a technical/ascertainment-bias
# confound check)? Lightweight -- a handful of joins + glm() on ~18k rows, no
# heavy compute.

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")
FPKM_THRESHOLD <- as.numeric(Sys.getenv("FPKM_THRESHOLD", unset = "1"))

fpkm <- read.csv(file.path(PHYLOGWAS_ROOT, "data/maxFPKM_v5.csv"), stringsAsFactors = FALSE)

# maxFPKM_v5.csv has 2,271 genes with >1 row (same Gene, different maxFPKM/tissue
# values -- confirmed not simple duplicates, likely from different pan-genome
# comparison contexts upstream; taxaShort is "B73" for every row, so that's not
# the differentiator). Documented resolution: aggregate via max across duplicate
# rows per gene -- conservative, consistent with the column's own "maxFPKM"
# naming ("does this gene show any high-expression signal, in any tissue").
tissueCols <- colnames(fpkm)[23:45]
stopifnot(tissueCols[1] == "Female.Spikelet.Collected.on.day.as.silk",
          tissueCols[length(tissueCols)] == "Secondary.Root.7.8.Days")

geneAgg <- aggregate(fpkm[, c("maxFPKM", tissueCols)], by = list(Gene = fpkm$Gene), FUN = max, na.rm = TRUE)
geneAgg$tissue_breadth <- rowSums(geneAgg[, tissueCols] > FPKM_THRESHOLD, na.rm = TRUE)
geneAgg <- geneAgg[, c("Gene", "maxFPKM", "tissue_breadth")]
cat(sprintf("Aggregated %d maxFPKM_v5.csv rows to %d unique B73 genes (FPKM > %.1f threshold)\n",
            nrow(fpkm), nrow(geneAgg), FPKM_THRESHOLD))

# The 330 candidates: union of the three per-envPC candidate lists, each
# already filtered on emp_p <= 0.001 (see revision-analyses.md).
candFiles <- c("envPC1", "envPC2", "envPC3")
candList <- lapply(candFiles, function(pc) {
  readLines(file.path(PHYLOGWAS_ROOT, sprintf("output/candidateOG_%s.txt", pc)))
})
names(candList) <- candFiles
allCandidates <- unique(unlist(candList))
cat(sprintf("Candidate OGs: %d total (envPC1=%d, envPC2=%d, envPC3=%d, %d unique)\n",
            length(allCandidates), length(candList$envPC1), length(candList$envPC2),
            length(candList$envPC3), length(allCandidates)))
stopifnot(length(allCandidates) == 330)

# OG -> ZmID mapping, from the same homolog-added result files (envPC1's mapping
# is used as the OG universe; ZmID should be consistent OG-to-gene regardless of
# which envPC file it's read from, since it comes from the same OrthoFinder run).
ogRes <- read.delim(file.path(PHYLOGWAS_ROOT, "output/finalModels_20251002/envPC_1/ASREML_res_empPadded_homologAdded.txt"),
                     stringsAsFactors = FALSE)
ogZm <- ogRes[!is.na(ogRes$ZmID) & ogRes$ZmID != "", c("OG", "ZmID")]
ogZm$Gene <- sub("_T[0-9]+$", "", ogZm$ZmID)  # strip transcript suffix to match maxFPKM_v5.csv's Gene column
ogZm <- merge(ogZm, geneAgg, by = "Gene", all.x = FALSE)  # drop OGs whose maize ortholog has no FPKM data
cat(sprintf("OGs with both a ZmID and FPKM data: %d (of %d OGs with any ZmID)\n",
            nrow(ogZm), length(unique(ogRes$OG[!is.na(ogRes$ZmID) & ogRes$ZmID != ""]))))

ogZm$candidate_status <- as.integer(ogZm$OG %in% allCandidates)
for (pc in candFiles) ogZm[[paste0("candidate_", pc)]] <- as.integer(ogZm$OG %in% candList[[pc]])

outDir <- file.path(PHYLOGWAS_ROOT, "output/revisionAnalyses/abundanceBreadth")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
write.csv(ogZm, file.path(outDir, "og_abundance_breadth_table.csv"), row.names = FALSE)

fitAndReport <- function(data, label) {
  m <- glm(candidate_status ~ scale(maxFPKM) + scale(tissue_breadth), data = data, family = binomial)
  s <- summary(m)$coefficients
  cat(sprintf("\n=== %s (n=%d, %d candidates) ===\n", label, nrow(data), sum(data$candidate_status)))
  print(round(s, 4))
  s
}

results <- list()
results[["pooled_330"]] <- fitAndReport(within(ogZm, {}), "Pooled (330 vs. rest)")
for (pc in candFiles) {
  tmp <- ogZm
  tmp$candidate_status <- tmp[[paste0("candidate_", pc)]]
  results[[pc]] <- fitAndReport(tmp, pc)
}

resDf <- do.call(rbind, lapply(names(results), function(nm) {
  d <- as.data.frame(results[[nm]])
  d$term <- rownames(d)
  d$model <- nm
  d
}))
rownames(resDf) <- NULL
write.csv(resDf, file.path(outDir, "glm_results.csv"), row.names = FALSE)
cat(sprintf("\nWrote %s and %s\n",
            file.path(outDir, "og_abundance_breadth_table.csv"),
            file.path(outDir, "glm_results.csv")))
