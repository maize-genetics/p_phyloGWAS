# Revision analysis 12L, step 2: per-assembly summary of the contigs that carry
# miniprot gene hits (per-hit tables from 12L_hitContigLength.sh).
#
# Scope: rank-1 hits (the copy each assembly contributes first) of OGs tested in
# the envPC1-3 models. All 727 assemblies are summarised; in_model_taxa flags the
# 707 used in the models.
#
# Per hit:
#   expected span = median genomic span (introns included) of that OG's rank-1
#                   hits in long-read assemblies, restricted to full-length
#                   (q_cov >= 0.9), gap-free hits
#   holds_gene    = host contig length >= expected span
#   edge          = hit ends within EDGE bp of a contig end
#   edge_trunc    = edge and q_cov < COV (cut off by the assembly break)
#
# Run from repo root:
#   /home/sh2246/miniconda3/envs/phyloGWAS-r/bin/Rscript notebook/12_revisionAnalyses/12L_hitContigLength_summary.R
suppressMessages({library(data.table); library(ggplot2)})

EDGE <- 50
COV  <- 0.9
out_dir <- "output/hitContigLength"
fig_dir <- file.path(out_dir, "figure"); dir.create(fig_dir, showWarnings = FALSE)

# ---- inputs ----
model_OG <- unique(unlist(lapply(1:3, function(k)
  fread(sprintf("output/finalModels_20251002/envPC_%d/ASREML_res.txt", k), select = "OG")$OG)))
model_taxa <- rownames(read.table("output/envData_707Poaceae_20250804.txt", header = TRUE))

meta <- fread("data/Poaceae_metadata_filtered_2025.08.28.tsv", select = c("assemblyID", "technology", "source"))
meta[, class := fifelse(grepl("pacbio|nanopore|\\bONT\\b|promethion|minion", technology, ignore.case = TRUE),
                        "long-read", "short-read")]

a12k <- fread("output/assemblyStats/assemblyStats_727.tsv", select = c("assemblyID", "contig_N50"))

files <- list.files(file.path(out_dir, "perHit"), pattern = "\\.tsv\\.gz$", full.names = TRUE)
hits <- rbindlist(lapply(files, function(f) {
  d <- fread(cmd = paste("zcat", shQuote(f)), showProgress = FALSE)
  d[rank == 1 & OG %in% model_OG]
}))
hits <- merge(hits, meta[, .(assemblyID, class)], by = "assemblyID", all.x = TRUE)
cat(sprintf("%d assemblies, %s rank-1 hits on %d model OGs\n",
            uniqueN(hits$assemblyID), format(nrow(hits), big.mark = ","), uniqueN(hits$OG)))

# ---- expected gene span per OG, from long-read assemblies ----
exp_span <- hits[class == "long-read" & q_cov >= COV & spans_gap == 0,
                 .(exp_span = as.numeric(median(hit_span)), n_ref = .N), by = OG]
hits <- merge(hits, exp_span, by = "OG", all.x = TRUE)
hits[, `:=`(holds_gene = contig_len >= exp_span,
            edge       = pmin(dist_left, dist_right) < EDGE)]
hits[, edge_trunc := edge & q_cov < COV]

# ---- per-assembly summary ----
summ <- hits[, .(
  n_hits                = .N,
  median_hitContig_len  = as.numeric(median(contig_len)),
  pct_holds_gene        = round(100 * mean(holds_gene, na.rm = TRUE), 2),
  pct_edge              = round(100 * mean(edge), 2),
  pct_edge_truncated    = round(100 * mean(edge_trunc), 2),
  pct_spans_gap         = round(100 * mean(spans_gap), 2),
  median_q_cov          = round(median(q_cov), 4),
  pct_q_cov_ge90        = round(100 * mean(q_cov >= COV), 2)
), by = .(assemblyID, class)]
summ <- merge(summ, a12k, by = "assemblyID", all.x = TRUE)
summ[, in_model_taxa := assemblyID %in% model_taxa]
setorder(summ, assemblyID)
fwrite(summ, file.path(out_dir, "hitContigLength_727.tsv"), sep = "\t")
fwrite(exp_span, file.path(out_dir, "expectedSpan_perOG.tsv"), sep = "\t")

# ---- console summary (model taxa) ----
m <- summ[in_model_taxa == TRUE]
q <- function(x) sprintf("%s [%s-%s]", format(median(x), big.mark = ","),
                         format(quantile(x, .1), big.mark = ","), format(quantile(x, .9), big.mark = ","))
cat("\nPer-assembly medians [10th-90th pct], 707 model taxa\n")
print(m[, .(n = .N,
            contig_N50           = q(contig_N50),
            hitContig_len        = q(median_hitContig_len),
            pct_holds_gene       = q(pct_holds_gene),
            pct_edge             = q(pct_edge),
            pct_edge_truncated   = q(pct_edge_truncated),
            pct_q_cov_ge90       = q(pct_q_cov_ge90)), by = class])
cat("\nPooled rank-1 hits, 707 model taxa\n")
print(hits[assemblyID %in% model_taxa, .(
  hits = .N,
  pct_holds_gene     = round(100 * mean(holds_gene, na.rm = TRUE), 1),
  pct_edge           = round(100 * mean(edge), 1),
  pct_edge_truncated = round(100 * mean(edge_trunc), 1),
  pct_q_cov_ge90     = round(100 * mean(q_cov >= COV), 1)), by = class])

cat("\nExpected gene span per OG (bp), from long-read assemblies\n")
print(round(quantile(exp_span$exp_span, c(.05, .1, .25, .5, .75, .9, .95, .99))))

# Is per-assembly truncation a real gene-space quality signal? Check it
# against BUSCO fragmentation, which is computed independently.
busco <- fread("data/Poaceae_metadata_filtered_2025.08.28.tsv", select = c("assemblyID", "fraged"))
sr <- merge(m[class == "short-read"], busco, by = "assemblyID")
cat("\nShort-read assemblies: Spearman rho of pct_edge_truncated with\n")
print(round(sapply(c(BUSCO_fragmented = "fraged", contig_N50 = "contig_N50",
                     median_hitContig_len = "median_hitContig_len"),
                   function(v) cor(sr$pct_edge_truncated, as.numeric(sr[[v]]),
                                   method = "spearman", use = "complete.obs")), 2))

# Truncation should concentrate in long, intron-rich genes if contig breaks
# (not missing gene space) are what's behind it.
span_bins <- c(0, 1e3, 2e3, 4e3, 8e3, 16e3, Inf)
byspan <- hits[assemblyID %in% model_taxa & !is.na(exp_span),
               .(hits = .N,
                 pct_holds_gene     = round(100 * mean(holds_gene), 1),
                 pct_edge_truncated = round(100 * mean(edge_trunc), 1)),
               by = .(class, span_bin = cut(exp_span, span_bins, dig.lab = 6))][order(class, span_bin)]
cat("\nPooled rank-1 hits by expected gene span, 707 model taxa\n")
print(byspan)
fwrite(byspan, file.path(out_dir, "hitContigLength_byGeneSpan.tsv"), sep = "\t")

# ---- figures (707 model taxa) ----
theme_set(theme_bw(base_size = 12))
cols <- c("long-read" = "#1b6ca8", "short-read" = "#d0892b")
p1 <- ggplot(m, aes(contig_N50, median_hitContig_len, colour = class)) +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_point(alpha = .6, size = 1.4) +
  scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_si("b"))) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_si("b"))) +
  scale_colour_manual(values = cols) +
  labs(x = "Contig N50 (all contigs)", y = "Median length of contigs carrying a gene hit",
       colour = NULL, title = "Gene-bearing contigs vs. assembly-wide contig N50",
       subtitle = "One point per assembly; dashed line = 1:1") +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "12L_hitContig_vs_contigN50.png"), p1, width = 7, height = 6.5, dpi = 300)

long <- melt(m, id.vars = c("assemblyID", "class"),
             measure.vars = c("pct_holds_gene", "pct_q_cov_ge90", "pct_edge", "pct_edge_truncated"))
long[, variable := factor(variable,
  levels = c("pct_holds_gene", "pct_q_cov_ge90", "pct_edge", "pct_edge_truncated"),
  labels = c("Host contig >= expected\ngene span", "Hit covers >= 90%\nof query",
             sprintf("Hit within %d bp\nof contig end", EDGE), "Edge-truncated\n(edge & cov < 90%)"))]
p2 <- ggplot(long, aes(class, value, fill = class)) +
  geom_boxplot(outlier.shape = NA, alpha = .5, width = .6) +
  geom_jitter(width = .15, size = .5, alpha = .35) +
  facet_wrap(~variable, nrow = 1) +
  scale_fill_manual(values = cols, guide = "none") +
  labs(x = NULL, y = "% of rank-1 hits per assembly", title = "Gene-space completeness per assembly")
ggsave(file.path(fig_dir, "12L_hitContig_metrics.png"), p2, width = 11, height = 5, dpi = 300)

p3 <- ggplot(m, aes(pct_edge_truncated, fill = class)) +
  geom_histogram(binwidth = 1, alpha = .8, position = "identity") +
  scale_fill_manual(values = cols) +
  labs(x = sprintf("%% of rank-1 hits edge-truncated (within %d bp of a contig end and < %d%% query coverage)",
                   EDGE, 100 * COV),
       y = "Assemblies", fill = NULL, title = "Edge truncation per assembly",
       subtitle = sprintf("%d model taxa", nrow(m))) +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "12L_edgeTruncated_hist.png"), p3, width = 8, height = 5, dpi = 300)

p4 <- ggplot(exp_span, aes(exp_span)) +
  geom_histogram(bins = 80, fill = cols[["long-read"]], alpha = .8) +
  geom_vline(xintercept = median(exp_span$exp_span), linetype = 2) +
  scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_si("b"))) +
  labs(x = "Expected genomic span per OG (median of full-length long-read hits, introns included)",
       y = "OGs", title = "Gene length distribution",
       subtitle = sprintf("%s model OGs; dashed line = median %s bp", format(nrow(exp_span), big.mark = ","),
                          format(round(median(exp_span$exp_span)), big.mark = ",")))
ggsave(file.path(fig_dir, "12L_expectedSpan_hist.png"), p4, width = 8, height = 5, dpi = 300)

p5 <- ggplot(byspan, aes(span_bin, pct_edge_truncated, fill = class)) +
  geom_col(position = position_dodge(.8), width = .75) +
  scale_fill_manual(values = cols) +
  labs(x = "Expected gene span (bp)", y = "% of rank-1 hits edge-truncated", fill = NULL,
       title = "Edge truncation by gene length", subtitle = "Pooled hits, 707 model taxa") +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "12L_edgeTruncated_byGeneSpan.png"), p5, width = 8, height = 5, dpi = 300)
cat("\nwrote", file.path(out_dir, "hitContigLength_727.tsv"), "and figures in", fig_dir, "\n")
