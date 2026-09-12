# Revision analyses

Plan for the five robustness analyses requested in review, before the next round of revision.
Work happens on the `revision-analyses` branch; `main` (tagged `v0.1.0`, the reviewed submission)
stays untouched until this is merged back the same way `refactor-cleanup` was.

## Key numbers, verified directly against the actual data (not assumed)

The base taxon population for every analysis below is **`output/envData_707Poaceae_20250804.txt`'s
707 taxa** — this is literally what `rownames(envDf)` means in `src/12_runPermulation_perOGModel.R`,
i.e. the exact set the original stage-08 models were fit on. The taxonomy file
(`data/Poaceae_taxonomic_groups_subtribeAdded_SKH20250402.csv`) has 805 rows, a superset; every
count below is restricted to the 707 that actually matter.

| Quantity | Value | Source |
|---|---|---|
| "330 top candidates" (full-set analysis) | 330 unique OGs | Union of `output/candidateOG_envPC{1,2,3}.txt` (160+109+68), filtered on empirical permutation `emp_p ≤ 0.001` in `notebook/11_candidateOGInvestigation/11A_candidateOGInvestigation.ipynb` |
| Andropogoneae taxa | 483 of 707 (484 as actually run, see below) | `tribe` column, taxonomy CSV, restricted to the 707 |
| Subtribes with ≥5 taxa | 24 (of 51 total) | `subtribe` column, restricted to the 707 |
| Unique species | 569 of 707 | `Species` column, restricted to the 707; 87 species have >1 assembly (224 assemblies total) |

One taxon in the 707, **`AN21TNTL0178`** (*Sorghum arundinaceum*), has no row in the taxonomy CSV
at all — it's excluded by default from every taxonomy-driven subset below (Andropogoneae, subtribe
jackknife, one-tip-per-species), which is exactly what produces the 569 figure above. It can be
added back with a one-row taxonomy-CSV fix later if desired (would make it 570 for the
one-tip-per-species analysis, and +1 Andropogoneae/Sorghinae tip for #1/#2).

**Rooting constraint (found via a real failed run, not assumed):**
`12_runPermulation_perOGModel.R` prunes the tree to the trait table's taxa *before* rooting
(`keep.tip()` then `root(spTre, "ASM1935983v1")`) — so `ASM1935983v1` (*Pharus latifolius*, the
actual designated outgroup, `clade="outgroup"`) must be present in every subset trait table, or the
script errors out at rooting before fitting a single OG. It's `subtribe`-NA and the sole assembly
for its species, so one-tip-per-species (569) and all 24 subtribe-drops already keep it naturally
— only the **Andropogoneae-only subset needs it force-included, making that population 484 taxa,
not 483** (mirrors the original full-707 analysis, which already includes this same outgroup taxon
as part of its own population — a negligible 1/484 addition, not a real contamination of the
intended comparison).

## Why this is mostly a data-subsetting exercise, not new modeling code

`src/12_runPermulation_perOGModel.R` (the per-OG ASReml model) derives its taxon set entirely from
the trait table passed as its first CLI argument, then internally prunes the (already-rooted)
species tree and master data table to match. **Analyses #1, #4's association half, and #5 need no
changes to this script** — only new, taxon-filtered trait tables and output directories, via its
existing 8-argument CLI interface (`dir.trait, dir.masterTab, dir.phyloK, dir.tree, responseVar,
dir.out, permulation, GLM`). The ASR step (`src/S05_envPC_analysis.R`) is equally subset-agnostic.
Only the power simulation needs a small code change, and only the abundance/breadth analysis needs
genuinely new data.

---

## 1. Andropogoneae-only association model

Build a trait table = `envData_707Poaceae_20250804.txt` filtered to the 483 Andropogoneae
`assemblyID`s, plus the forced rooting outgroup (484 total, see above). Rerun
`12_runPermulation_perOGModel.R` unchanged (new `dir.trait`, new `dir.out`) for envPC_1/2/3,
permulation on.

**Comparison to the full-set run:**
- Spearman rank correlation of `-log10(raw p)` across all OGs tested in both runs (join on `OG`;
  raw `p`, not `emp_p`, since `emp_p` isn't computed for most OGs).
- Concordance of the 330: what fraction are also `emp_p ≤ 0.001` candidates in the Andropogoneae-only run.

## 2. Power simulation on the Andropogoneae tree

Prune `output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk` to the 484-taxon Andropogoneae
subset (483 + rooting outgroup); rerun the ASR step on that pruned tree/trait subset to get an
Andropogoneae-specific transition-node set (the equivalent of `output/powerSimulation/realTransition.json`).

Small code change to `src/powerSimulation_XY_revised.R`: make `TREE_PATH`/`KMAT_PATH`/`TRANS_PATH`
overridable via `Sys.getenv()`, defaulting to the current hardcoded values (same pattern already
used for `PHYLOGWAS_ROOT`) — keeps `08E_power_simulation.sh`'s existing invocation working
unchanged while allowing an env-var-driven Andropogoneae rerun. Same simulation config otherwise
(effect sizes, `n_sims=500`, scenarios), so it's a like-for-like power-vs-sample-size comparison.

## 3. Abundance/breadth confound check

For each maize ortholog in `data/maxFPKM_v5.csv`: pull `maxFPKM` as the abundance predictor, and
compute tissue-expression breadth = count of per-tissue FPKM columns clearing a detectability
threshold (proposed FPKM > 1; adjust if a different convention is preferred). Join to OG via the
`ZmID` columns already in `output/finalModels_20251002/envPC_*/ASREML_res_empPadded_homologAdded.txt`.

Fit `candidate_status (0/1, in the 330) ~ maxFPKM + tissue_breadth` (logistic regression), per
envPC and/or pooled. Report coefficients, p-values, and a plot.

## 4. One-tip-per-species downsampling

Within the 707: 87 species have multiple assemblies (224 assemblies total among them — e.g. Zea
mays, Triticum aestivum, Hordeum vulgare at 5 each). For each, keep the assembly ranking best under
the existing high-error-filter QC criteria (the same logic behind
`data/Poaceae_metadata_highErrorFiltered_2025.10.08.tsv`); the 483 single-assembly species pass
through unchanged; `AN21TNTL0178` excluded per the default above. **Result: 569 taxa, 138 dropped.**

- Rerun the ASR step on this pruned set — check whether transition-node calls stay stable versus
  full taxon sampling.
- Rerun `12_runPermulation_perOGModel.R` unchanged on this taxon subset for envPC_1/2/3 — same
  rank-correlation + 330-concordance comparison as analysis #1.

## 5. Leave-one-subtribe-out jackknife

**Scope reduced after real timing data** (see below): loop over the **top 8 of the 24** subtribes
with ≥5 taxa within the 707, ranked by taxon count:

Andropogoninae (150), Anthistiriinae (79), Oryzinae (35), Tripsacinae (33), Saccharinae (33),
Triticinae (28), Aveninae (24), Ischaeminae (23).

This naturally covers Andropogoneae's own largest subtribes (Andropogoninae, Anthistiriinae,
Tripsacinae, Saccharinae, Ischaeminae) plus the largest non-Andropogoneae lineages (Oryzinae,
Triticinae, Aveninae) — the drops with the most statistical leverage; the remaining 16 smaller
subtribes (Sorghinae down to Germainiinae, 5-23 taxa each) are deferred, not abandoned — can be
run later if more compute becomes available. The other 16:  Sorghinae (23), Apludinae (22),
Rottboelliinae (19), Ratzeburgiinae (17), Chrysopogoninae (16), Hordeinae (14), Rhytachninae (14),
Eleusininae (12), Cenchrinae (12), Arundinariinae (10), Arthraxoninae (10), Boivinellinae (8),
Panicinae (8), Brachypodiinae (7), Loliinae (5), Germainiinae (5).

For each, build a trait table excluding that subtribe's assemblies (everyone else — including the
39 NA-subtribe taxa and `AN21TNTL0178` — stays in), rerun `12_runPermulation_perOGModel.R`
unchanged for envPC_1/2/3.

**Per dropped subtribe**, report the same two stability metrics as #1 (rank correlation of
`-log10(p)`, and 330-candidate concordance), summarized across the 8 drops (table or plot of
concordance vs. subtribe size). This shows whether stability degrades for any particular lineage's
removal, and — since Andropogoneae's own subtribes are among the 8 — directly complements analysis
#1 from the opposite direction (drop-Andropogoneae's-subtribes vs. keep-only-Andropogoneae).

**Compute cost (verified, not estimated)**: a single envPC's *base model only* (no permutation) on
the 484-taxon Andropogoneae subset took **~12.5 hours at an 8-core budget** (20,634 OGs fit;
timestamps Sep 11 17:25 → Sep 12 06:00). The full 24-subtribe × 3-envPC = 72-run version of this
analysis would be ~900+ hours (5+ weeks) sequential at 8 cores — not feasible on a machine shared
with other running jobs. Reduced to top-8 × 3 = 24 runs (~same order as #1+#4 combined) as a result.
8-core budget chosen because a concurrent unrelated job (`p_lifeHistoryPoaceae`'s HyPhy RELAX run)
was already using 30 of this machine's 40 cores; revisit the core budget once that frees up.

---

## Open items before implementation starts

- Confirm the FPKM > 1 detectability threshold for tissue-expression breadth (#3), or supply a
  preferred alternative.
- Decide whether `AN21TNTL0178` should be added to the taxonomy CSV (with its correct
  tribe/subtribe) before these analyses run, or left excluded as the default above.
