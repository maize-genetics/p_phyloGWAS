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

**Deviation from plan:** permulation was run for envPC_1 initially, but convergence on the
183-taxon-smaller subset turned out to be far too low to be useful — the `topOG_list` OGs (already
selected for having a strong signal under the *real* trait) converge on the order of 3-6% of the
time per permuted-trait iteration (vs ~100% on the original full-707 runs), so completing even one
envPC's 1000-iteration permulation loop would take 35+ hours and still leave most OGs `NA` under a
sane minimum-converged-permutations threshold. Per author's direction, replaced the emp_p-based
330-concordance check with a raw-`p<0.001` threshold comparison (a looser, uncorrected proxy) —
still informative, just not the same statistic. This also exposed and fixed a real, generic bug:
`emp_p`'s old denominator was the fixed permutation count (1000) rather than the number that
actually converged, which would have silently deflated `emp_p` (false significance) whenever
permulation convergence was low — see `12_runPermulation_perOGModel.R`'s fix commit. All 3 envPCs'
base-model fits (permulation off) ran cleanly: 94.3% / 97.4% / 96.5% converged.

**Results (raw `p`, no permulation):**

| envPC | OGs compared (both converged) | Spearman ρ(-log10 p) | full p<0.001 → subset p<0.001 (recall / precision) | full p<0.001 → subset p<0.05 (recall) |
|---|---|---|---|---|
| envPC_1 | 18,075 | 0.759 | 92/275 (33.5%) / 92/151 (60.9%) | 252/275 (91.6%) |
| envPC_2 | 19,044 | 0.792 | 82/204 (40.2%) / 82/146 (56.2%) | 196/204 (96.1%) |
| envPC_3 | 18,846 | 0.790 | 88/209 (42.1%) / 88/163 (54.0%) | 205/209 (98.1%) |

The low p<0.001→p<0.001 recall is mostly a **power** effect, not signal disagreement: 92-98% of the
full-set's significant OGs are still at least nominally significant (p<0.05) in the Andropogoneae
subset. Moderately strong, consistent rank correlation (ρ≈0.76-0.79) across all three envPCs.

## 2. Power simulation on the Andropogoneae tree

Prune `output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk` to the 484-taxon Andropogoneae
subset (483 + rooting outgroup); rerun the ASR step on that pruned tree/trait subset to get an
Andropogoneae-specific transition-node set (the equivalent of `output/powerSimulation/realTransition.json`).

Small code change to `src/powerSimulation_XY_revised.R`: make `TREE_PATH`/`KMAT_PATH`/`TRANS_PATH`
overridable via `Sys.getenv()`, defaulting to the current hardcoded values (same pattern already
used for `PHYLOGWAS_ROOT`) — keeps `08E_power_simulation.sh`'s existing invocation working
unchanged while allowing an env-var-driven Andropogoneae rerun. Same simulation config otherwise
(effect sizes, `n_sims=500`, scenarios), so it's a like-for-like power-vs-sample-size comparison.

**ASR step** (`12D_andropogoneaeTransitionNodes.R`): the Andropogoneae-restricted reconstruction
yields **24 toH + 18 toL = 42 transition nodes**, against **27 + 19 = 46** on the full 707-taxon
tree. Only envPC1 is computed, since `powerSimulation_XY_revised.R` reads only
`fromJSON(...)[[1]]`. Outputs: `output/revisionAnalyses/powerSimulation/andropogoneae_realTransition.json`
and `andropogoneae_tree.nwk` (the 484-tip pruned tree).

**Results** — power difference (Andropogoneae − full) at the 100% convergence level, across the
full β range. Full run: `output/figure/power_df_merged.csv`; Andropogoneae run:
`output/revisionAnalyses/powerSimulation/power_df_merged.csv` (176 rows each = 4 scenarios × 4
convergence levels × 11 β).

| Scenario | max power loss | max power gain | β range where curves differ at all |
|---|---|---|---|
| empirical transition | −0.036 (β=0.05) | +0.016 (β=0.00) | both saturate at 1.000 from β≥0.30 |
| random transition − 1x | −0.032 (β=0.15) | +0.010 (β=0.25) | both ≥0.99 from β≥0.25 |
| random transition − 2x | −0.054 (β=0.10) | +0.010 (β=0.00) | both ≥0.98 from β≥0.20 |
| random transition − 1/2x | −0.006 (β=0.00) | **+0.092** (β=0.20) | Andropogoneae consistently *higher* |

Type I error (β=0) stayed near nominal in both runs across all four scenarios: 0.038–0.066.

**Interpretation**: restricting the tree to Andropogoneae costs **essentially no power**. Every
difference in the empirical-transition scenario and the 1x/2x controls sits inside Monte-Carlo
noise at `n_sims=500` (≈±0.04 near power 0.5). The 1/2x control is the one systematic departure,
and it runs in the *opposite* direction to a sample-size penalty — Andropogoneae is up to +0.09
higher. That is mechanistically consistent with the transition counts above: power in this design
tracks the number of independent transition events (42 vs 46, a 9% difference), not the number of
tips (484 vs 707, a 32% difference). The 1/2x scenario draws `n_trans × 0.5` nodes — 21 for
Andropogoneae vs 23 for the full tree — into a tree with far fewer deep, tip-rich clades, which
makes the drawn clades more evenly sized and the resulting contrast slightly cleaner.

## 3. Abundance/breadth confound check

For each maize ortholog in `data/maxFPKM_v5.csv`: pull `maxFPKM` as the abundance predictor, and
compute tissue-expression breadth = count of per-tissue FPKM columns clearing a detectability
threshold (proposed FPKM > 1; adjust if a different convention is preferred). Join to OG via the
`ZmID` columns already in `output/finalModels_20251002/envPC_*/ASREML_res_empPadded_homologAdded.txt`.

Fit `candidate_status (0/1, in the 330) ~ maxFPKM + tissue_breadth` (logistic regression), per
envPC and/or pooled. Report coefficients, p-values, and a plot.

**As run** (`12C_abundanceBreadthConfound.R`): breadth = count of per-tissue FPKM columns above
`FPKM_THRESHOLD` (default **1**, env-overridable), giving a 0–23 tissue count (median 21). After
joining `maxFPKM_v5.csv` to the OG→ZmID mapping, **16,773 OG–gene pairs** carry expression data,
of which **288 of the 330 candidates** are recoverable (147 envPC1 / 88 envPC2 / 60 envPC3 — the
shortfall is candidates with no maize ortholog or no expression record). Both predictors are
z-scaled, so coefficients are per-SD on the logit scale. Outputs:
`output/revisionAnalyses/abundanceBreadth/{og_abundance_breadth_table.csv,glm_results.csv}`.

**Results**:

| Model | `scale(maxFPKM)` β (p) | `scale(tissue_breadth)` β (p) |
|---|---|---|
| pooled (all 330) | −0.045 (p=0.705) | +0.023 (p=0.696) |
| envPC1 | −0.115 (p=0.612) | **+0.181 (p=0.048)** |
| envPC2 | +0.010 (p=0.897) | +0.010 (p=0.928) |
| envPC3 | −0.206 (p=0.629) | **−0.270 (p=0.021)** |

**Interpretation**: no evidence that candidate status is confounded by protein abundance or
expression breadth. `maxFPKM` is null everywhere (all p ≥ 0.61, including pooled). Breadth is null
pooled and in envPC2, and marginal in envPC1 and envPC3 — but those two point in **opposite
directions** (envPC1 positive, envPC3 negative). Across 8 tests, two hits at p=0.048 and p=0.021
is roughly what noise produces (neither survives any multiple-testing correction across the 8),
and a genuine detectability/abundance confound would have to push candidate status consistently one
way, not split by sign across envPCs. The pooled model — the direct test of the reviewer's
question — is unambiguously null on both predictors.

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

**Results — ASR transition-node stability** (`12E_oneTipPerSpeciesTransitionNodes.R`; permulation
dropped here too, for the same reason as analysis #1 above):

Transition-clade counts (toH / toL), full-707 vs one-tip-per-species (569):

| envPC | Full (707) | One-tip (569) |
|---|---|---|
| envPC_1 | 27 / 19 | 22 / 20 |
| envPC_2 | 21 / 22 | 18 / 22 |
| envPC_3 | 27 / 30 | 24 / 31 |

Tip-membership stability (each one-tip clade matched to its best-overlapping full-set clade by
Jaccard similarity, restricted to the shared 569-taxon overlap):

| envPC | toH mean Jaccard (≥0.5 match rate) | toL mean Jaccard (≥0.5 match rate) |
|---|---|---|
| envPC_1 | 0.955 (21/22) | 0.750 (15/20) |
| envPC_2 | 0.944 (17/18) | 0.864 (19/22) |
| envPC_3 | 0.878 (21/24) | 0.871 (27/31) |

Transition-node calls are stable under this downsampling: counts are close across the board, and
75-95.5% of one-tip-per-species clades substantially overlap (Jaccard ≥ 0.5) a full-set clade in
every envPC × direction combination — expected, since dropping only redundant duplicate-species
assemblies doesn't change the underlying species-level tree topology ASR depends on.

**Results — association model** (raw `p`, no permulation, same deviation as analysis #1):

| envPC | OGs compared (both converged) | Spearman ρ(-log10 p) | full p<0.001 → subset p<0.001 (recall / precision) | full p<0.001 → subset p<0.05 (recall) |
|---|---|---|---|---|
| envPC_1 | 18,099 | 0.793 | 100/275 (36.4%) / 100/172 (58.1%) | 266/275 (96.7%) |
| envPC_2 | 19,095 | 0.914 | 78/204 (38.2%) / 78/92 (84.8%) | 203/204 (99.5%) |
| envPC_3 | 18,981 | 0.861 | 78/211 (37.0%) / 78/107 (72.9%) | 210/211 (99.5%) |

Noticeably stronger than the Andropogoneae comparison (ρ 0.79-0.91 vs 0.76-0.79) — expected, since
one-tip-per-species keeps 569 taxa spanning the whole tree rather than restricting to one clade,
preserving far more phylogenetic breadth. 97-99.5% of the full-set's p<0.001 hits remain at least
nominally significant (p<0.05) here. Base-model fits: 96.6% / 98.5% / 98.4% converged.

## 5. Leave-one-subtribe-out jackknife

**Final scope**: the **10 subtribes with >20 taxa** within the 707, **envPC_1 only** (not all 3),
base-model fit with permulation off (same rationale as #1/#4: see the emp_p fix commit — permulation
convergence collapses on any taxon-reduced subset, and this analysis's question — does dropping one
lineage destabilize the raw association ranking — doesn't need it):

Andropogoninae (150), Anthistiriinae (79), Oryzinae (35), Tripsacinae (33), Saccharinae (33),
Triticinae (28), Aveninae (24), Ischaeminae (23), Sorghinae (23), Apludinae (22).

This covers Andropogoneae's own largest subtribes (Andropogoninae, Anthistiriinae, Tripsacinae,
Saccharinae, Ischaeminae, Sorghinae) plus the largest non-Andropogoneae lineages (Oryzinae,
Triticinae, Aveninae, Apludinae). The remaining 14 subtribes with 5-20 taxa are deferred, not
abandoned. (Scope evolved during the session: originally top-8 by count under a permulation-
inclusive, 8-core estimate of 5+ weeks; that estimate turned out to be based on the wrong bottleneck
-- see below -- so the real constraint is much lighter and the cutoff was raised to >20 taxa, adding
Sorghinae and Apludinae, at 32 cores and permulation off.)

For each, built a trait table excluding that subtribe's assemblies (everyone else — including the
39 NA-subtribe taxa and `AN21TNTL0178` — stays in), ran `12_runPermulation_perOGModel.R` (base model
only, `12B_runSubsetModels.sh`'s new `PERMULATION=FALSE ENVPCS="1"` overrides).

**Compute cost (verified, not estimated)**: at 32 cores, permulation off, each run took ~3.0-3.5
hours (199 min average across the 10). All 10 completed cleanly: 92.3-99.9% convergence, zero
`vpc.char`/crash errors. Total wall-clock: ~33 hours sequential. (The original "~12.5 hours per
envPC at 8 cores, 5+ weeks for the full 72-run version" estimate was based on including the
1000-iteration permulation loop, which -- per the emp_p fix commit -- turned out to be both far
slower and far less informative than assumed on taxon-reduced subsets; dropping it changed the
feasibility picture entirely.)

**Results — rank correlation and candidate recall vs. the full-707-taxon baseline** (envPC_1):

| Subtribe dropped | Taxa removed | ρ(-log10 p) | full p<0.001 → drop p<0.001 | full p<0.001 → drop p<0.05 |
|---|---|---|---|---|
| Andropogoninae | 150 | 0.782 | 62.6% | 97.8% |
| Anthistiriinae | 79 | 0.781 | 60.2% | 92.7% |
| Saccharinae | 33 | 0.849 | 58.5% | 97.8% |
| Ischaeminae | 23 | 0.933 | 74.5% | 100.0% |
| Triticinae | 28 | 0.952 | 75.6% | 98.5% |
| Aveninae | 24 | 0.963 | 78.2% | 99.6% |
| Tripsacinae | 33 | 0.970 | 85.5% | 100.0% |
| Apludinae | 22 | 0.972 | 85.0% | 100.0% |
| Sorghinae | 23 | 0.976 | 84.7% | 100.0% |
| Oryzinae | 35 | 0.993 | 88.0% | 100.0% |

Full-707-model output: `output/finalModels_20251002/envPC_1/ASREML_res.txt`. Per-subtribe outputs:
`output/revisionAnalyses/perSubsetModels/drop_<Subtribe>/envPC_1/`. Summary table:
`output/revisionAnalyses/jackknife_n20_envPC1_summary.csv`.

**Interpretation**: the model is broadly robust to any single subtribe's removal — 7 of 10 drops
give ρ > 0.93 and every drop recovers ≥92.7% of the original p<0.001 hits at a lenient p<0.05
threshold. The two largest subtribes (Andropogoninae, Anthistiriinae) are, unsurprisingly, the most
influential (ρ≈0.78) — but even they retain strong rank correlation and near-complete lenient
recall. One notable exception to a pure taxon-count story: Saccharinae (33 taxa) is nearly as
influential as the two much larger subtribes (ρ=0.849), while Tripsacinae -- the exact same size
(33 taxa) -- is one of the *most* stable drops (ρ=0.970). This asymmetry suggests Saccharinae
specifically (not just "a 33-taxon subtribe") carries disproportionate weight for this signal, worth
a mention alongside the size-driven pattern rather than folding into it.

### 5b. How many transition events does each jackknife actually remove?

`12J_jackknifeTransitionLoss.R`. The full-707 envPC1 ASR is computed **once** and its 46 transition
clades (27 toH + 19 toL, spanning 149 transition-carrying tips) are then held **fixed**; each
jackknife is scored purely structurally — what fraction of each fixed clade's tips the drop removes.

**Why the ASR is not re-run per jackknife.** Re-running it conflates genuine loss of
transition-bearing lineages with re-estimation effects: dropping taxa shifts the
`quantile(trait, .3)`/`quantile(trait, .7)` cutoffs and the `fastAnc` CIs, so edges that previously
missed the threshold get newly called. Scored that way, three subtribes appear to **gain**
transitions (Oryzinae +3, Triticinae +2, Aveninae +2) despite removing *zero* transition-carrying
tips — a pure artifact. Holding the ASR fixed removes that confound, and since nothing is
re-detected there is no "retained vs lost" judgement to make: a clade is wiped, partially eroded,
or untouched.

| Subtribe dropped | taxa | clades wiped (toH/toL) | partial | untouched | tips lost | % of transition tips | ρ |
|---|---|---|---|---|---|---|---|
| Andropogoninae | 150 | **14** (8/6) | 1 | 31 | 51 | **34.2%** | 0.782 |
| Anthistiriinae | 79 | **11** (8/3) | 1 | 34 | 35 | **23.5%** | 0.781 |
| Sorghinae | 23 | 2 (2/0) | 0 | 44 | 6 | 4.0% | 0.976 |
| Saccharinae | 33 | 1 (0/1) | 0 | 45 | 9 | 6.0% | 0.849 |
| Apludinae | 22 | 1 (1/0) | 0 | 45 | 7 | 4.7% | 0.972 |
| Ischaeminae | 23 | 1 (1/0) | 0 | 45 | 3 | 2.0% | 0.933 |
| Tripsacinae | 33 | 1 (1/0) | 0 | 45 | 2 | 1.3% | 0.970 |
| Aveninae | 24 | 0 | 0 | 46 | 0 | 0.0% | 0.963 |
| Oryzinae | 35 | 0 | 0 | 46 | 0 | 0.0% | 0.993 |
| Triticinae | 28 | 0 | 0 | 46 | 0 | 0.0% | 0.952 |

Correlation with the jackknife rank correlation ρ: `tips_lost` **Pearson −0.889** (Spearman −0.607),
`clades_wiped` −0.867 (−0.488), versus plain `taxa_dropped` −0.800 (−0.476). Output:
`output/revisionAnalyses/jackknife_transitionLoss_fixedASR_envPC1.csv`.

**Transition clades are almost perfectly nested within subtribes.** Across all 460 clade × subtribe
combinations there is exactly **one** partially-eroded clade (Anthistiriinae); every other clade is
either wiped outright or completely untouched. Both transition clades and subtribes are
monophyletic, so a clade is nearly always wholly inside or wholly outside a dropped subtribe. This
is why `clades_wiped` and `tips_lost` carry almost the same information here, and why "wiped" is the
natural unit.

**Interpretation**: transition loss predicts destabilization better than taxon loss, and the two
biggest drops are cleanly explained — Andropogoninae and Anthistiriinae remove 34% and 24% of all
transition-carrying tips (14 and 11 of the 46 clades outright), and have the two lowest ρ. This is
the same mechanism analysis #2 identified independently from the power simulation: what matters is
the number of independent transition events, not the number of tips. Two different analyses
converging on that conclusion is worth stating explicitly in the manuscript.

It also reframes, without fully resolving, the **Saccharinae** anomaly from §5. Saccharinae removes
6.0% of transition-carrying tips — the third-most of any subtribe — so its third-lowest ρ is at
least *rank-consistent* with the predictor, rather than the flat contradiction it appeared to be on
taxon count alone. But the slope is still steep: Apludinae removes 4.7% and sits at ρ=0.972 versus
Saccharinae's 0.849, and Saccharinae has the lowest strict recall of all ten drops (58.5%) despite
wiping only a single toL clade. Its residual from a transition-loss-only fit is −0.078, more than
double the next largest — so a Saccharinae-specific effect is real. The `12J` spread metrics below
give the most plausible account of it.

### 5b-ii. Environmental spread of each dropped subtribe

Second predictor family, computed in the same `12J` loop: how much envPC1 variance each subtribe
spans, on the model's own 707-taxon trait basis (`output/envData_707Poaceae_20250804.txt`).

| Subtribe dropped | taxa | envPC1 SD | envPC1 range | eccentricity | % total SS | frac < q30 | frac > q70 | ρ |
|---|---|---|---|---|---|---|---|---|
| Anthistiriinae | 79 | 0.76 | 2.92 | 0.03 | 6.4 | 0.25 | 0.24 | 0.781 |
| Andropogoninae | 150 | 0.69 | 3.55 | 0.42 | 13.6 | 0.11 | 0.40 | 0.782 |
| **Saccharinae** | 33 | **0.96** | **3.54** | 0.41 | 4.9 | 0.42 | 0.12 | 0.849 |
| Ischaeminae | 23 | 0.70 | 2.75 | 0.80 | 3.6 | 0.04 | 0.52 | 0.933 |
| Triticinae | 28 | 0.67 | 3.03 | 1.30 | 8.4 | 0.82 | 0.00 | 0.952 |
| Aveninae | 24 | 0.58 | 2.11 | 0.92 | 4.0 | 0.62 | 0.00 | 0.963 |
| **Tripsacinae** | 33 | 0.59 | 1.84 | 0.23 | 1.8 | 0.18 | 0.27 | 0.970 |
| Apludinae | 22 | 0.70 | 2.98 | 0.18 | 1.6 | 0.14 | 0.23 | 0.972 |
| Sorghinae | 23 | 0.80 | 2.21 | 0.19 | 2.1 | 0.30 | 0.35 | 0.976 |
| Oryzinae | 35 | 0.37 | 1.20 | 1.06 | 6.2 | 0.00 | 0.86 | 0.993 |

Correlations with ρ, all predictors together:

| Predictor | Pearson | Spearman | R² (univariate) |
|---|---|---|---|
| % transition tips lost | **−0.889** | −0.607 | **0.790** |
| clades wiped | −0.867 | −0.488 | — |
| taxa dropped | −0.800 | −0.476 | 0.641 |
| envPC1 range | −0.698 | −0.673 | 0.488 |
| % of total envPC1 SS | −0.641 | −0.564 | — |
| envPC1 SD (within) | −0.488 | −0.395 | — |
| eccentricity (\|mean − grand mean\|) | **+0.395** | +0.152 | — |

**The matched pair.** Saccharinae and Tripsacinae have identical taxon count (33) *and* identical
clades wiped (1), yet ρ 0.849 vs 0.970 and strict recall 58.5% vs 85.5%. What separates them is
environmental spread: Saccharinae has the **highest within-subtribe envPC1 SD of all ten** (0.96 vs
0.59) and nearly **twice** the range (3.54 vs 1.84) — its range essentially matches Andropogoninae's
(3.55) on a fifth of the taxa. That is the cleanest available explanation for the §5 anomaly:
Saccharinae is the most environmentally heterogeneous clade in the set, so removing it costs more
contrast than its size or transition count implies.

**Spread does not replace transition loss as the global predictor**, though. Univariate R² is 0.488
for range vs 0.790 for transition loss; adding range to transition loss lifts R² to 0.851 but at
p=0.133 with n=10, which is suggestive rather than established. Report it as the explanation for the
Saccharinae residual specifically, not as a second confirmed main effect.

**Eccentricity runs the opposite way** (r = **+0.395**): subtribes whose mean envPC1 sits far from
the grand mean are *more* stable to drop, not less. The threshold fractions show why — Oryzinae is
0% below q30 / 86% above q70 and Triticinae is 82% / 0%, so both are environmentally extreme but
sit wholly on one side of the transition thresholds and therefore generate almost no within-clade
transitions. Environmental **extremeness** and environmental **heterogeneity** are distinct
properties with opposite consequences here, and the manuscript should keep them separate rather
than collapsing both into "environmentally unusual".

---

## 6. ESM2 zero-shot mutagenesis validation (arabidopsis + maize)

Integrated from the standalone `protein_LLMscoring_pkg` (all 46 files verified against its
shipped `MD5SUMS` before copying). Validates the ESM2 zero-shot score — the `scaledESM`
predictor used by the stage-08 per-OG model — by asking whether it separates deleterious from
neutral missense mutants relative to their wild-type protein. The revision extends the original
arabidopsis benchmark to a small, **experimentally validated maize** set, showing the score works
equally well on a second species.

| dataset | genes | deleterious | neutral | PR AUC | ROC AUC |
|---|---|---|---|---|---|
| arabidopsis (TAIR10 subset) | 996 | 2,906 | 1,574 | 0.923 | 0.862 |
| maize (B73v5, experimentally validated) | 35 | 54 | 53 | 0.956 | 0.965 |

Overlay figure: `output/revisionAnalyses/mutagenesisValidation/figure/At_vs_maize_PR_ROC.png`.

**Scripts** (drivers continue the 12A-12F series; reusable workers live in `src/`):

| Script | Role |
|---|---|
| `notebook/12_revisionAnalyses/12G_buildDeleteriousSet.py` | maize: WT + deleterious FASTA from the validated mutation table, with a per-mutation validation report |
| `notebook/12_revisionAnalyses/12H_fetchHapmapRegions.sh` | maize: download HapMap3.2.1 282-panel VCFs, extract the 35 gene regions |
| `notebook/12_revisionAnalyses/12I_buildNeutralSet.py` | maize: pick common missense SNPs as the neutral set, assemble the final FASTA |
| `src/4_ESM_logits.py` | step 1: per-residue ESM2 probabilities → `.npz` (**already in this repo, byte-identical to the package copy — not duplicated**) |
| `src/10_logit2zeroShot_MutagenesisValidation.R` | step 2: `.npz` → per-sequence zero-shot score (see the divergence note below) |
| `src/mutagenesisValidation_visualization.R` | step 3: gene-level contrast, PR/ROC, figures, tables |
| `src/mutagenesisValidation_compareDatasets.R` | overlay PR/ROC of several datasets on one figure |

**Divergence from the pipeline's own `src/10_logit2zeroShot.R` — read before reusing either.**
The benchmark copy is a hardened rewrite of the same script and is *not* a drop-in replacement.
Both are kept deliberately: `src/10_logit2zeroShot.R` is the provenance record for the published
stage-07 → stage-08 results and stays frozen; the `_MutagenesisValidation` copy is what this
benchmark runs. Differences:

| | `10_logit2zeroShot.R` (pipeline) | `10_logit2zeroShot_MutagenesisValidation.R` |
|---|---|---|
| logit ↔ sequence matching | positional (`testfasta[[j]]`) | **by name** (`seqByID[[seqID]]`) |
| gene/group key for z-scaling | hardcoded split on `":"` | `groupRegex` 4th arg, default `"_\| "` |
| output `row.names` | `TRUE` | `FALSE` |
| `RETICULATE_PYTHON` | hardcoded `/usr/bin/python` | overridable via env var |

The grouping difference is the load-bearing one: stage-07 OG seqIDs look like
`OG0006571:10wheat_assembly_jagger:0`, so the pipeline's `":"` split yields the orthogroup, while
the benchmark default would yield `OG0006571:10wheat`. Swapping the files without passing
`groupRegex=":"` would silently corrupt the `scaledESM` predictor. **Open item:** the
positional-vs-name matching difference is a latent correctness question for the pipeline copy, not
just a style change — if `4_ESM_logits.py`'s npz key order ever differs from the input FASTA record
order, the pipeline version scores sequences against the wrong protein. Worth checking against real
stage-07 inputs; tracked separately from this integration.

**Environment**: `envs/environment-mutagenesisValidation-full.yml` (full conda export of the tested
`pytorch-2-1` env), with minimal pins in `envs/requirements-mutagenesisValidation-{python,R}.txt`
(torch 2.1.1, transformers 4.40.0, biopython 1.83, numpy 1.26.4; R needs reticulate, Biostrings,
limma, PRROC). Kept as a separate env rather than merged into `environment-py.yml` — these pins are
specific to the benchmark and merging risks disturbing a working solve. GPU helps only step 1.

**What is and isn't tracked in git** (see `.gitignore`'s 12G-12I block and `DATA.md`):
the ESM2 logit `.npz` files (259 MB arabidopsis, 6.3 MB maize) and the large public B73v5
references (30 MB peptide FASTA, 12 MB GFF3) are **not** tracked. The zero-shot score tables **are**,
so steps 3+ (figures, PR/ROC) reproduce from a clean clone without a GPU or the npz. The HapMap
region extract is tracked because regenerating it means downloading 5.75 GB of VCFs.

**Known caveats** (carried over from the source package, unchanged): 6 of the 35 maize genes have
no eligible neutral site (`Zm00001eb153340, 190120, 194210, 215710, 234730, 373430`), so their
deleterious mutants enter PR/ROC without a within-gene neutral comparator; `Zm00001eb336930` drew
two SNPs in the same codon (N77S, N77D); and with 107 mutants total the maize AUCs have wide
confidence intervals. `Zm00001eb252380_P001 A139T` is applied as **A140T** (an off-by-one in the
source table) and flagged `OK_CORRECTED_FROM_A139T` in `maize_deleterious_validation.txt`.

---

## Open items before implementation starts

- Confirm the FPKM > 1 detectability threshold for tissue-expression breadth (#3), or supply a
  preferred alternative.
- Decide whether `AN21TNTL0178` should be added to the taxonomy CSV (with its correct
  tribe/subtribe) before these analyses run, or left excluded as the default above.
