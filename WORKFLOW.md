# WORKFLOW.md — Analysis Pipeline, Stage by Stage

The pipeline runs as a numbered sequence of stages (`01` → `11`). Each row below maps in
the data it needs (cross-referenced to `DATA.md`), the scripts/notebooks that do the work,
and what it hands off to the next stage.

**Current data readiness:** every **T1** item in `DATA.md` is now physically present in
this repo (`data/`, `output/`, or `src/`) or reachable via the noted external source
(HuggingFace, GBIF/BIEN live query, etc.) — see `DATA.md` for exact paths. **T2** items are
already provided as ready-made intermediates, so you don't have to regenerate them unless
you want to re-derive them from scratch.

---

## Pipeline order

| Stage | Purpose | Data in (`DATA.md`) | Key scripts | Output → feeds |
|---|---|---|---|---|
| **01** `genomeAssembly` | Short-read genome assembly (megahit) | Raw sequencing reads (external; not tracked in this repo) | `notebook/01_genomeAssembly/README.md` — external pipeline: [bucklerlab/p_reelgene](https://bitbucket.org/bucklerlab/p_reelgene/src/master/short_read_assembly/) | Raw assembly FASTAs → `data/assemblies/` → 02 |
| **02** `metadataCuration` | Merge PanAnd/LIMS QC metadata + manual species-ID curation into a filtered genome list | Poaceae accession metadata; per-assembly QC statistics | `02A_metadataProcessing.ipynb`, `02B_furtherFilter.ipynb` | Filtered metadata table → 03, 05, 08 |
| **03** `orthogroup` | OrthoFinder (32 representative genomes) → OG filtering → ancestral AA sequence reconstruction (for the filtered OGs) → miniprot cross-mapping to Zm/Pv/At/Angiosperms353 → TableS5 summary; OG-name translation | 32 representative genome assemblies + Helixer annotations; rice→OG mapping; Angiosperms353 Oryza reference; maize v5 mRNA reference; miniprot GFF annotations | `03A_buildHelixerOG.sh` (build+OrthoFinder), `03B_OGFilter.ipynb` (filter), `03C_ancestralSeqReconstruction.sh` (ancestral seq + miniprot ID matching), `03D_miniProtResult_eval.ipynb` (miniprot eval), `03E_TableS5Generation.ipynb` (TableS5), `03F_OGtranslation.R` (OG-name translation) | Filtered OG list + ancestral sequences + OG→maize/Pv/At/Angiosperms353 mapping → 04, 07, 08, 09, 11 |
| **04** `msaGeneration` | Per-OG CDS multiple sequence alignment (mafft) | Filtered OG list + CDS sequences (from 03) | `notebook/04_msaGeneration/README.md` — `mafft --ep 0 --genafpair --maxiterate 1000 <input> > <output>` | Per-OG MSAs (`output/OrthofinderMAFFT/*_mafft.fa` — **no longer retained on disk**, consumed by 05A's gap-stripping step; regenerate by rerunning this stage) → 05 (gap-stripping/gene trees), 07 (dN/dS calculation needs the MSA directly), 09 (HyPhy RELAX needs the MSA directly) |
| **05** `phylotreeConstruction` | Gap-strip CDS MSAs → extract angiosperm353 per-gene sequences + genetic distance → RAxML gene trees → ASTRAL-Pro species tree → filter/visualize/annotate species tree → phylogenetic K (relatedness) matrix | Gap-stripped CDS MSAs (from 04); angiosperm353 OG-name list; Poaceae metadata (for tree filtering) | `05A_treeConstruction` (gap-strip, RAxML, ASTRAL), `src/S04_angiosperm353_extractAndDist.R` (invoked from 05A), `05B_neutralPhylogenyVisualization.ipynb` (filter/visualize/phyloK) | Species tree + phyloK matrix → 08 (predictor); per-OG gene trees → 09 (HyPhy RELAX runs on the gene tree from 05) |
| **06** `envirotyping` | Species occurrence coordinates → WorldClim/soil rasters → habitat summary → envPC1–3 (PCA) → visualize distributions/tree overlay → ancestral state reconstruction | Species-name list; GBIF/BIEN occurrence records; WorldClim + soil rasters; environmental metadata; derived occurrence dataset (Zenodo) | `06A_spCoordEnvData.sh` (coords → env data → envPC), `src/08_pulling_geo_data.R`/`src/09_pulling_envData.r` (invoked from 06A), `src/S05_envPC_analysis.R` (envPC computation, invoked from 06A), `06B_visualizationEnvAdapt.ipynb` (visualize + ASR) | envPC1–3 table (Fig. 1) → 08; ASR transition nodes → power simulation (08E) |
| **07** `summaryStats` | Per-OG premature-stop/frameshift calling, tip-to-outgroup dN/dS calculation, ESM2 & PlantCAD zero-shot scores | miniprot GFF annotations; seqIDmapping tables; OrthoFinder protein MSAs; gap-stripped CDS MSAs (from 05A); ESM2 weights; PlantCAD weights | `07Aa` (frameshift), `07Ba` (premature stop, local), `07Bb` (dN/dS, local); `07Ca` (ESM2, SCINET GPU); `07Da`/`07Db`/`07Dc` (PlantCAD, SCINET GPU) | Per-OG activity scores + dN/dS table (Fig. 4) → 08 |
| **08** `linearModeling` | Master data table → proteome a.a./genome GC estimation → genome-wide feature association (Fig. 3) → per-OG phylogenetic mixed model + permulation (Fig. 5) → power simulation (Fig. 2) | dN/dS table (from 07, used as a predictor); OG→maize mapping; maize v5 expression (FPKM) | `08A_masterDataTableGeneration.ipynb`, `08B_genomicFeatureEstimation.ipynb`, `08C_genomicFeatureAssociation.ipynb`, `08D_perOGmodel.sh`, `08E_power_simulation.sh` | Candidate-OG lists + model results → 09, 11 |
| **09** `molEvolution` | MSA cleaning (from 04) → RAxML gene trees (from 05) → foreground/background branch labeling → HyPhy RELAX selection-intensity tests per trait | OG→maize mapping | `09A_HyPhyPipeline.sh`, `09B_RELAX_resultSummary.ipynb` | RELAX result tables → 11 |
| **10** `aprioriCandidate` | OG→gene-ID mapping (Helixer) via miniprot → convert manually-extracted DE gene lists per abiotic-stress study to OG IDs → consistent (multi-study) per-stress candidate gene sets | Per-species CDS FASTAs (stress genes) + DEG study metadata; rice→OG mapping | `10A_DEG_IDconversion.sh` (ID mapping), `10B_consistentStressResponsiveGenes.ipynb` (DE-gene ID conversion + consistency filtering) | Gene-ID mapping + consistent per-stress DEG sets → 11 |
| **11** `candidateOGInvestigation` | Integrate ASReml (08) + RELAX (09) + gene-ID mapping (10) + expression/GO evidence → final candidate OG list (Fig. 5) | DeepGO GO annotation; Maize v5 GO annotation | `11A_candidateOGInvestigation.ipynb` (core: p-value distributions, DE/GO enrichment, 3-layer overlap), `11B_candidateGeneResidueModeling.ipynb` (per-residue modeling for 2 validated candidates) | 27 high-confidence candidate OGs (final — see "Note on 11", corrected from a previously-reported 17) |

`slurm/`, `XX_archived/`, and `*/archived/` subfolders hold SLURM job templates and
superseded/exploratory notebooks — not part of the active sequence above.

**Note on 10/11:** `10A_DEG_IDconversion.sh` and the notebook formerly named
`10B_candidateOGInvestigation.ipynb` used to be presented as one stage. They're distinct
enough in purpose (ID mapping vs. final candidate-list integration) to warrant separate
numbers and separate folders — the notebook has been moved and renamed to
`notebook/11_candidateOGInvestigation/11_candidateOGInvestigation.ipynb` (its in-progress
edits carried over with the move, untouched otherwise). **Update:** since renamed again, to
`11A_candidateOGInvestigation.ipynb`, to make room alongside `11B_candidateGeneResidueModeling.ipynb`.

**Note on notebook archiving:** `10B_stressInducedGene_enrichment.ipynb` had been
mistakenly filed under `notebook/10_aprioriCandidate/archived/` despite being active
(not superseded) — moved back up to `notebook/10_aprioriCandidate/`.

**Note on 03A–03F:** originally `03A_buildHelixerOG.sh` bundled OG construction with
ancestral-sequence reconstruction, and `03B_OGFilter.ipynb` bundled OG filtering with
TableS5 generation and a standalone gffcompare-based ID-matching validation check. Per
author review, the real dependency order is OG construction → filter → ancestral
reconstruction (for the filtered OGs only) → TableS5, so these have been split into their
own files (03A/03C from the old 03A; 03B/03D from the old 03B) and the old
`03C_miniProtResult_eval.ipynb`/`03D_OGtranslation.R` renumbered to 03E/03F accordingly.
The ID-matching validation check (output read by nothing downstream) was moved to
`notebook/03_orthogroup/archived/03B_gffCompIDMatching.ipynb` rather than renumbered, since
it isn't part of the active pipeline. **Update:** 03D and 03E were later swapped from the
above (miniProtResult_eval is now 03D, TableS5Generation is now 03E) so that evaluating the
miniprot cross-mapping (03D) directly follows where it's generated (03C), with the TableS5
summary wrap-up (03E) last.

**Note on miniprot cross-mapping (03C):** the miniprot ID-matching block (OG→maize/Pv/At/
Angiosperms353 mapping) originally sat at the end of `05A_treeConstruction`, unrelated to
that stage's actual tree-building work. Per author review, it's been moved into 03C
(right after the ancestral-sequence fasta it depends on is generated) along with dropping
several superseded pre-"_v2" mapping lines that referenced an older, no-longer-produced
ancestral-seq file. `data/OGToZm_mapping_v2.txt`'s write target was also corrected from
`output/` to `data/` in the move, matching where every active consumer (03E, 08B, 11)
actually reads it from and where the current (2024-09-16) file lives. (The archived
`09B_GWASRelaxEnrichment.ipynb` also reads it, but that's a superseded side-analysis, not
the active `09B_RELAX_resultSummary.ipynb`; and `src/12_runPermulation_perOGModel.R`'s own
read of this file is entirely commented-out dead code, not an active consumer — both
excluded from the "Used by" list above and in `DATA.md`.)

**Note on 05A/05B:** `05B_neutralPhylogenyVisualization.ipynb` originally mixed three
unrelated things: (1) extraction of angiosperm353 per-gene sequences from the gap-stripped
per-OG MSAs plus a K81 genetic-distance calc — this was actually prep that 05A's RAxML step
needed but never had a producer for (`output/geneTree_angiosperm353/*.fa` had no source
anywhere in the repo); (2) the core ASTRAL species-tree filter/visualize/phyloK block; and
(3) supplemental cross-project analyses (a 14-species divergence-time figure, a Zea-only
subtree, and a full comparison against an external collaborator's independently-built tree).
Per author review: (1) moved to new `src/S04_angiosperm353_extractAndDist.R`, invoked from
05A right before its RAxML step; (2) stays in a trimmed `05B_neutralPhylogenyVisualization.ipynb`;
(3) moved to `notebook/05_phylotreeConstruction/archived/05B_supplementalTreeComparisons.ipynb`
(unparameterized, matching the archived-code convention). Two real bugs fixed in the process:
`05B` read `output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk` as its species-tree
input, but 05A's active `astral-pro` call was writing to `output/PoaceaeTree_angiosperm353.nwk`
— a filename with no real historical output and no downstream reader anywhere in the repo, while
the actual, always-real output every consumer reads is `_astral.nwk`-suffixed. Corrected 05A's
`astral-pro` target to `output/PoaceaeTree_angiosperm353_astral.nwk`, and verified via a real
rerun of the now-fixed 05B (800 tips reproduced) that this is genuinely what 05B expects to read
downstream of that step. Also added the missing `mkdir -p` calls 05A needed before `cd`-ing into
`output/geneTree_allOGs/`/`output/geneTree_angiosperm353/` (neither directory was ever created
by the script itself — a real crash-on-fresh-run bug, not just a naming issue). **Separately
flagged, not fixed in this pass:** `output/geneTree_allOGs/` has no producer anywhere in this
repo for its actual per-OG gene-tree `.fa` inputs (unlike `output/geneTree_angiosperm353/`, whose
`.fa` inputs are written by the new `src/S04_angiosperm353_extractAndDist.R`) — the `mkdir` fix
only prevents an outright crash on a fresh run, it doesn't close this gap. Needs its own
investigation before 05A's "all-OGs" gene-tree branch can be considered reproducible end-to-end.
05B's filtered/labeled tree outputs were writing "astral3"-suffixed filenames matching the *inactive*
`astral-pro3` line in 05A (commented out) rather than the non-"3", `_20250407`-dated naming
every real downstream consumer (`08C`/`src/12_runPermulation_perOGModel.R`, `06B`, `08B`, `11`)
actually reads — corrected to `output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk`
and `..._astral_spLabeled_20250407.nwk`/`angiosperm353_astral_spLabeled_20250407.png` (dated,
matching the same `20250407` batch as `phyloK_728Poaceae_astral_20250407.txt`).

**Update:** two more real bugs found and fixed. (1) `src/S01_phyloK.R`'s `phyloK()` had a
calculation error — `out = 2*out/max_brlen` doubled the shared-branch-length ratio (max value
across the matrix was 2, when a proper ratio should max out at 1); corrected to
`out = out/max_brlen`. This is exactly what the `phyloKMat = phyloKMat/2 # to correct the error`
line downstream in `08B`/`08C` had been silently compensating for — now that the function
itself is fixed, that line has been removed from both (see stage 08's note). (2) `05B`'s
metadata read for tip-set filtering pointed at `data/Poaceae_metadata_highErrorFiltered_2025.10.08.tsv`
— a **stage-08 output** (08A's PMS/frameshift QC filter), created after stage 05 in pipeline
order. Rerunning `05B` as it stood would have silently pruned to a different, smaller tip set
(687 instead of 727) than every historical downstream file — corrected to
`data/Poaceae_metadata_filtered_2025.08.28.tsv`, verified to reproduce the exact historical
727-tip set. `output/phyloK_728Poaceae_astral_20250407.txt` was regenerated with both fixes
(correct tree/metadata reproducing the same 727 tips, corrected `phyloK()` formula) and the
real file replaced — confirmed the new matrix is exactly half the old one, element-wise, to
floating-point precision (max abs diff between old and `2 × new` is `4e-15` across all
727×727 entries).

**Note on 06A/06B:** `06B_visualizationEnvAdapt.ipynb` originally mixed four things: (1) envPC
computation (PCA over per-species environmental-feature quantiles, writing
`envData_707Poaceae_*`/HyPhy adaptation-list files); (2) the core KG3-climate/tree-overlay/
ancestral-state-reconstruction visualization; (3) a large "life history paper" side-analysis
(annual/perennial transitions, rhizome association) reusing the same envPC/tree objects; and
(4) a few stray cells (OG dN/dS counts) unrelated to envirotyping entirely. Per author review:
(1) moved to new `src/S05_envPC_analysis.R`, invoked from `06A_spCoordEnvData.sh` right after
the (newly wired-up) `src/09_pulling_envData.r` call; (2) stays in a trimmed
`06B_visualizationEnvAdapt.ipynb`; (3) moved to
`archived/06B_supplementalLifeHistory.ipynb` (unparameterized, matching the archived-code
convention — this archived code depends on objects computed in the original monolithic
notebook and isn't runnable standalone); (4) deleted (not part of any figure). **Update:**
`06B_supplementalLifeHistory.ipynb` has since moved again, out of this repo entirely, to
`p_lifeHistoryPoaceae/notebook/06_envirotyping/supplementalLifeHistory.ipynb` — see "Note on
the life-history/rhizome thread's extraction" near the end of this file. `06D_supplFig_envPCpipeline.R` — split out on the assumption that `library(raster)`/
`library(terra)` would namespace-mask `dplyr`/`ggplot2` functions used elsewhere in the
notebook — is merged back into `06B` and the file removed; tested empirically (real-data
sandbox run) and no such conflict actually occurs with the packages currently used elsewhere
in the notebook, so both are attached normally (fully-qualifying `raster::`/`terra::` calls
alone isn't sufficient — `subset()`/`points()` need the packages attached to dispatch their
S3/S4 methods on raster/terra objects). `src/09_pulling_envData.r`
is now `commandArgs`-wrapped like `src/08_pulling_geo_data.R` (input coords file, output env-data
path), and a new `GIS_DATA_ROOT` env var (default `/workdir/sh2246/p_evolBNI/data/GIS_env_data`,
mirroring `PHYLOGWAS_ROOT`'s pattern) parameterizes its 3 GIS raster sources that live outside
this repo entirely, in the sibling `p_evolBNI` project. Two real bugs fixed in the process:
`09_pulling_envData.r`'s GSDE variable-name extraction assumed a hardcoded path depth
(`strsplit2(...)[,8]`) that only happened to work for the original `p_evolBNI` path — fixed to
`basename(...)`; and `06D`'s occurrence-species filter used a positional column index
(`bien_data_clean[,3]`) that pointed at `scientificName` in its old frozen input file but at
`decimalLatitude` (silently matching nothing) in the live `coordinates_clean.csv` it now reads
— fixed to reference the column by name. A new small intermediate,
`output/KG3_perSpecies_20250804.txt`, was introduced so `06B`'s KG3 analysis (which needs the
per-species dominant Köppen class, only derivable from the raw per-occurrence env data) doesn't
need to re-read the 194MB raw env-data file itself — `S05` derives and persists it once,
alongside `envData_707Poaceae_*`.

**Update:** confirmed `output/KG3_perSpecies_20250804.txt` and
`output/PoaceaeTree_angiosperm353_astral_filtered_withEnvData_20250804.nwk` were genuinely
missing on disk (only older `_20250414`/`_20250416`-dated versions of the latter existed,
i.e. `S05_envPC_analysis.R` had only been partially rerun since the stage-06 split) — reran
the script for real (all T1 inputs already present; ~30s) and confirmed both now exist and
`06B` can read them. **A real non-determinism bug surfaced during verification:**
`process_synthetic()` (`src/process_synthetic_fun.R`) calls `missMDA::estim_ncpPCA`/
`imputePCA` to impute missing environmental-feature values before PCA — both are stochastic
(k-fold CV / random imputation) and the script never called `set.seed()`, unlike every other
stochastic step in this repo (`08C`, `src/05_4dMSA_sampling.R`, etc., all seeded with `123`).
Rerunning produced a `envData_707Poaceae_20250804.txt` that matched the previous version's
707-species membership exactly but differed by small floating-point deltas in the 16 species
that needed imputation, which then propagated into the cold/warm-adapted-assemblies
percentile-cutoff lists (1-2 species shifted across the 30th/70th-percentile threshold).
Fixed by adding `set.seed(123)` immediately before the `process_synthetic()` call — verified
across 3 consecutive reruns that every output (`envData_707Poaceae_*`, KG3, the withEnvData
tree, and all 6 percentile-based adaptation lists) is now byte-identical run to run.

**Note on 07:** stage 07 computes 4 independent per-OG/per-sequence stats, each merged into one
master table in `08A_masterDataTableGeneration.ipynb`. Per author review:
- **Premature stop** (`07Ba_PMS_run.sh`, replacing the SCINET-only `07Ba_PMS_SCINET.sh`, now
  archived): `src/11_find_premature_stops.py` was already correct — this is now a plain local
  GNU-parallel runner over `output/orthofinderProteinMSAs_fullset_20250710/` (the current, full
  OG set, superseding the old main/`_additionalOGs` SCINET-batch split), writing
  `output/combined_PMS_20250421.txt` (kept as the exact filename `08A` reads, even though the
  file itself is regenerable — matches this project's convention of keeping historically-read
  filenames stable). `07Bb_PMS_run.sh` was a byte-identical duplicate of
  `notebook/slurm/PMS_run.sh` and was removed rather than archived.
- **Frameshift** (`07Aa_frameShiftMutation.ipynb`): reads `output/seqIDmapping.txt` /
  `_additionalOGs.txt` (migrated into `output/` from `p_phyloGWAS_archived`, never brought in
  during the earlier data-consolidation pass) against `Frameshift`-flagged miniprot GFF
  entries. MPID (miniprot's own per-assembly protein ID) is **not** globally unique — it's
  assigned independently within each of the two separate miniprot batches, so the same
  `MPID:assembly` key can mean different genes in the main vs. additionalOGs batch. Per author
  review, the two batches are kept as two fully separate passes (their own GFF dir +
  seqIDmapping-file pairing each), each writing its own output
  (`output/frameShiftMutation.txt` / `_additionalOGs.txt`) — not merged — matching what `08A`
  already does itself (`FS = rbind(fread(frameShiftMutation.txt), fread(frameShiftMutation_additionalOGs.txt))`).
  Verified: the rewritten notebook reproduces both real historical output files byte-for-byte.
- **dN/dS to reference** (new `07Bb_dNdS_run.sh`): `src/07B_getOmega2Ref.R` was already
  correct/portable (CLI args) but had never been wired to a driver anywhere in the repo — this
  new script loops it (via GNU parallel) over `output/CDSMSAPerOG_gs/*.gs.fa` (05A's
  gap-stripped MSAs; the raw MAFFT output this step originally read no longer exists on disk,
  consumed by 05A's own gap-stripping and not retained, but `07B_getOmega2Ref.R`'s
  reference-based gap-stripping is a no-op on already-stripped input) with
  `--ref ASM1935983v1` (the same outgroup rooting the species tree in 05/06), reproducing
  `data/fullSetOGs_240903.txt` — verified byte-for-byte identical against the existing file
  for two real OGs.
- **LLM zero-shot scores**: `src/10_logit2zeroShot.R` (not `src/4_ESM_logits_to_zero_shot.py`,
  archived — it has a real bug where every sequence's result overwrites the same dict key,
  silently keeping only the last one) is the authoritative ESM2 zero-shot conversion script,
  matching what `07Ca_ESM_SCINET.sh` already calls. `src/5_PlantCAD_logits.py` (new) adapts
  `src/4_ESM_logits.py`'s structure for PlantCaduceus (`kuleshov-group/PlantCaduceus_l32`), a
  DNA-sequence Caduceus/Mamba model, not a standard transformer: lowercase 4-letter nucleotide
  vocab instead of the 20-amino-acid one, no leading special token (confirmed empirically,
  unlike ESM's `<cls>`), and `trust_remote_code=True` for both tokenizer and model. Verified as
  far as possible without a GPU (tokenizer + model loading + code path all confirmed correct up
  to the point PlantCaduceus's `mamba_ssm` backend requires actual CUDA — it has no CPU
  fallback, confirmed by a real forward-pass attempt failing inside its Triton kernel) — a full
  run needs to happen on SCINET/GPU, matching this stage's existing convention for LLM scoring.
  New `src/6_PlantCAD_logit2zeroShot.R` (DNA analog of `src/10_logit2zeroShot.R`: same
  mean-log-ratio zero-shot approach, lowercase 4-letter vocab and `readDNAStringSet` instead
  of `readAAStringSet`) is driven by a new `07Da_PlantCAD_SCINET.sh`/`07Db_PlantCAD_run.sh`/
  `07Dc_PlantCADzeroshot_run.sh` series mirroring 07C's structure, with one difference: since
  PlantCAD's input is a single combined nucleotide FASTA rather than one file per OG, these
  run as 2 single SLURM jobs (chained via `--dependency=afterok`) instead of array jobs over
  per-OG cmd files — `07Cb_ESM_run.sh`/`07Cc_zeroshot_run.sh` turned out to be byte-identical
  duplicates of the generic `notebook/slurm/ESM_run.sh`/`zeroshot_run.sh` templates, so this
  series doesn't recreate that duplication. Verified `src/6_PlantCAD_logit2zeroShot.R`'s logic
  end-to-end against synthetic logits (real PlantCAD logits need a GPU) — correct scores and
  per-OG z-scaling; the only failure hit while testing (`reticulate`/numpy 2.x incompatibility
  under this machine's old `reticulate` build) is a pre-existing fragility already present in
  the original `10_logit2zeroShot.R`'s hardcoded `RETICULATE_PYTHON` path, not something new.
- `07Ab_prematureStopCodon.ipynb` was actually a rhizome/life-history enrichment test on
  stop-codon presence, unrelated to this stage's own premature-stop scoring — archived as
  `archived/07Ab_prematureStopCodon_rhizomeLifeHistory.ipynb`. **Update:** since moved out of
  this repo entirely, to `p_lifeHistoryPoaceae/notebook/07_summaryStats/prematureStopCodon_rhizomeLifeHistory.ipynb`
  — see "Note on the life-history/rhizome thread's extraction" near the end of this file.

**Note on 08A–08E:** stage 08 was previously `08A_masterDataTableGeneration.ipynb`,
`08B_genomicFeatureAssociation.ipynb`, `08C_perOGmodel.sh`, `08D_power_simulation.sh`. Per
author review, `08B` genuinely mixed two things: proteome amino-acid composition + genome GC
**estimation and investigation** (per-taxa/per-OG a.a. composition, FPKM-abundance filtering,
`canprot` physicochemical features, per-taxa/per-OG GC content — writing the 3 real feature
tables `poaceae_aaComposition_20251008.txt`, `poaceae_dnaComposition_20251008.txt`,
`genomicFeatureData_20251008.txt` + variants), versus the phylogenetic **modeling** built on
top of those tables (variance partitioning, Fig4a/b, ASReml mixed model vs. envPCs,
permulation, empirical p-values, Fig4c). Split into `08B_genomicFeatureEstimation.ipynb` and
`08C_genomicFeatureAssociation.ipynb`, each self-contained (its own copy of the metadata/
phyloK/envPC/tree loading header) rather than sharing in-memory state — `08C` re-reads the
feature tables `08B` writes, exactly as it already did before the split. The two shell
drivers were renumbered `08C_perOGmodel.sh` → `08D_perOGmodel.sh` and
`08D_power_simulation.sh` → `08E_power_simulation.sh` to keep a single A–E letter sequence
(the latter also gained the `PHYLOGWAS_ROOT`/`cd` header every other driver in this pass has —
`src/powerSimulation_XY_revised.R` already resolved `PHYLOGWAS_ROOT` itself, just needed it
exported first). A trailing "Supplemental analysis for binomial modeling (life history)"
block, reusing the same feature tables to test lifeHistory/rhizome associations, was archived
verbatim as `archived/08X_binomialLifeHistoryModel.ipynb` (unparameterized, matching the
established archived-code convention — it also has a pre-existing bug, undefined
`correctionFactor`/`correctionFactor2` variables, left as-is since archived code isn't fixed).
**Update:** since moved out of this repo entirely, to
`p_lifeHistoryPoaceae/notebook/08_linearModeling/binomialLifeHistoryModel.ipynb` — see "Note
on the life-history/rhizome thread's extraction" near the end of this file. The archived
`08B_GOenrichmentForAnnualism.ipynb` (a separate GO-enrichment side-analysis for the same
annual/perennial contrast, not previously called out in this note) moved there too.

Three real bugs fixed in the process (all in `08A`, all pre-existing — i.e. present before this
pass, not introduced by it): (1) `08A`'s GO-enrichment block referenced `mappingFileMerged`
without ever building it (an undefined variable) — reconstructed from `data/OGToZm_mapping_v2.txt`
(the only mapping file both used downstream and actually present on disk; a second file, Pv
mapping, that 10B/11 also reference doesn't exist anywhere in this repo — flagged, not fixed,
since 10B/11 haven't been reached in this cleanup pass yet); (2) the same block then indexed the
result positionally (`mappingFileMerged[...,3]`) — fixed to reference the `ZmID` column by name,
the same class of fix already applied to 06's occurrence-filter column index; (3) a
`conservedOG`-based GO-enrichment branch (cells computing `tgd2`/`GO_res_table2`/`topGOTab2`)
referenced a `conservedOG` list that was only ever defined inside a commented-out block — dead,
unreachable code (Fig2 only ever consumes `topGOTab`, built from the real `lostOG` list), so
removed rather than fixed, along with an unrelated broken diagnostic print
(`length(geneLosscountPerTaxaa)`, a typo'd variable name) and a stale `annual_assemblies_20250423.txt`
write superseded by (and never read instead of) the newer `annual_assemblies_20260213.txt`
(see below — `08A` no longer writes this itself). A fourth bug, in what's now `08B`: the
a.a.-physiochemical-properties cells referenced `aa.comp.busco`/`aa.feat.busco`, never assigned
anywhere in the notebook (dead BUSCO-based investigation, abandoned before completion) — these
lines were removed since nothing downstream consumes them either. One inconsistency flagged when
the split first landed — `08A` read `phyloK_728Poaceae_astral_20250407.txt` as-is, while `08B`/
`08C` applied an extra `phyloKMat = phyloKMat/2 # to correct the error` right after loading the
same file — has since been resolved by the author: the `/2` correction is no longer needed now
that the phyloK-generating function itself has been fixed upstream, so it was removed from `08B`
(matching `08A`'s and `08C`'s uncorrected treatment).

**Update:** the "writing out assembly ID for annual-perennial contrast" block (writing
`data/annual_assemblies_20260213.txt`, `perennial_assemblies_20260213.txt`,
`perennialRhizome_assemblies_20260303.txt`) was moved out of `08A` entirely, into
`archived/08X_binomialLifeHistoryModel.ipynb`. These 3 files only ever fed life-history/
rhizome-specific analyses — this same archived binomial model, and the life-history/rhizome
sections of `09A_HyPhyPipeline.sh` — not 08A's core (trait-agnostic) master-table pipeline.
Per author: the life-history/rhizome thread as a whole is intended to move to its own,
separate repo once this repository's cleanup is finished. **Update:** now extracted from
`09A_HyPhyPipeline.sh` into `archived/09A_lifeHistoryRhizome_HyPhyPipeline.sh` — see the
"Note on 09A/09B" below.

**Note on 09A/09B:** `09A_HyPhyPipeline.sh` mixed the real envPC-based HyPhy RELAX pipeline
(cold/warm, drought/wet, sand/clay branch tests on candidate OGs from `08D`) with a
life-history/rhizome-specific thread (separate MSA/tree rebuilds restricted to perennial or
nonrhizomatous assemblies, target-tip extraction, tree labeling, and RELAX tests for 5
annual-perennial/rhizome model variants). Per author, the life-history/rhizome portion is
extracted verbatim into `archived/09A_lifeHistoryRhizome_HyPhyPipeline.sh` (unparameterized,
matching the archived-code convention; depends on the trimmed `09A`'s steps 1/2/5 having
already created `output/CDSMSAPerOG_HyPhy_20250203/`, `output/geneTree_allOGs_20250203/`,
and the `output/HyPhyResult/` parent directory) — same rationale as the `08A`→`08X` move:
this whole thread is slated to move to its own repo once cleanup here is finished. The
trimmed `09A_HyPhyPipeline.sh` now runs only the envPC/cold/warm/drought/wet/sand/clay
branch-test pipeline. **Update:** `archived/09A_lifeHistoryRhizome_HyPhyPipeline.sh` has since
moved out of this repo entirely, to
`p_lifeHistoryPoaceae/notebook/09_molEvolution/lifeHistoryRhizome_HyPhyPipeline.sh` — but it
still must be run from *inside a p_phyloGWAS checkout* (it still `cd`s into `PHYLOGWAS_ROOT`
and reads this repo's un-vendored upstream intermediates); see "Note on the life-history/
rhizome thread's extraction" near the end of this file.

`09B_RELAX_resultSummary.ipynb`'s core job is reading HyPhy RELAX json output and writing
per-trait result tables (`RELAX_resultTable_*.txt`). Per author review: dropped several
diagnostic-only cells with no saved output (an `install.packages("rjson")` call — package
availability belongs in environment setup, not an inline notebook call; two ad hoc QQ-style
plots; bare `OG`/`table()` prints; a single-JSON deep-dive with a histogram + t-test; a
log2(k)-statistic histogram) — none of these feed the result-table writes. A trailing
"GWAS (08D) vs. HyPhy RELAX overlap/enrichment" block (Fisher's-exact overlap tests between
envPC1 GWAS hits and RELAX candidates across top-N cutoffs, bubble plots, GWAS/RELAX
correlation plots, topGO enrichment on "Relax" vs. "Intensify" OG sets) was a distinct
side-analysis, not part of the core json→txt conversion — moved to
`archived/09B_GWASRelaxEnrichment.ipynb` (unparameterized, depends on `testRes`/`testRes2`
built earlier in `09B`, isn't runnable standalone, matching the archived-code convention).

**Update:** `09B`'s 3 remaining life-history/rhizome cells (`testRes7`/`testRes8`/`testRes9` +
their `jsonFileList`s and `fwrite`s, reading `output/HyPhyResult/{lifeHistory,rhizome_20260515,
lifeHistory_20260515}/` and writing `RELAX_resultTable_{lifeHistory_20250216,rhizome_20260524,
lifeHistory_20260524}.txt`) have been extracted out of this still-active notebook into
`p_lifeHistoryPoaceae/notebook/09_molEvolution/LH_RELAX_resultConversion.ipynb` — the one
core-embedded piece of the life-history/rhizome thread's extraction (see "Note on the
life-history/rhizome thread's extraction" near the end of this file). `09B` now only builds
the core cold/warm/envPC2/envPC3 result tables.

**Note on 10A/10B:** `10B_stressInducedGene_enrichment.ipynb`'s core job is converting
manually-extracted DE gene lists (one per abiotic-stress study, tracked in
`data/DEG_study_metadata.csv`) to OG IDs via 10A's miniprot-based mapping files, then
keeping only the genes flagged consistently (in >2 study-species pairs) per stress
category, writing `output/candidateGenes/consistentEnvResponsiveGenes.json` for stage 11.
Renamed to `10B_consistentStressResponsiveGenes.ipynb` to reflect this — the old name
described what's now the archived half. Per author, kept the diagnostic plots (raw-vs-
translated list-length scatter, species×stress coverage heatmap, per-stress UpSet plots,
observed-vs-random-background consistency barplots, final UpSet of the 4 consistent
gene sets) and archived a large trailing envPC-GWAS/RAO/developmental-gene-overlap
enrichment analysis to `archived/10B_stressGeneEnrichmentAnalysis.ipynb` (unparameterized,
depends on objects built earlier in the active notebook, matching the archived-code
convention).

Three real bugs fixed in the kept half: (1) the metadata load pointed at
`p_phyloGWAS_archived/data/DEG_study_metadata.csv`, which doesn't exist there at all —
`data/DEG_study_metadata.csv` (and every DE-gene-list file it references) lives in this
repo, matching `DATA.md`'s own documented location; (2) the observed-vs-random-background
consistency check (kept per above) sampled from a `background` variable that was only
ever assigned much later, inside the now-archived enrichment section — would have left
`background` undefined once that section moved out; per author, now loads the stage-03
filtered-OG universe (`output/poaceaeHelixerOG_filtered_20250331.txt`) as the null-sampling
background instead; (3) immediately after writing the real output, a follow-up cell
re-read a stale copy of the same file from the sibling `p_phyloGWAS_archived` project
instead of just reusing the object already in memory — removed the unnecessary,
wrong-source re-read, keeping only the UpSet visualization. A separate, now-broken
worked-example block (cells validating the ID-conversion heuristic against one
hardcoded study row) was dropped rather than patched: that row's CDS.fa column is empty
in the current metadata file, and nothing downstream used its output anyway — the real,
generalized conversion (with a proper empty-CDS.fa guard) was already correct. Verified:
real-data run of the fully-fixed active notebook reproduces the existing
`output/candidateGenes/consistentEnvResponsiveGenes.json` exactly (same 4 stress
categories, same gene-set membership, same counts — 2212/2874/1292/2524 for
cold/heat/drought/waterlogging) — confirms the metadata-path fix doesn't change any
actual result, just where it's correctly read from.

**Update:** `10A_DEG_IDconversion.sh`'s already-disabled "pot-in-pot dataset" block (a
Helixer-annotation-based variant of the same miniprot ID-conversion step, for a different
assembly set under `output/candidateGenes_remob/assemblies_pp/`) — commented out by the
author rather than actively run — moved to
`archived/10A_potInPot_DEG_IDconversion.sh`, kept commented out exactly as found, matching
the archived-code convention rather than left disabled in place.

**Note on 11:** `11_candidateOGInvestigation.ipynb` (120 cells) integrated the paper's core
candidate-identification story — phylogenetic-mixed-model p-value distributions (Fig 5a/5b),
enrichment of mixed-model hits against a-priori DE genes and GO terms (BP/CC, Fig 5c/5d), and
the 3-layer overlap (empirical p < 0.001 AND a-priori DE AND HyPhy RELAX significant) yielding
the final candidate OGs — with two side analyses: per-residue association modeling for two
already-validated candidates (PGAM/OG0018915, EXP5A/OG0001399), and an untracked, ad hoc lookup
checking where a collaborator-supplied aminotransferase gene list falls relative to the
analysis. Per author, the two side analyses don't belong in the core notebook. The per-residue
modeling — real, substantial supplementary analysis on named candidates, not throwaway — moved
to a new **active** sibling notebook, `11B_candidateGeneResidueModeling.ipynb` (confirmed via
full-notebook grep to have zero in-memory dependency on the core notebook — self-contained
aside from the `phyloKMat` load carried over from the old shared header). The aminotransferase
lookup — a closed, one-off collaborator request — moved to
`archived/11X_aminotransferaseLookup.ipynb` (unparameterized, depends on objects built earlier
in the core notebook, matching the archived-code convention). Also dropped the old notebook's
setup cell entirely: it loaded `metadata`/`envData`/the species tree/`commonID`, none of which
were referenced anywhere else in the notebook — only the `phyloKMat` load it also contained was
actually used, and only by the per-residue modeling now split out.

**A real, consequential bug was found and fixed** in the core notebook's 3-layer overlap: the
HyPhy-significance layer for envPC2 (drought/wet) and envPC3 (sand/clay) was built as
`union(hyphy_drought$OG[...], hyphy_drought$wet[hyphy_wet$logp>sigCutOff])` and the sand/clay
equivalent — but `RELAX_resultTable_drought`/`_clay` have no `wet` column (real columns are
just `LRT, p-value, k, logp, OG`), so `hyphy_drought$wet`/`hyphy_clay$wet` silently evaluated to
`NULL` in R. This meant the wet- and clay-significant OGs were **entirely missing** from the
envPC2/envPC3 HyPhy layer — confirmed `hyphy_wet$OG`/`hyphy_clay$OG` were never referenced
anywhere else in the notebook. Fixed to `union(hyphy_drought$OG[hyphy_drought$logp>sigCutOff],
hyphy_wet$OG[hyphy_wet$logp>sigCutOff])` and the sand/clay equivalent, mirroring the already-
correct cold/warm pattern one cell above it. Verified against the real RELAX result tables:
`hyphyCandidate2` (envPC2) grows from 19 to 28 OGs, `hyphyCandidate3` (envPC3) from 5 to 13 OGs.
**This changes the final candidate-OG count from the previously-reported 17 to 27**
(envPC1 layer unaffected at 10; envPC2's 3-layer overlap goes 6→10; envPC3's goes 1→7) — this
is a real change to a reported scientific result, not a formatting fix; flagged prominently for
the author to review against the manuscript.

**Update:** the notebook had never actually written out the final candidate OGs with their
per-candidate statistics (phylogenetic mixed model p-value, RELAX significance, a-priori DE
conditions) — the "candidate gene identification" cells only ever printed comma-joined gene-ID
strings inline, nothing was saved. Added new cells building `output/finalCandidateOGs_summary.txt`
(one row per final candidate: `OG, envPC, ZmID, OsID, phylo_p, phylo_empP,
RELAX_dir1_condition/_logp/_k, RELAX_dir2_condition/_logp/_k, aprioriDE_conditions` — the two
RELAX "dir" columns hold whichever pair of opposite-direction tests apply to that OG's envPC
axis, e.g. cold/warm for envPC1, so the table stays one unified schema across all three axes)
and ran it for real — confirmed 27 rows, matching the corrected count above. Also updated
`src/fig6a_sankey_v2.py`'s hardcoded `SELECTION`/`EXPRESSION` dicts (the "selection shift" and
"differential expression" Sankey columns) from the pre-fix values (30/19/5 and 10/6/1) to the
corrected ones (30/28/13 and 10/10/7, re-verified directly from real data) and re-ran the
script to confirm it still renders correctly.

**Update:** the core notebook has since been renamed `11_candidateOGInvestigation.ipynb` →
`11A_candidateOGInvestigation.ipynb`, matching the `A`/`B` lettering already used by its
sibling `11B_candidateGeneResidueModeling.ipynb`.

**Note on figure numbering:** the manuscript's figure order (Fig1–6 + FigS1–S6) had shifted
independently of the code since much of it was last run, leaving several `output/figure/*.png`
filenames stale relative to their actual manuscript identity. Per author review, reconciled
across 7 notebooks and 3 scripts, keeping every existing `_v2`/`_v3`/etc. version suffix intact
(only the figure *number*, and one capitalization, changed — an older, unversioned filename may
have no producing code left at all, so the version marker itself was never in question): 08A's
`Fig3c_v2.png`/`Fig3d.png` (LLM score validation) → `Fig4a_v2.png`/`Fig4b.png`; 08C's
`Fig4a_v3.png`/`Fig4b_v4.png`/`Fig4c_v6.png` (physicochemical panels) → `Fig3a_v3.png`/
`Fig3b_v4.png`/`Fig3c_v6.png`; 08A's `Fig2.png` (supplemental QC histogram) → `FigS3.png`;
`suppFig_genomicPipeline_*` (03B/03D/08A) → `FigS5_genomicPipeline_*`; `suppFig_envPCPipeline_*`
(06B / `src/S05_envPC_analysis.R`) → `FigS1_envPCPipeline_*`; `powerSimulation_XY_revised.R`/
`_vis.R`'s `Fig8_revised.png` → `Fig2_revised.png`; `src/fig6a_sankey_v2.py`'s
`fig6a_sankey_v2.png`/`.svg` → `Fig6a_sankey_v2.png`/`.svg` (capitalization only). 06B's
`FigSX_transition.png` (ancestral-state transition-node plot) was confirmed **not in the
paper at all** — renamed to the non-numbered `envPC_ASR_transitionNodes.png` rather than
given a stale placeholder number. Also re-enabled a disabled save call for `Fig1a_v2.png`
(06B cell 35 — the plot was still being built/shown, just never written to disk, apparently
by accident) and dropped one stale dead-code cell (an unsaved duplicate plot under a leftover
"Figure 5B" comment in 11). Every existing real output file was physically renamed to match
(backing up any pre-existing file that already sat at the target name).

`output/figure/` also carries a number of known-orphaned files — superseded reruns/versions
with no producing code anywhere in the active pipeline, left in place (gitignored, low risk,
not worth deleting without being asked): `Fig1a.png`, `Fig1a_map.png`, `Fig1c.png`,
`Fig1c_tree.png` (superseded by their `_v2` versions), `Fig3b.png`, `Fig3c.png` (superseded by
`_v4`/`_v6`), `Fig4a.png` (superseded by `_v2`), `Fig4b_v2.png`, `Fig4b_v3.png` (superseded by
the current unversioned `Fig4b.png`), `Fig4c.png`/`Fig4c_v2–v5.png`/`Fig4c_S1.png`/
`Fig4_S1.png` (all superseded by 08C's `Fig3c_v6.png`), `Fig5d.png` (superseded by `_v2`),
`FigS6.pdf`/`FigS6_check.pdf` (superseded by `FigS6_v2.pdf`), `fig6a_sankey.png`/`.svg`
(pre-`_v2`, lowercase), `suppFig/suppFig_genomicPipeline.pptx` (stray, never code-generated),
`suppFig/suppFig_envPCPipeline_d.png`/`_e.png` (bare, no `_v2` ever existed for these two).

**Note on the life-history/rhizome thread's extraction:** the recurring annual/perennial
life-history and rhizome (underground-stem presence/absence) side-thread — flagged throughout
stages 06–10's notes above as "slated to move to its own repo" — has now moved, to a new
sibling repo, `/local/workdir/sh2246/p_lifeHistoryPoaceae` (fresh `git init`, no history
preserved). Moved: `archived/06B_supplementalLifeHistory.ipynb`,
`archived/07Ab_prematureStopCodon_rhizomeLifeHistory.ipynb`,
`archived/08X_binomialLifeHistoryModel.ipynb`, `archived/08B_GOenrichmentForAnnualism.ipynb`,
`archived/09A_lifeHistoryRhizome_HyPhyPipeline.sh`,
`archived/10E_candidateOGInvestigation_lifeHistory.ipynb`, `src/01_runGWAS.R` (removed from
this repo's script library entirely), and 3 non-functional early prototypes
(`XX_archived/example_GWAS_MVP_MixLog.R`, `XX_archived/06D_lifehistoryOnMap.R`,
`XX_archived/18_pubTraitDB.ipynb`, kept there for lineage only). The one piece still embedded
in active code — `09B`'s 3 life-history/rhizome cells — was extracted into a new notebook
there (see "Note on 09A/09B" above).

The new repo vendors ~1GB of small core data/output files as real copies (its own README.md
lists every one); the 48GB `output/masterDataTable_PAVFill_20251001.txt` is symlinked rather
than copied (both repos read the one real file in this repo's `output/`); ~74GB of upstream
HyPhy intermediates (`output/CDSMSAPerOG_gs/`, `output/CDSMSAPerOG_HyPhy_20250203/`,
`output/geneTree_allOGs_20250203/`, and their life-history-specific derived siblings) were
**not** vendored — only the ~360MB of already-computed `output/HyPhyResult/{lifeHistory*,
rhizome*}/` results moved. `08A_masterDataTableGeneration.ipynb`'s `rhizomeDat`/`lifeHistory`
merge (cells 2/17) stays in this repo, since it feeds the shared master table both the core
analysis and the moved thread read from — left as harmless passenger columns, confirmed
unused by `08B`/`08C` themselves. Every known gap the new repo's own inventory surfaced
(several files missing entirely, including a whole third sibling repo `p_panAndOGASR`; a
handful of orphan data files with no producer script found anywhere; one 0-byte HyPhy RELAX
json) is documented in that repo's README.md rather than fixed or silently worked around here.

---

## Tools used across stages

Software/tools are tracked here, not in `DATA.md` (which is data only). Conda environments
that provide these are in `envs/` (see README.md's "Environment setup" section); versions
below are confirmed (from a hardcoded invocation path in the script itself, or directly with
the author) — not guessed.

| Tool | Version | Used in stage(s) | Provided by |
|---|---|---|---|
| megahit | — (external repo) | 01 (assembly; external pipeline, see stage README) | not applicable — run outside this repo |
| patch-scaffolds | 0.2.1 | 04 (extract per-assembly orthologous CDS/protein sequences via miniprot GFFs, see stage README) | vendored source/jar, `src/patch-scaffolds` (requires Java 17+; compile/build instructions in that directory) |
| mafft | 7.520 | 04 (per-OG MSA generation) | `envs/environment-tools.yml` |
| OrthoFinder | 2.5.4 | 03 (orthogroup construction) | `envs/environment-tools.yml` |
| miniprot | 0.13 | 03 (cross-mapping), 05 (tree construction), 10 (DEG ID conversion) | `envs/environment-tools.yml` |
| RAxML | 8.2.12 | 05 (gene trees), 09 (gene trees for RELAX) | vendored source, `src/standard-RAxML` (compile from source; not a conda package) |
| astral-pro (`aster` package) | 1.16 | 05 (main species tree) | `envs/environment-tools.yml` |
| ASTER (astral-pro3, vendored) | commit `6df009e` | 05 (secondary exploratory tree only; the invocation in 05A is commented out) | vendored source, `src/ASTER` (compile from source; not a conda package) |
| HyPhy | 2.5.49 | 09 (RELAX) | `envs/environment-tools.yml` |
| seqkit, gffread, GNU parallel | 0.15.0, unpinned, unpinned | utility use across multiple stages | `envs/environment-tools.yml` |
| R | 4.2 | 02, 06, 08, 09, 11 (stats/envirotyping/modeling/visualization) | `envs/environment-r.yml` |
| ASReml-R (licensed) | 4.2.0.302 | 08 (per-OG model, power simulation) | not conda — see README.md's "Licensed software" section |
| Python, PyTorch, transformers | 3.11, 2.1.1, 4.40.0 | 07 (ESM2/PlantCAD zero-shot scoring) | `envs/environment-py.yml` |

---

## Running a stage

1. Check `DATA.md` for that stage's row(s) — confirm the **Path** exists (T1 items should
   already be in place; T2 items can be taken as-is or regenerated by rerunning the
   upstream stage that produces them).
2. Activate the conda env this stage's tools live in — check the **Provided by** column
   above (e.g. `conda activate phyloGWAS-tools` for stages 03/04/09/10, `phyloGWAS-r` for
   02/06/08/09/11, `pytorch-2-1` for 07). Open the notebook(s)/script(s) listed under
   **Key scripts** for that stage.
3. Set `export PHYLOGWAS_ROOT=/path/to/your/clone` before running any stage **script**
   (`.sh`/`.R`/`.r`) — they all read this var, falling back to `/workdir/sh2246/p_phyloGWAS`
   if unset, so nothing changes if you're on this machine. **Notebooks (`.ipynb`) still
   hardcode `/workdir/sh2246/p_phyloGWAS/` directly** — that pass hasn't been done yet;
   edit those paths by hand if running a notebook from a different mount point.
4. SCINET-only steps (`07Ca` ESM2, `07Da`/`07Db`/`07Dc` PlantCAD) require a SLURM allocation
   on USDA SCINET Atlas (`buckler_lab_panand` account) and a GPU node, and are not runnable on
   this machine directly — their outputs are already present locally (see `DATA.md`), except
   PlantCAD's, which hasn't been run for real yet (no combined single-fasta CDS input exists
   on disk — see the `TODO` in `07Da_PlantCAD_SCINET.sh`).
5. Stage 01 (assembly) is not runnable from this repo at all — it's maintained in a
   separate repository ([bucklerlab/p_reelgene](https://bitbucket.org/bucklerlab/p_reelgene/src/master/short_read_assembly/));
   this repo picks up downstream of its output (`data/assemblies/`).
