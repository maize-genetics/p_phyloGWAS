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
| **04** `msaGeneration` | Per-OG CDS multiple sequence alignment (mafft) | Filtered OG list + CDS sequences (from 03) | `notebook/04_msaGeneration/README.md` — `mafft --ep 0 --genafpair --maxiterate 1000 <input> > <output>` | Per-OG MSAs (`output/OrthofinderMAFFT/*_mafft.fa`) → 05 (gap-stripping/gene trees), 07 (dN/dS calculation needs the MSA directly), 09 (HyPhy RELAX needs the MSA directly) |
| **05** `phylotreeConstruction` | Gap-strip CDS MSAs → RAxML gene trees → ASTRAL-Pro species tree → phylogenetic K (relatedness) matrix | Gap-stripped CDS MSAs (from 04) | `05A_treeConstruction`, `05B_neutralPhylogenyVisualization.ipynb` | Species tree + phyloK matrix → 08 (predictor); per-OG gene trees → 09 (HyPhy RELAX runs on the gene tree from 05) |
| **06** `envirotyping` | Species occurrence coordinates → WorldClim/soil rasters → habitat summary → envPC1–3 | Species-name list; GBIF/BIEN occurrence records; WorldClim + soil rasters; environmental metadata; derived occurrence dataset (Zenodo) | `06B_spCoordEnvData.sh`, `06C_visualizationEnvAdapt.ipynb`, `06D_supplFig_envPCpipeline.R` | envPC1–3 table (Fig. 1) → 08 |
| **07** `summaryStats` | Per-OG premature-stop/frameshift calling, tip-to-outgroup dN/dS calculation, ESM2 & PlantCAD zero-shot scores | miniprot GFF annotations; OrthoFinder protein MSAs; ESM2 weights; PlantCAD weights | `07Aa`/`07Ab`; `07Ba`/`07Bb` (SCINET); `07Ca`/`07Cb`/`07Cc` (SCINET GPU) | Per-OG activity scores + dN/dS table (Fig. 4) → 08 |
| **08** `linearModeling` | Master data table → genome-wide feature association (Fig. 3) → per-OG phylogenetic mixed model + permulation (Fig. 5) → power simulation (Fig. 2) | dN/dS table (from 07, used as a predictor); OG→maize mapping; maize v5 expression (FPKM) | `08A_masterDataTableGeneration.ipynb`, `08B_genomicFeatureAssociation.ipynb`, `08C_perOGmodel.sh`, `08D_power_simulation.sh` | Candidate-OG lists + model results → 09, 11 |
| **09** `molEvolution` | MSA cleaning (from 04) → RAxML gene trees (from 05) → foreground/background branch labeling → HyPhy RELAX selection-intensity tests per trait | OG→maize mapping | `09A_HyPhyPipeline.sh`, `09B_RELAX_resultSummary.ipynb` | RELAX result tables → 11 |
| **10** `aprioriCandidate` | OG→gene-ID mapping (Helixer) via miniprot | Per-species CDS FASTAs (stress genes) + DEG study metadata; rice→OG mapping | `10A_DEG_IDconversion.sh` | Gene-ID mapping → 11 |
| **11** `candidateOGInvestigation` | Integrate ASReml (08) + RELAX (09) + gene-ID mapping (10) + expression/GO evidence → final candidate OG list (Fig. 6 Sankey) | DeepGO GO annotation; Maize v5 GO annotation | `notebook/11_candidateOGInvestigation/11_candidateOGInvestigation.ipynb` | 17 high-confidence candidate OGs (final) |

`slurm/`, `XX_archived/`, and `*/archived/` subfolders hold SLURM job templates and
superseded/exploratory notebooks — not part of the active sequence above.

**Note on 10/11:** `10A_DEG_IDconversion.sh` and the notebook formerly named
`10B_candidateOGInvestigation.ipynb` used to be presented as one stage. They're distinct
enough in purpose (ID mapping vs. final candidate-list integration) to warrant separate
numbers and separate folders — the notebook has been moved and renamed to
`notebook/11_candidateOGInvestigation/11_candidateOGInvestigation.ipynb` (its in-progress
edits carried over with the move, untouched otherwise).

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
`output/` to `data/` in the move, matching where every active consumer (03E, 08B, 09B, 11,
`src/12_runPermulation_perOGModel.R`) actually reads it from and where the current
(2024-09-16) file lives.

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
| ASTER (astral-pro3, vendored) | commit `6df009e` | 05 (secondary exploratory tree only, in 05B) | vendored source, `src/ASTER` (compile from source; not a conda package) |
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
   02/06/08/09/11, `phyloGWAS-py` for 07). Open the notebook(s)/script(s) listed under
   **Key scripts** for that stage.
3. Set `export PHYLOGWAS_ROOT=/path/to/your/clone` before running any stage **script**
   (`.sh`/`.R`/`.r`) — they all read this var, falling back to `/workdir/sh2246/p_phyloGWAS`
   if unset, so nothing changes if you're on this machine. **Notebooks (`.ipynb`) still
   hardcode `/workdir/sh2246/p_phyloGWAS/` directly** — that pass hasn't been done yet;
   edit those paths by hand if running a notebook from a different mount point.
4. SCINET-only steps (07Ba/07Bb, 07Ca/07Cb/07Cc) require a SLURM allocation on USDA
   SCINET Atlas (`buckler_lab_panand` account) and are not runnable on this machine
   directly — their outputs are already present locally (see `DATA.md`).
5. Stage 01 (assembly) is not runnable from this repo at all — it's maintained in a
   separate repository ([bucklerlab/p_reelgene](https://bitbucket.org/bucklerlab/p_reelgene/src/master/short_read_assembly/));
   this repo picks up downstream of its output (`data/assemblies/`).
