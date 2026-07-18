# WORKFLOW.md — Analysis Pipeline, Stage by Stage

The pipeline runs as a numbered sequence of stages (`01` → `10`). Each row below maps in
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
| **03** `orthogroup` | OrthoFinder (32 representative genomes) → orthogroups → ancestral AA sequence reconstruction → miniprot cross-mapping → OG filtering | Genome assemblies; dN/dS table; rice→OG mapping; OG→maize mapping; miniprot GFF annotations | `03A_buildHelixerOG.sh`, `03B_OGFilter.ipynb`, `03C_miniProtResult_eval.ipynb`, `03D_OGtranslation.R` | Filtered OG list + ancestral sequences → 04, 07 |
| **04** `msaGeneration` | Per-OG CDS multiple sequence alignment (mafft) | Filtered OG list + CDS sequences (from 03) | `notebook/04_msaGeneration/README.md` — `mafft --ep 0 --genafpair --maxiterate 1000 <input> > <output>` | Per-OG MSAs (`output/OrthofinderMAFFT/*_mafft.fa`) → 05 |
| **05** `phylotreeConstruction` | Gap-strip CDS MSAs → RAxML gene trees → ASTRAL-Pro species tree → phylogenetic K (relatedness) matrix | Angiosperms353 Oryza reference; maize v5 mRNA reference (regenerates OG→maize mapping and the species-name list as a side effect) | `05A_treeConstruction`, `05B_neutralPhylogenyVisualization.ipynb` | Species tree + phyloK matrix → 08, 09 |
| **06** `envirotyping` | Species occurrence coordinates → WorldClim/soil rasters → habitat summary → envPC1–3 | Species-name list; GBIF/BIEN occurrence records; WorldClim + soil rasters; environmental metadata; derived occurrence dataset (Zenodo) | `06B_spCoordEnvData.sh`, `06C_visualizationEnvAdapt.ipynb`, `06D_supplFig_envPCpipeline.R` | envPC1–3 table (Fig. 1) → 08 |
| **07** `summaryStats` | Per-OG premature-stop/frameshift calling, tip-to-outgroup dN/dS, ESM2 & PlantCAD zero-shot scores | miniprot GFF annotations; OrthoFinder protein MSAs; ESM2 weights; PlantCAD zero-shot scores | `07Aa`/`07Ab`; `07Ba`/`07Bb` (SCINET); `07Ca`/`07Cb`/`07Cc` (SCINET GPU) | Per-OG activity scores (Fig. 4) → 08 |
| **08** `linearModeling` | Master data table → genome-wide feature association (Fig. 3) → per-OG phylogenetic mixed model + permulation (Fig. 5) → power simulation (Fig. 2) | dN/dS table; OG→maize mapping; maize v5 expression (FPKM); ASReml-R (licensed) | `08A_masterDataTableGeneration.ipynb`, `08B_genomicFeatureAssociation.ipynb`, `08C_perOGmodel.sh`, `08D_power_simulation.sh` | Candidate-OG lists + model results → 09, 10 |
| **09** `molEvolution` | MSA cleaning → RAxML gene trees → foreground/background branch labeling → HyPhy RELAX selection-intensity tests per trait | Maize v5 GO annotation; OG→maize mapping; HyPhy, RAxML (tools) | `09A_HyPhyPipeline.sh`, `09B_RELAX_resultSummary.ipynb` | RELAX result tables → 10 |
| **10** `aprioriCandidate` | OG→gene-ID mapping (Helixer) → integrate ASReml + RELAX + expression evidence → final candidate OG list (Fig. 6 Sankey) | Per-species CDS FASTAs (stress genes); rice→OG mapping; DeepGO GO annotation; DEG study metadata; maize/switchgrass GO tables | `10A_DEG_IDconversion.sh`, `10B_candidateOGInvestigation.ipynb` | 17 high-confidence candidate OGs (final) |

`slurm/`, `XX_archived/`, and `*/archived/` subfolders hold SLURM job templates and
superseded/exploratory notebooks — not part of the active sequence above.

---

## Tools used across stages

Cross-referencing `DATA.md`'s "Licensed / external tools" table to where each is actually
invoked:

| Tool | Used in stage(s) |
|---|---|
| megahit | 01 (assembly; external pipeline, see stage README) |
| mafft | 04 (per-OG MSA generation) |
| miniprot | 03 (cross-mapping), 05 (tree construction), 10 (DEG ID conversion) |
| RAxML | 05 (gene trees), 09 (gene trees for RELAX) |
| HyPhy | 09 (RELAX) |
| ASReml-R (licensed) | 08 (per-OG model, power simulation) |
| seqkit, gffread, GNU parallel | utility use across multiple stages |

---

## Running a stage

1. Check `DATA.md` for that stage's row(s) — confirm the **Path** exists (T1 items should
   already be in place; T2 items can be taken as-is or regenerated by rerunning the
   upstream stage that produces them).
2. Open the notebook(s)/script(s) listed under **Key scripts** for that stage.
3. Hardcoded paths in notebooks/scripts assume the repo lives at
   `/workdir/sh2246/p_phyloGWAS/` — adjust if running from a different mount point.
4. SCINET-only steps (07Ba/07Bb, 07Ca/07Cb/07Cc) require a SLURM allocation on USDA
   SCINET Atlas (`buckler_lab_panand` account) and are not runnable on this machine
   directly — their outputs are already present locally (see `DATA.md`).
5. Stage 01 (assembly) is not runnable from this repo at all — it's maintained in a
   separate repository ([bucklerlab/p_reelgene](https://bitbucket.org/bucklerlab/p_reelgene/src/master/short_read_assembly/));
   this repo picks up downstream of its output (`data/assemblies/`).
