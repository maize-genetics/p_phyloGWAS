# DATA.md — Input-Data Manifest

This manifest inventories every external input needed to reproduce the phylogenetic-GWAS
pipeline in this repository, from raw genome assemblies through to the final figures. It
records what each input is, where it currently lives (verified against disk on this
machine), where it canonically comes from, and which pipeline stage(s) consume it.

Repo root: `/local/workdir/sh2246/p_phyloGWAS`. Sibling repo referenced below:
`/local/workdir/sh2246/p_evolBNI`. Note that in code (notebooks/scripts) these are
hardcoded with paths of the form `/workdir/sh2246/...` (no `/local` prefix) — that's a
machine-mount detail, not a typo.

**Status: reviewed.** The author annotated every row of the first draft of this manifest
(via a "Reviewer action" column) and the corrections below have been folded in: unneeded
rows removed, wrong-version references corrected, and the files the review flagged for
consolidation have been physically moved from `p_phyloGWAS_archived` into this repo
(`data/`, `output/`, or `src/`, as specified per item) — including updating every active
notebook/script that hardcoded the old path. See "Resolved in this pass" at the bottom for
the full list of what changed and why.

## Tier legend

- **[T1]** — needed to reproduce the full pipeline end-to-end starting from raw genome
  assemblies through to the final figures.
- **[T2]** — needed only to reproduce downstream stages from the provided intermediate
  files already checked into `output/` (i.e. you can skip the raw/upstream inputs and
  start partway through).

One general note up front: the raw sequence/assembly release (Section 1) is currently
under a **publication-timed embargo** per the manuscript's Data Availability Statement —
this is expected, intentional, and not a repo bug. The gap self-resolves at publication,
when the assemblies are deposited to SRA / Ag Data Commons.

---

## 1. Genome assemblies (raw) [T1]

| File/dataset | What it is | Current location | Canonical source | Used by (stage) |
|---|---|---|---|---|
| 727 × `*.fa.gz` genome assembly files (504 GB total) | Raw genome assemblies, one per accession, across the Poaceae sampling panel | `data/assemblies/` (verified: directory exists, 727 files, 504 GB) | SRA + USDA Ag Data Commons — **embargoed until manuscript publication** | Upstream of all stages (orthology inference, phylogenetics, molecular evolution, GWAS-style modeling) |

Confirmed by the author: keep as-is; these will be uploaded to Ag Data Commons at
publication.

---

## 2. Repo `data/` files (present) [T1/T2]

All paths below were verified present with `test -f`/`ls` under `data/`. Rows the author
flagged as unneeded (superseded versions, or belonging to a different project sharing this
workdir) have been removed from this manifest — the underlying files were **not** deleted
from disk, since some are shared with other projects.

| File/dataset | What it is | Current location | Canonical source | Used by (stage) |
|---|---|---|---|---|
| `Poaceae_metadata_highErrorFiltered_2025.10.08.tsv`, `Poaceae_metadata_filtered_2025.08.28.tsv` | Per-accession sample/taxonomy metadata. Only these two versions are actually used; `highErrorFiltered_2025.10.08` is a QC-error-rate-filtered subset of `filtered_2025.08.28`. Older dated versions (`2024.08.21`, `2025.04.02/03/04/07/16`) are superseded and not tracked here. | `data/` | Compiled in-house from assembly panel + QC (see stage 02A) | Metadata processing, taxon selection |
| `fullSetOGs_240903.txt` | **This is the dN/dS table for the full orthogroup set**, not just an OG list (corrected per author; the file is not renamed on disk since it's read by name — `data/fullSetOGs_240903.txt` — by stage 08A's master-table notebook, and a bare rename would break that without also updating the consumer) | `data/` | Derived in-house (dN/dS calculation across all OGs) | Orthology/dN-dS-based stages (03B, 08A, 08B, 09B, 10B) |
| `maxFPKM_v5.csv` | Maize v5 max-expression (FPKM) table | `data/` | MaizeGDB / Phytozome (maize B73 v5 expression atlas) | Expression-informed candidate filtering |
| `Zm00001eb.1.fulldata_goList.txt` | Maize B73 v5 GO annotation list | `data/` | MaizeGDB/Phytozome (maize v5 annotation) | GO enrichment (09B, 10B) |
| `Pv_GOTable.txt` | Switchgrass (*Panicum virgatum*) GO table | `data/` | Phytozome (*P. virgatum*) | GO enrichment |
| `Zmays_cds.fa` | Maize reference CDS FASTA | `data/` | MaizeGDB/Phytozome | Sequence-based stages (orthology, molecular evolution) |
| `env_metadata.txt` | Environmental metadata table | `data/` | Compiled in-house (envirotyping pipeline) | Envirotyping (stage 06) |
| `busco_gene_OG_20250828.txt` | BUSCO gene ↔ orthogroup mapping (QC/orthology cross-check only, not part of the main analysis) | `data/` | Derived from BUSCO run (in-house) | QC / orthology cross-check |
| `maize_busco.tsv` | Maize BUSCO completeness table (QC only, not part of the main analysis) | `data/` | Derived from BUSCO run (in-house) | QC |

**The actual phylogenetic backbone used is not this section's external Kew tree** — see
the correction below.

> ~~`angiosperm353_ASTRALPRO3.greaterthan50genespersample.tre`~~ — **removed.** The author
> confirmed this published Angiosperms353/Kew reference tree is *not* what was used for the
> analysis backbone. The actual phylogenetic backbone is a pipeline-generated intermediate,
> not an external download: **`output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk`**
> (produced by stage 05, from the 353 gene trees reconciled via ASTRAL-Pro). Nothing to
> "gather" here — it's regenerated by re-running stage 05 end-to-end.

**Rows removed** (confirmed by the author as unneeded — either superseded versions,
belonging to a different project sharing this workdir, or not part of this analysis; files
left in place on disk, not deleted):
`tabasco_Poaceae_20250402.csv` / `tabasco_Poaceae800.csv` (redundant with the metadata
table above), `Poaceae_taxonomic_groups_subtribeAdded_SKH20250402.csv` (redundant with the
metadata table above), `annual_assemblies_*.txt`, `perennial_assemblies_20260213.txt`,
`perennialRhizome_assemblies_20260303.txt`, `nonrhizomatous_assemblies_20260415.txt` (all:
different project), `grassBase_cleaned_rhizome.txt` (different project),
`ATgene_sorghum_maize_syntelog_orthogroup_mapping.tsv` (an aminotransferase gene
investigation set, not relevant to this paper), `JonathanModule/` SD5/SD8 (different
project — grass senescence network modules), `Zmarina_cds.fa`,
`Zm-LH244-REFERENCE-BAYER-1.0.fa` (not needed).

---

## 3. Moved into this repo from `p_phyloGWAS_archived/` [T1/T2]

These files previously lived only in the sibling `p_phyloGWAS_archived` (not under this
repo's version control). Per the author's review, they have now been **physically moved**
(not copied) into `p_phyloGWAS`, and every active notebook/script that hardcoded the old
archived path has been updated to point at the new location.

| File/dataset | What it is | New location (moved) | Used by (stage) |
|---|---|---|---|
| `fullQCData_20240819.txt` | Per-assembly QC statistics | `data/fullQCData_20240819.txt` | 02A metadata processing |
| `DEG_mappingFiles/` (incl. `OG_mapping_Oryza_sativa.IRGSP-1.0.cds.all.txt`) | Rice → orthogroup mapping files | `data/DEG_mappingFiles/` | 03B, 10A, 10B |
| `OGToZm_mapping_v2.txt` | Orthogroup → maize gene mapping | `data/OGToZm_mapping_v2.txt` | 03B, 05A (regenerates it), 08B, 09B, 10B |
| `orthofinderMiniProt/`, `orthofinderMiniProt_additionalOGs/` (801 + 801 GFF files, 13 GB) | Per-orthogroup miniprot GFF annotations | `output/orthofinderMiniProt/`, `output/orthofinderMiniProt_additionalOGs/` | 03C, 07Aa |
| `poaceaeHelixerOG_ancSeq_gapRemoved_v2_uppercase_output_CC_f.gmt` | DeepGO GO-annotation GMT file | `data/poaceaeHelixerOG_ancSeq_gapRemoved_v2_uppercase_output_CC_f.gmt` | 10B |
| `DEG_study_metadata.csv` | Curated DEG/abiotic-stress study metadata | `data/DEG_study_metadata.csv` | Stage 10 (candidate gene / expression evidence) |
| `DEG_CDS_fasta/` (23 files, 626 MB) | Per-species CDS FASTAs for a-priori stress-responsive gene sets | `data/DEG_CDS_fasta/` | 10A |
| `spNameMetadata_20240819.txt` | Species-name list for coordinate lookup | `data/spNameMetadata_20240819.txt` | 06B (02A also regenerates this file directly) |
| `Angiosperms353_orysaSequences.fasta` | Angiosperms353 Oryza reference sequences | `data/Angiosperms353_orysaSequences.fasta` | 05A (angiosperm353 cross-mapping) |
| `Zea_mays_v5_mrna.fa` (128 MB) | Maize B73 v5 mRNA sequences | `data/Zea_mays_v5_mrna.fa` | 05A (OG→Zm cross-mapping) |
| `asreml_4.2.0.302_linux-intel64_R4.2.0.tar.gz` (22 MB) | ASReml-R v4.2.0 installer archive | `src/asreml_4.2.0.302_linux-intel64_R4.2.0.tar.gz` (matches how other vendored tools, e.g. RAxML, are organized — see Section 9) | 08C, 08D |

**Active-code references fixed as part of this move** (each hardcoded the old
`p_phyloGWAS_archived/...` path; all now point at the new local path above):
`notebook/02_metadataCuration/02A_metadataProcessing.ipynb`,
`notebook/03_orthogroup/03B_OGFilter.ipynb` (2 occurrences),
`notebook/03_orthogroup/03C_miniProtResult_eval.ipynb` (3 cells),
`notebook/08_linearModeling/08B_genomicFeatureAssociation.ipynb`,
`notebook/09_molEvolution/09B_RELAX_resultSummary.ipynb` (also fixed its
`p_panAndOGASR`-sourced GO-list read to use the local `data/Zm00001eb.1.fulldata_goList.txt`
mirror instead), `notebook/10_aprioriCandidate/10B_candidateOGInvestigation.ipynb` (edited
only the specific affected lines; its other in-progress edits were left untouched).
`notebook/05_phylotreeConstruction/05A_treeConstruction`,
`notebook/07_summaryStats/07Aa_frameShiftMutation.ipynb`, and
`notebook/10_aprioriCandidate/10A_DEG_IDconversion.sh` already used relative/local paths
and needed no change.

`bien_coordinates_clean_2023.12.04.csv` — this file was already inside this repo's own
`output/envData/speciesRange/` (an earlier pass had mis-attributed it to the archived
sibling; confirmed resolved by the author). One additional consumer,
`notebook/06_envirotyping/06D_supplFig_envPCpipeline.R`, was found during verification
still hardcoding the nonexistent archived path for this file (a pre-existing bug,
unrelated to this pass's moves) — fixed to point at the local `output/` copy.

---

## 4. `p_evolBNI/` GIS rasters (present, archived independently) [T1]

Sibling repo `p_evolBNI`, directory `p_evolBNI/data/GIS_env_data/`. Used by envirotyping
stages 06B/06D via `src/09_pulling_envData.r` (and variants).

| File/dataset | What it is | Current location | Canonical source | Used by (stage) |
|---|---|---|---|---|
| `wc2.1_2.5m_elev.tif` | WorldClim 2.1 elevation raster (2.5 arc-min) | `p_evolBNI/data/GIS_env_data/WorldClim_raw_2.5m_files/wc2.1_2.5m_elev/wc2.1_2.5m_elev.tif` | WorldClim 2.1 | 06B/06D |
| `wc2.1_2.5m_bio_1.tif` (+ other bioclim variables) | WorldClim 2.1 bioclimatic variable raster(s) | `p_evolBNI/data/GIS_env_data/WorldClim_raw_2.5m_files/wc2.1_2.5m_bio/wc2.1_2.5m_bio_1.tif` | WorldClim 2.1 | 06B/06D |
| `HYSOGs250m.tif` | Global Hydrologic Soil Group raster (250 m) | `p_evolBNI/data/GIS_env_data/Global_Hydrologic_Soil_Group_1566/Global_Hydrologic_Soil_Group_1566/data/HYSOGs250m.tif` | NASA / ORNL DAAC Global Hydrologic Soil Groups (1566) | 06B/06D |
| `GSDE_raw_nc_files/*.nc` | Global Soil Dataset for Earth system modeling (33 per-variable NetCDF files) | `p_evolBNI/data/GIS_env_data/GSDE_raw_nc_files/` | GSDE | 06B/06D |

Confirmed by the author: **left in place, not moved.** `GIS_env_data/` will be archived
independently and linked directly from the paper, so it doesn't need to be consolidated
into `p_phyloGWAS`.

---

## 5. Off-machine siblings and legacy path references [T1]

| Item | Status | Detail |
|---|---|---|
| `p_geneModelEvaluation/` | Confirmed absent on this machine; **resolved** | `src/S00_basicFunction.R` was vendored into this repo (confirmed by the author), and the one remaining active reference to the missing sibling — in `notebook/09_molEvolution/09B_RELAX_resultSummary.ipynb` — has now been fixed to `source("/workdir/sh2246/p_phyloGWAS/src/S00_basicFunction.R")`, matching the pattern already used by `10B`/`10E`. |
| `Zm00001eb.1.fulldata_goList.txt`, `Pv_GOTable.txt` | No action needed (confirmed by the author) | Legacy notebooks reference these under the missing `p_panAndOGASR/data/`, but both files are already mirrored into `data/` (Section 2) — the active pipeline reads the local copies. |

**Removed:** the `p_panAndOGASR/` row — the author confirmed none of the files it hosted
(`sixDigitCode2AnnotID`, `phyloP_*` results, `OG_4-8annual17perennial.csv`, referenced only
by the legacy `07Ab_prematureStopCodon.ipynb`) are useful for the current analysis. This
sibling repo is not present on this machine and nothing in the active pipeline needs it.

---

## 6. External databases / live queries [T1]

| Item | What it is | Access pattern | Canonical source | Used by (stage) |
|---|---|---|---|---|
| GBIF occurrence records | Species occurrence points | Live query via `rgbif::occ_search`; local intermediate at `output/metadataFormalOut/` | GBIF.org; deposited on Zenodo (see below) | 06B (`src/08_pulling_geo_data.R`) |
| BIEN occurrence records | Species occurrence points | Live query via `BIEN::BIEN_occurrence_species`; local intermediate at `output/metadataFormalOut/` | BIEN; deposited on Zenodo (see below) | 06B (`src/08_pulling_geo_data.R`) |
| Derived occurrence dataset (citable) | Cleaned/combined GBIF+BIEN occurrence dataset used in the manuscript — specifically `formal_envData_20240820.txt` | N/A (deposited dataset) | **Zenodo 14968186** (per manuscript) | 06B/06D |
| `get_spatial_fun.R` | Helper function for spatial/environmental extraction | Vendored into this repo at `src/get_spatial_fun.R` (fetched per the author's request, so it no longer depends on a live fetch at runtime) | Originally `https://raw.githubusercontent.com/gcostaneto/envirotypeR/main/R/get_spatial_fun.R` | `src/09_pulling_envData.r` |

**Resolved:** the Zenodo DOI discrepancy flagged previously — `src/08_pulling_geo_data.R`'s
code comment cited `10.5281/zenodo.14967966`, which didn't match the manuscript's stated
`14968186`. Fixed to match the manuscript's version.

---

## 7. Model weights [T1]

| Model | What it is | Access pattern | Canonical source | Used by (stage) |
|---|---|---|---|---|
| ESM2 (`facebook/esm2_t33_650M_UR50D`) | Protein language model, used for variant-effect logits | Downloaded from HuggingFace Hub at runtime (GPU) | HuggingFace (`facebook/esm2_t33_650M_UR50D`) | 07C (`src/4_ESM_logits.py`) |
| PlantCAD (`kuleshov-group/PlantCaduceus_l32`) | Plant genomic language model | Not loaded directly in the traced active scripts — PlantCAD scores enter modeling only as a precomputed predictor column | HuggingFace (`kuleshov-group/PlantCaduceus_l32`); precomputed scores in `combined_plantCAD_zeroShotScores_20250421.txt` | Modeling stage (predictor input only) |

Confirmed by the author: no action needed for either — both are pulled from HuggingFace as
needed.

---

## 8. SCINET-only raw inputs [T1]

| File/dataset | What it is | Current location | Canonical source | Used by (stage) |
|---|---|---|---|---|
| `orthofinderProteinMSAs/` (+ `_additionalOGs`, now merged into the same directory) | OrthoFinder protein multiple-sequence alignments, by orthogroup | **Already present locally** at `output/orthofinderProteinMSAs/` (6.7 GB, confirmed by the author) — not SCINET-only after all | Originally generated by A. Schulz on SCINET Atlas (OrthoFinder run); `/project/90daydata/buckler_lab_panand/aimee.schulz/panand/output/orthofinderProteinMSAs[_additionalOGs]/` | 07Bb → `src/11_find_premature_stops.py`; 07Cb → `src/4_ESM_logits.py` |

**Deferred, per the author:** `src/11_find_premature_stops.py` and `src/4_ESM_logits.py`
(called from `07Bb`/`07Cb`) still hardcode the SCINET path rather than the local
`output/orthofinderProteinMSAs/` copy that's already present. The author asked to update
this code "later" — intentionally not changed in this pass.

---

## 9. Licensed / external tools [T1]

| Tool | Version | Current location | Notes | Used by (stage) |
|---|---|---|---|---|
| ASReml-R | 4.2.0.302 (linux-intel64, R4.2.0) | Installer moved to `src/asreml_4.2.0.302_linux-intel64_R4.2.0.tar.gz` (Section 3) | **Licensed software (VSNi)** — must be separately licensed/installed by the reproducing party | `src/12_runPermulation_perOGModel.R` (08C), `src/powerSimulation_XY_revised.R` (08D) |
| HyPhy | 2.5.49 | Verified present at `/programs/hyphy-2.5.49/bin/hyphy` | Assumed pre-installed on the compute environment | Molecular evolution (selection tests) |
| RAxML | (standard-RAxML build) | Verified present at `src/standard-RAxML/raxmlHPC` (vendored into this repo) | — | Phylogenetic tree inference |
| seqkit | 0.15.0 | Not verified on this machine (assumed installed per environment spec) | Assumed pre-installed | Sequence manipulation utilities |
| miniprot | — | Not verified on this machine (assumed installed per environment spec) | Assumed pre-installed | Protein-to-genome alignment (03C) |
| gffread | — | Verified present at `/programs/bin/cufflinks/gffread` | — | Annotation/GFF processing |
| GNU parallel | — | Verified present at `/programs/bin/parallel` | — | Pipeline parallelization across many scripts |

Confirmed by the author: no action needed on any of the above beyond the ASReml relocation
already reflected in Section 3. An `environment.yml`/`renv.lock` capturing these versions
is deferred to a later cleanup pass (per the author's own note).

---

## 10. Copied into `data/` (from `p_phyloGWAS_archived/`) [T1/T2]

The remaining files the author marked `copy` (in addition to the items already listed in
Section 3) have been moved into `data/`. Rows marked "not needed" have been removed from
this manifest (source files left in place, not deleted — several are for other projects).

| File | New location | Used by (stage) |
|---|---|---|
| `spNameMetadata_20240819.txt` | see Section 3 | 06B |
| `Angiosperms353_orysaSequences.fasta` | see Section 3 | 05A |
| `Zea_mays_v5_mrna.fa` | see Section 3 | 05A |

**Removed** (confirmed not needed by the author): `allOGs_withHeaders_dnds.txt` (old
version), `Perennitality_Numeric_Poaceae.csv` (different project),
`rhizomePAV_traits.txt` (different project),
`Araport11_seq_20220914_representative_gene_model.fa.gz` (not needed),
`Yeaman2024_RAO_AtID.txt` (not needed).

---

## Resolved in this pass

- **9 files/directories physically moved** from `p_phyloGWAS_archived` into
  `p_phyloGWAS` (`data/`, `output/`, or `src/` per item — Section 3), plus 3 more
  (Section 10), for **12 total moves** (~14 GB). All active notebook/script references to
  the old archived paths were updated to match (7 files edited:
  `02A_metadataProcessing.ipynb`, `03B_OGFilter.ipynb`, `03C_miniProtResult_eval.ipynb`,
  `08B_genomicFeatureAssociation.ipynb`, `09B_RELAX_resultSummary.ipynb`,
  `10B_candidateOGInvestigation.ipynb`).
- **`get_spatial_fun.R` vendored** into `src/` instead of being fetched live from GitHub at
  runtime.
- **Zenodo DOI corrected** in `src/08_pulling_geo_data.R` (`14967966` → `14968186`).
- **`p_geneModelEvaluation` path fixed** in `09B_RELAX_resultSummary.ipynb` to use the
  already-vendored local `src/S00_basicFunction.R`.
- **`p_panAndOGASR/` row removed** — confirmed not needed for the current analysis.
- **Wrong phylogeny reference corrected** — the external Kew Angiosperms353 tree was not
  actually used; the real backbone is the self-generated
  `output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk`.
- **~13 unneeded rows removed** across Sections 2 and 10 (superseded file versions, or
  files belonging to other projects sharing this workdir) — underlying files were **not**
  deleted from disk.
- **`fullSetOGs_240903.txt` description corrected** (it's a dN/dS table, not an OG list) —
  not renamed on disk, since it's read by name by stage 08A and a bare rename would break
  that notebook.

## Still deferred (explicitly, per the author)

- Updating `src/11_find_premature_stops.py` / `src/4_ESM_logits.py` to read the
  now-local `output/orthofinderProteinMSAs/` instead of the SCINET path (Section 8).
- `environment.yml`/`renv.lock` capturing tool versions (Section 9).
