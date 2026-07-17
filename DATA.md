# DATA.md — Input-Data Manifest

This manifest inventories every external input needed to reproduce the phylogenetic-GWAS
pipeline in this repository, from raw genome assemblies through to the final figures. It
records what each input is, where it currently lives (verified against disk on this
machine), where it canonically comes from, and which pipeline stage(s) consume it.

Repo root: `/local/workdir/sh2246/p_phyloGWAS`. Sibling repos referenced below:
`/local/workdir/sh2246/p_phyloGWAS_archived` and `/local/workdir/sh2246/p_evolBNI`. Note
that in code (notebooks/scripts) these are hardcoded with paths of the form
`/workdir/sh2246/...` (no `/local` prefix) — that's a machine-mount detail, not a typo.

## Tier legend / how to use

- **[T1]** — needed to reproduce the full pipeline end-to-end starting from raw genome
  assemblies through to the final figures.
- **[T2]** — needed only to reproduce downstream stages from the provided intermediate
  files already checked into `output/` (i.e. you can skip the raw/upstream inputs and
  start partway through).

One general note up front: the raw sequence/assembly release (Section 1) is currently
under a **publication-timed embargo** per the manuscript's Data Availability Statement —
this is expected, intentional, and not a repo bug. The gap self-resolves at publication,
when the assemblies are deposited to SRA / Ag Data Commons.

## How to annotate

Every table below now ends with a **`Reviewer action`** column for marking corrections
directly in this file, rather than one-by-one in chat. Convention:

- **Blank** — accept the row as-is (default; no change needed).
- **`remove`** — flag the row as unnecessary; optionally add `— <reason>`.
- **`wrong version — use <correct file/path>`** — flag a version/path correction.
- Anything else — free-text notes, questions, uncertainty are all fine too.

Once you've marked up the rows you care about, share the file back (or just tell me which
rows changed) and the corrections will be folded into a cleaned `DATA.md`, with the
published artifact re-rendered to match.

---

## 1. Genome assemblies (raw) [T1]

| File/dataset | What it is | Current location | Canonical source | Used by (stage) | Reviewer action |
|---|---|---|---|---|---|
| 727 × `*.fa.gz` genome assembly files (504 GB total) | Raw genome assemblies, one per accession, across the Poaceae sampling panel | `data/assemblies/` (verified: directory exists, 727 files, 504 GB) | SRA + USDA Ag Data Commons — **embargoed until manuscript publication** | Upstream of all stages (orthology inference, phylogenetics, molecular evolution, GWAS-style modeling) | |

---

## 2. Repo `data/` files (present) [T1/T2]

All paths below were verified present with `test -f`/`ls` under `data/`.

| File/dataset | What it is | Current location | Canonical source | Used by (stage) | Reviewer action |
|---|---|---|---|---|---|
| `Poaceae_metadata_*.csv/tsv` series (multiple dated versions, incl. `Poaceae_metadata_highErrorFiltered_2025.10.08.tsv`, `Poaceae_metadata_filtered_2025.08.28.tsv`) | Per-accession sample/taxonomy metadata, successive filtered/QC'd versions | `data/` | Compiled in-house from assembly panel + QC (see stage 02A) | Metadata processing, taxon selection | |
| `tabasco_Poaceae_20250402.csv`, `tabasco_Poaceae800.csv` | TABASCO(-derived) trait/taxon reference tables | `data/` | Compiled in-house | Taxon/trait curation | |
| `Poaceae_taxonomic_groups_subtribeAdded_SKH20250402.csv` | Taxonomic group assignments (subtribe-level) | `data/` | Compiled in-house | Taxonomic grouping across stages | |
| `annual_assemblies_20250423.txt`, `annual_assemblies_20260213.txt` | Trait lists: annual life-history assemblies | `data/` | Compiled in-house | Trait-based modeling (stage 08+) | |
| `perennial_assemblies_20260213.txt` | Trait list: perennial assemblies | `data/` | Compiled in-house | Trait-based modeling | |
| `perennialRhizome_assemblies_20260303.txt` | Trait list: perennial + rhizomatous assemblies | `data/` | Compiled in-house | Trait-based modeling | |
| `nonrhizomatous_assemblies_20260415.txt` | Trait list: non-rhizomatous assemblies | `data/` | Compiled in-house | Trait-based modeling | |
| `fullSetOGs_240903.txt` | Full orthogroup (OG) list | `data/` | OrthoFinder output (in-house pipeline) | Orthology-based stages (03B, 08B, 09B, 10B) | |
| `grassBase_cleaned_rhizome.txt` | Cleaned rhizome trait data | `data/` | Derived from GrassBase | Trait curation | |
| `maxFPKM_v5.csv` | Maize v5 max-expression (FPKM) table | `data/` | MaizeGDB / Phytozome (maize B73 v5 expression atlas) | Expression-informed candidate filtering | |
| `Zm00001eb.1.fulldata_goList.txt` | Maize B73 v5 GO annotation list | `data/` | MaizeGDB/Phytozome (maize v5 annotation) | GO enrichment (10B) | |
| `Pv_GOTable.txt` | Switchgrass (*Panicum virgatum*) GO table | `data/` | Phytozome (*P. virgatum*) | GO enrichment | |
| `ATgene_sorghum_maize_syntelog_orthogroup_mapping.tsv` | Arabidopsis–sorghum–maize syntelog/orthogroup mapping | `data/` | Compiled in-house from public synteny resources | Cross-species candidate mapping | |
| `angiosperm353_ASTRALPRO3.greaterthan50genespersample.tre` | Angiosperms353-based species tree (ASTRAL-Pro3) | `data/` | Royal Botanic Gardens, Kew (Angiosperms353 project) | Phylogenetic backbone (stage 05+) | |
| `JonathanModule/Supplementary Dataset Grass Senescence - SD5. Leaf network modules.tsv`, `.../SD8. Und network modules.tsv` | Grass senescence gene co-expression network modules | `data/JonathanModule/` (verified: dir present, both files present) | Published supplementary datasets (senescence network study) | Candidate gene / module analysis | |
| `Zmays_cds.fa` | Maize reference CDS FASTA | `data/` | MaizeGDB/Phytozome | Sequence-based stages (orthology, molecular evolution) | |
| `Zmarina_cds.fa` | *Zostera marina* reference CDS FASTA (outgroup) | `data/` | Phytozome | Sequence-based stages | |
| `Zm-LH244-REFERENCE-BAYER-1.0.fa` | Maize LH244 reference genome FASTA | `data/` | Bayer/MaizeGDB reference release | Reference-guided analyses | |
| `env_metadata.txt` | Environmental metadata table | `data/` | Compiled in-house (envirotyping pipeline) | Envirotyping (stage 06) | |
| `busco_gene_OG_20250828.txt` | BUSCO gene ↔ orthogroup mapping | `data/` | Derived from BUSCO run (in-house) | QC / orthology cross-check | |
| `maize_busco.tsv` | Maize BUSCO completeness table | `data/` | Derived from BUSCO run (in-house) | QC | |

---

## 3. `p_phyloGWAS_archived/` inputs (present on this machine, undeclared) [T1/T2]

These live in the sibling repo `p_phyloGWAS_archived`, **not** under version control of
this repo — a fresh `git clone` of `p_phyloGWAS` will **not** bring these along. All paths
verified present with `test -f`/`test -d` except where noted.

| File/dataset | What it is | Current location | Canonical source | Used by (stage) | Reviewer action |
|---|---|---|---|---|---|
| `fullQCData_20240819.txt` | Per-assembly QC statistics | `p_phyloGWAS_archived/output/metadata_processing/fullQCData_20240819.txt` (verified) | Generated in-house (QC pipeline) | 02A metadata processing | |
| `OG_mapping_Oryza_sativa.IRGSP-1.0.cds.all.txt` | Rice → orthogroup mapping | `p_phyloGWAS_archived/data/DEG_mappingFiles/OG_mapping_Oryza_sativa.IRGSP-1.0.cds.all.txt` (verified) | Derived from rice IRGSP-1.0 CDS + OrthoFinder | 03B, 10B | |
| `OGToZm_mapping_v2.txt` | Orthogroup → maize gene mapping | `p_phyloGWAS_archived/output/OGToZm_mapping_v2.txt` (verified) | Derived in-house | 03B, 08B, 09B, 10B | |
| `gffcmp.Zm-B73-REFERENCE-NAM-5.0_helixer.gff.tmap` | gffcompare tmap output (Helixer vs. B73 NAM5 annotation) | `p_phyloGWAS_archived/data/annotations/gffcmp.Zm-B73-REFERENCE-NAM-5.0_helixer.gff.tmap` (verified) | Generated in-house via gffcompare | 03B | |
| `orthofinderMiniProt/*.gff`, `orthofinderMiniProt_additionalOGs/*.gff` | Per-orthogroup miniprot GFF annotations | `p_phyloGWAS_archived/output/orthofinderMiniProt/` (801 files, verified) and `p_phyloGWAS_archived/output/orthofinderMiniProt_additionalOGs/` (801 files, verified) | Generated in-house via miniprot | 03C | |
| `poaceaeHelixerOG_ancSeq_gapRemoved_v2_uppercase_output_CC_f.gmt` | DeepGO GO-annotation GMT file | `p_phyloGWAS_archived/output/deepGO_poaceaeHelixerOG/poaceaeHelixerOG_ancSeq_gapRemoved_v2_uppercase_output_CC_f.gmt` (verified) | Generated in-house via DeepGO | 10B | |
| `DEG_study_metadata.csv` | Curated DEG/abiotic-stress study metadata | `p_phyloGWAS_archived/data/DEG_study_metadata.csv` (verified) | Compiled in-house from published stress-response studies | 10 (archived path) | |
| `DEG_CDS_fasta/` | Per-species CDS FASTAs for a-priori stress-responsive gene sets | `p_phyloGWAS_archived/data/DEG_CDS_fasta/` (verified dir, 23 files) | Compiled in-house from source genome annotations | 10A | |
| `bien_coordinates_clean_2023.12.04.csv` | Cleaned BIEN occurrence coordinates | **NOT FOUND** at the archived path claimed (`p_phyloGWAS_archived/output/envData/speciesRange/...`); instead verified present at `/local/workdir/sh2246/p_phyloGWAS/output/envData/speciesRange/bien_coordinates_clean_2023.12.04.csv` — i.e. it is already in this repo's own `output/`, not in the archived sibling | Derived from BIEN (`BIEN::BIEN_occurrence_species`) | 06D | |
| `asreml_4.2.0.302_linux-intel64_R4.2.0.tar.gz` | ASReml-R v4.2.0 installer archive | `p_phyloGWAS_archived/asreml_4.2.0.302_linux-intel64_R4.2.0.tar.gz` (verified) | VSNi (licensed software, see Section 9) | 08C, 08D | |

**Correction to the verified inventory:** the `bien_coordinates_clean_2023.12.04.csv` path
given in the original inventory (under `p_phyloGWAS_archived/output/...`) does not exist.
The file was located instead already inside this repo, under
`output/envData/speciesRange/`. No action needed — it is not actually a missing
dependency, just mis-attributed to the wrong repo in the original note.

---

## 4. `p_evolBNI/` GIS rasters (present, undeclared) [T1]

Sibling repo `p_evolBNI`, directory `p_evolBNI/data/GIS_env_data/`. Used by envirotyping
stages 06B/06D via `src/09_pulling_envData.r` (and variants). Not under this repo's
version control.

| File/dataset | What it is | Current location | Canonical source | Used by (stage) | Reviewer action |
|---|---|---|---|---|---|
| `wc2.1_2.5m_elev.tif` | WorldClim 2.1 elevation raster (2.5 arc-min) | `p_evolBNI/data/GIS_env_data/WorldClim_raw_2.5m_files/wc2.1_2.5m_elev/wc2.1_2.5m_elev.tif` (verified) | WorldClim 2.1 | 06B/06D | |
| `wc2.1_2.5m_bio_1.tif` (+ other bioclim variables) | WorldClim 2.1 bioclimatic variable raster(s) | `p_evolBNI/data/GIS_env_data/WorldClim_raw_2.5m_files/wc2.1_2.5m_bio/wc2.1_2.5m_bio_1.tif` (verified) | WorldClim 2.1 | 06B/06D | |
| `HYSOGs250m.tif` | Global Hydrologic Soil Group raster (250 m) | `p_evolBNI/data/GIS_env_data/Global_Hydrologic_Soil_Group_1566/Global_Hydrologic_Soil_Group_1566/data/HYSOGs250m.tif` (verified; note the doubled directory nesting) | NASA / ORNL DAAC Global Hydrologic Soil Groups (1566) | 06B/06D | |
| `GSDE_raw_nc_files/*.nc` | Global Soil Dataset for Earth system modeling (per-variable NetCDF, e.g. `BD1/BD1.nc`, `BS1/BS1.nc`, `TK1/TK1.nc`, ...) | `p_evolBNI/data/GIS_env_data/GSDE_raw_nc_files/` (verified: 33 `.nc` files present across per-variable subdirectories) | GSDE (Global Soil Dataset for Earth system modeling) | 06B/06D | |

---

## 5. Off-machine siblings and legacy path references [T1]

| Item | Status | Detail | Reviewer action |
|---|---|---|---|
| `p_panAndOGASR/` | **Confirmed absent on this machine** (`test -d` fails) | Referenced only from `notebook/07_summaryStats/07Ab_prematureStopCodon.ipynb`, which is a **legacy/alternate** stop-codon analysis notebook. Reads `data/sixDigitCode2AnnotID`, `output/phyloP_res_filtered_v2.txt` / `output/phyloP_og.filtered.txt`, `output/OG_4-8annual17perennial.csv`, and others, all under the missing `p_panAndOGASR/`. The **active** premature-stop-codon path is 07Ba/07Bb (run on SCINET, Section 8), not 07Ab. Flag as needs-attention (dead legacy path), not a pipeline blocker. | |
| `p_geneModelEvaluation/` | **Confirmed absent on this machine** (`test -d` fails) | `src/S00_basicFunction.R` is sourced from `/workdir/sh2246/p_geneModelEvaluation/src/S00_basicFunction.R` in several notebooks: `notebook/09_molEvolution/09B_RELAX_resultSummary.ipynb` (active), `notebook/08_linearModeling/archived/08B_perOGModeling.ipynb` (archived), `notebook/10_aprioriCandidate/archived/10B_stressInducedGene_enrichment.ipynb` and `archived/10C_nutrientCyclingGene_enrichment.ipynb` (archived). **However**, the file has already been vendored into this repo at `src/S00_basicFunction.R` (verified present), and the active `notebook/10_aprioriCandidate/10B_candidateOGInvestigation.ipynb` and `notebook/10_aprioriCandidate/archived/10E_candidateOGInvestigation_lifeHistory.ipynb` already point at the local `/workdir/sh2246/p_phyloGWAS/src/S00_basicFunction.R` copy. Remaining work: update the still-outstanding reference in the active `09B_RELAX_resultSummary.ipynb` to point at the local vendored copy instead of the missing sibling repo. | |
| `Zm00001eb.1.fulldata_goList.txt`, `Pv_GOTable.txt` | No action needed | Legacy notebooks (e.g. `07Ab_prematureStopCodon.ipynb`) reference these under the missing `p_panAndOGASR/data/`, but both files are already mirrored into `data/` in this repo (Section 2) — the **active** pipeline is unaffected since it reads the local copies. | |

---

## 6. External databases / live queries [T1]

| Item | What it is | Access pattern | Canonical source | Used by (stage) | Reviewer action |
|---|---|---|---|---|---|
| GBIF occurrence records | Species occurrence points | Live query via `rgbif::occ_search` | GBIF.org | 06B (`src/08_pulling_geo_data.R`) | |
| BIEN occurrence records | Species occurrence points | Live query via `BIEN::BIEN_occurrence_species` | Botanical Information and Ecology Network (BIEN) | 06B (`src/08_pulling_geo_data.R`) | |
| Derived occurrence dataset (citable) | Cleaned/combined GBIF+BIEN occurrence dataset used in the manuscript | N/A (deposited dataset) | **Zenodo 14968186** (per manuscript) | 06B/06D | |
| `get_spatial_fun.R` | Helper function for spatial/environmental extraction, fetched at runtime | Live fetch via `source()`/download at runtime | `https://raw.githubusercontent.com/gcostaneto/envirotypeR/main/R/get_spatial_fun.R` | `src/09_pulling_envData.r` | |

**Discrepancy to reconcile:** `src/08_pulling_geo_data.R` (line ~91) contains a code
comment citing `10.5281/zenodo.14967966`, which does not match the manuscript's stated
Zenodo DOI of **14968186**. These are two different, very close Zenodo record IDs — worth
confirming which is correct and updating the stale one (likely the in-code comment).

---

## 7. Model weights [T1]

| Model | What it is | Access pattern | Canonical source | Used by (stage) | Reviewer action |
|---|---|---|---|---|---|
| ESM2 (`facebook/esm2_t33_650M_UR50D`) | Protein language model, used for variant-effect logits | Downloaded from HuggingFace Hub at runtime (GPU) | HuggingFace (`facebook/esm2_t33_650M_UR50D`) | 07C (`src/4_ESM_logits.py`) | |
| PlantCAD (`kuleshov-group/PlantCaduceus_l32`) | Plant genomic language model | Not loaded directly in the traced active scripts — PlantCAD scores enter modeling only as a precomputed predictor column | HuggingFace (`kuleshov-group/PlantCaduceus_l32`); precomputed scores in `combined_plantCAD_zeroShotScores_20250421.txt` | Modeling stage (predictor input only) | |

---

## 8. SCINET-only raw inputs [T1]

Live only on USDA SCINET Atlas, SLURM account `buckler_lab_panand`. Not reachable/verifiable
from this machine.

| File/dataset | What it is | Current location | Canonical source | Used by (stage) | Reviewer action |
|---|---|---|---|---|---|
| `orthofinderProteinMSAs/` | OrthoFinder protein multiple-sequence alignments, by orthogroup | `/project/90daydata/buckler_lab_panand/aimee.schulz/panand/output/orthofinderProteinMSAs/` (SCINET Atlas; not accessible/verifiable from this machine) | Generated by A. Schulz (OrthoFinder run) | 07Bb → `src/11_find_premature_stops.py` (premature stop codons); 07Cb → `src/4_ESM_logits.py` (ESM scoring) | |
| `orthofinderProteinMSAs_additionalOGs/` | Same as above, for a supplementary set of orthogroups | `/project/90daydata/buckler_lab_panand/aimee.schulz/panand/output/orthofinderProteinMSAs_additionalOGs/` (SCINET Atlas; not accessible/verifiable from this machine) | Generated by A. Schulz (OrthoFinder run) | 07Bb, 07Cb | |

---

## 9. Licensed / external tools [T1]

| Tool | Version | Current location | Notes | Used by (stage) | Reviewer action |
|---|---|---|---|---|---|
| ASReml-R | 4.2.0.302 (linux-intel64, R4.2.0) | Installer verified present at `p_phyloGWAS_archived/asreml_4.2.0.302_linux-intel64_R4.2.0.tar.gz` (Section 3) | **Licensed software (VSNi)** — must be separately licensed/installed by the reproducing party | `src/12_runPermulation_perOGModel.R` (08C), `src/powerSimulation_XY_revised.R` (08D) | |
| HyPhy | 2.5.49 | Verified present at `/programs/hyphy-2.5.49/bin/hyphy` | Assumed pre-installed on the compute environment | Molecular evolution (selection tests) | |
| RAxML | (standard-RAxML build) | Verified present at `src/standard-RAxML/raxmlHPC` (vendored into this repo) | — | Phylogenetic tree inference | |
| seqkit | 0.15.0 | Not verified on this machine (assumed installed per environment spec) | Assumed pre-installed | Sequence manipulation utilities | |
| miniprot | — | Not verified on this machine (assumed installed per environment spec) | Assumed pre-installed | Protein-to-genome alignment (03C) | |
| gffread | — | Verified present at `/programs/bin/cufflinks/gffread` | — | Annotation/GFF processing | |
| GNU parallel | — | Verified present at `/programs/bin/parallel` | — | Pipeline parallelization across many scripts | |

---

## 10. "Read from `data/` but not in `data/`" pointers [T1/T2]

These notebooks reference `data/<file>`, but the referenced file is currently only
present in the sibling repo `p_phyloGWAS_archived/data/` (all verified present there).
To reproduce these stages, copy the listed file from the given archived path into this
repo's `data/` directory.

| Referenced as | Used by (stage) | Verified source path to copy from | Reviewer action |
|---|---|---|---|
| `data/spNameMetadata_20240819.txt` | 06B | `p_phyloGWAS_archived/data/spNameMetadata_20240819.txt` | |
| `data/allOGs_withHeaders_dnds.txt` | 06C | `p_phyloGWAS_archived/data/allOGs_withHeaders_dnds.txt` | |
| `data/Perennitality_Numeric_Poaceae.csv` | 07Ab | `p_phyloGWAS_archived/data/Perennitality_Numeric_Poaceae.csv` | |
| `data/rhizomePAV_traits.txt` | 07Ab | `p_phyloGWAS_archived/data/rhizomePAV_traits.txt` | |
| `data/Angiosperms353_orysaSequences.fasta` | 05A | `p_phyloGWAS_archived/data/Angiosperms353_orysaSequences.fasta` | |
| `data/Zea_mays_v5_mrna.fa` | 05A | `p_phyloGWAS_archived/data/Zea_mays_v5_mrna.fa` | |
| `data/Araport11_seq_20220914_representative_gene_model.fa.gz` | 05A | `p_phyloGWAS_archived/data/Araport11_seq_20220914_representative_gene_model.fa.gz` | |
| `data/Yeaman2024_RAO_AtID.txt` | 05A/10 | `p_phyloGWAS_archived/data/Yeaman2024_RAO_AtID.txt` | |

---

## Summary of path issues found

- **NOT FOUND (mis-attributed, not missing):** `p_phyloGWAS_archived/output/envData/speciesRange/bien_coordinates_clean_2023.12.04.csv` does not exist; the file was found instead at `p_phyloGWAS/output/envData/speciesRange/bien_coordinates_clean_2023.12.04.csv` (i.e. already inside this repo). See Section 3.
- **Needs-attention (dead legacy path, not a blocker):** `07Ab_prematureStopCodon.ipynb` depends on the absent `p_panAndOGASR/` sibling; the active premature-stop-codon path is 07Ba/07Bb on SCINET. See Section 5.
- **Needs-attention (one file to vendor):** `09B_RELAX_resultSummary.ipynb` still sources `S00_basicFunction.R` from the absent `p_geneModelEvaluation/` sibling, even though the file is already vendored locally at `src/S00_basicFunction.R` and other notebooks (10B, 10E) already use the local copy. See Section 5.
- **Discrepancy to reconcile:** Zenodo DOI in `src/08_pulling_geo_data.R` code comment (`10.5281/zenodo.14967966`) vs. manuscript-cited DOI (`14968186`). See Section 6.
- **Not verifiable from this machine (expected):** SCINET Atlas paths in Section 8; `seqkit`/`miniprot` binaries in Section 9 (assumed installed per environment spec, not checked here).
