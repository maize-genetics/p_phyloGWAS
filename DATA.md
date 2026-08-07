# DATA.md — Input-Data Manifest

Every input needed to reproduce the phylogenetic-GWAS pipeline, reviewed and reduced to
what's actually necessary (unneeded/superseded items removed — see git history for the
full review trail). Four columns only: what it is, where it lives, which stage reads it,
and whether it's needed for a full raw-to-figures run or only to resume from provided
intermediates.

**Tier:** **T1** = needed to run the full pipeline end-to-end from raw genomes. **T2** =
an intermediate that a script in this repo can regenerate — provided as a shortcut so you
don't have to rerun the upstream stage yourself.

Repo root: `/local/workdir/sh2246/p_phyloGWAS`. Sibling repo: `/local/workdir/sh2246/p_evolBNI`.
Raw sequence/assembly release (row 1) is under a publication-timed embargo per the
manuscript's Data Availability Statement — expected, not a repo bug.

## What stays on GitHub

Every file listed here is small enough (~8.7 MB combined) to just live in this repo directly —
`.gitignore` carves out explicit exceptions for exactly these, on top of its blanket `data/`/
`output/` ignore: `data/Poaceae_metadata_2024.08.21.csv`, `data/Poaceae_metadata_2025.04.03_CH.csv`,
`data/Poaceae_metadata_2025.04.07.csv`, `data/Poaceae_metadata_filtered_2025.04.07.csv`,
`data/Poaceae_metadata_filtered_2025.08.28.tsv`, `data/Poaceae_metadata_highErrorFiltered_2025.10.08.tsv`,
`data/fullQCData_20240819.txt`, `data/spNameMetadata_20240819.txt`, `data/env_metadata.txt`,
`data/DEG_study_metadata.csv`, `data/OGToZm_mapping_v2.txt`, `data/poaceaeRepAssemblies_32.txt`,
`data/poaceaeHelixerOG_ancSeq_gapRemoved_v2_uppercase_output_CC_f.gmt`,
`data/busco_gene_OG_20250828.txt`, `data/maize_busco.tsv`,
`output/poaceaeHelixerOG_filtered_20250331.txt`.

Everything else in this manifest is one of:
- **Large T1 raw/reference data, deposited externally** — marked "external, link TBD" below
  until the author provides the actual link/DOI. Destination by category: genome assemblies →
  agDataCommons; GIS raster data → Zenodo; `data/Zmays_cds.fa`, `data/Zea_mays_v5_mrna.fa`,
  `data/maxFPKM_v5.csv`, and `data/Zm00001eb.1.fulldata_goList.txt` → MaizeGDB (not planned to
  be hosted anywhere by this project at all, just cite the source — resolved, no link needed);
  `data/Angiosperms353_orysaSequences.fasta` → also resolved, cites the source GitHub repo/paper
  directly in its own row rather than needing a placeholder; ESM2/PlantCAD → HuggingFace
  (already resolved, pulled at runtime, see "Model weights" below); GBIF/BIEN occurrence
  records → live API query, or the existing Zenodo `14968186` shortcut.
  `data/DEG_CDS_fasta/` is a separate case — re-derived/re-hosted copies of already-public
  genome data rather than a single published dataset with one link; see its own row for the
  per-file provenance ledger and the open "whether to deposit at all" question.
- **T2 intermediates** — not hosted anywhere; regenerate by rerunning the stage/script named
  in each row's Tier column. (The handful of small T2 files listed above are the exception —
  cheap enough to just keep in the repo as a convenience even though they're regenerable.)

---

## Genomes & reference sequences

| Data | Path | Used by | Tier |
|---|---|---|---|
| 727 raw genome assemblies (504 GB) | `data/assemblies/` (external — agDataCommons, link TBD) | all stages | T1 |
| 32 representative genome assemblies + Helixer gene annotations | genomes are a named subset of the row-1 `data/assemblies/` set — see `data/poaceaeRepAssemblies_32.txt` (tracked on GitHub). One substitution from the original OrthoFinder run: `Pp-Kellogg1297-DRAFT-PanAnd-1.0` (draft, used to build the OG queries) → `Poa_pratensis_v1` (the same species' genome, published since) — reproduce with the published genome. Helixer annotations are **not deposited or otherwise hosted by this project** — regenerate them yourself by running [Helixer](https://github.com/weberlab-hhu/Helixer) on these 32 assemblies. | 03 (OrthoFinder input) | T1 |
| Maize reference CDS | external — MaizeGDB (`data/Zmays_cds.fa` locally; not tracked in this repo, not planned for any external deposit by this project — just cite MaizeGDB as the source) | orthology, molecular evolution | T1 |
| Angiosperms353 Oryza reference sequences | `data/Angiosperms353_orysaSequences.fasta` — a subset of just the `Orysa` (*Oryza sativa*) entries pulled out of the full Angiosperms353 universal target file ([mossmatters/Angiosperms353](https://github.com/mossmatters/Angiosperms353), [Johnson et al. 2019, *Syst Biol* 68(4):594–606](https://academic.oup.com/sysbio/article/68/4/594/5237557)), which covers 353 genes across angiosperms via 5–15 exemplar sequences per gene (k-medoids clustering, chosen to keep 95% of angiosperm sequences within 30% divergence of some target instance). Restricting to a single grass exemplar for grass-only cross-mapping is reasonable — it should sit closer, on average, to any other Poaceae genome than the file's ~300 non-grass exemplars, which exist to cover lineages irrelevant here. **Known coverage gap**: `Orysa` covers 332 of the 353 genes (94%); the 21 missing genes have *no* grass-lineage exemplar anywhere in the source file at all — `Sorbi` (sorghum, the file's only other grass genome) is also absent from those same 21 genes, so this is very likely a real absence of any monocot/grass homolog identified for those genes during the original probe design, not an artifact of picking rice specifically. | 03C (miniprot cross-mapping) | T1 |
| Maize B73 v5 mRNA reference | external — MaizeGDB (`data/Zea_mays_v5_mrna.fa` locally; not tracked in this repo, not planned for any external deposit by this project — just cite MaizeGDB as the source) | 03C (miniprot cross-mapping) | T1 |
| Per-species CDS FASTAs (a-priori stress genes) | `data/DEG_CDS_fasta/` — per-file provenance (source study, species, and the exact `CDS.fa` filename) is documented row-by-row in `data/DEG_study_metadata.csv`'s `study`/`species`/`CDS.fa` columns. Each file was either downloaded as a pre-generated CDS FASTA from its source, or self-extracted via `gffread` from that source's own genome+GFF — not novel data generated by this project. **Deposit status undecided**: since these are re-derived/re-hosted copies of already-public reference genome data rather than original data, whether (and where) to redistribute them separately is still an open question, not yet a "link TBD" with a settled destination. | 10A | T1 |
| DEG study metadata (manually-extracted DE gene lists per stress study, with paths to each list under `output/candidateGenes/`) | `data/DEG_study_metadata.csv` (tracked on GitHub) | 10B | T1 |

## Metadata & trait data

| Data | Path | Used by | Tier |
|---|---|---|---|
| Poaceae accession metadata (initial manual curation) | `data/Poaceae_metadata_2024.08.21.csv` (tracked on GitHub) | 02A | T1 |
| Poaceae accession metadata (further-filtering inputs/outputs) | `data/Poaceae_metadata_2025.04.03_CH.csv` (input), `data/Poaceae_metadata_2025.04.07.csv`, `data/Poaceae_metadata_filtered_2025.04.07.csv` (outputs) — all tracked on GitHub | 02B, read downstream by 03B/05B | T2 (regenerated by 02B) |
| Poaceae accession metadata (post-tree-filtering / QC-filtering) | `data/Poaceae_metadata_filtered_2025.08.28.tsv`, `data/Poaceae_metadata_highErrorFiltered_2025.10.08.tsv` — both tracked on GitHub | outputs of 08A; read by 05B, 06B, 08A/08B/08C, 03D, `src/S05_envPC_analysis.R` | T2 (regenerated by 08A) |
| Per-assembly QC statistics | `data/fullQCData_20240819.txt` (tracked on GitHub) | 02A | T1 |
| Species-name list | `data/spNameMetadata_20240819.txt` (tracked on GitHub) | 06A | T2 (regenerated by 02A) |
| Environmental metadata | `data/env_metadata.txt` (tracked on GitHub) | stage 06 | T1 |

## Orthogroup & gene-mapping data

| Data | Path | Used by | Tier |
|---|---|---|---|
| dN/dS table (full orthogroup set) | `data/fullSetOGs_240903.txt` | output of stage 07 (tip-to-outgroup dN/dS calculation); used as a predictor in 08's phylogenetic mixed model | T2 (regenerated by stage 07) |
| Rice → orthogroup mapping | `data/DEG_mappingFiles/` | 03E, 10A, 10B, 11 | T2 (regenerated by 10A) |
| Orthogroup → maize gene mapping | `data/OGToZm_mapping_v2.txt` (tracked on GitHub) | 03E, 08A, 08B, 11 | T2 (regenerated by 03C) |
| Per-OG miniprot GFF annotations | `output/orthofinderMiniProt/`, `output/orthofinderMiniProt_additionalOGs/` | 03D, 07Aa | T2 (regenerated via miniprot) |
| MPID → OG:assembly:index mapping (2 independent miniprot batches - MPID not comparable across them) | `output/seqIDmapping.txt`, `output/seqIDmapping_additionalOGs.txt` | 07Aa | T2 (regenerated via miniprot; migrated from `p_phyloGWAS_archived` — not part of the earlier consolidation pass) |
| DeepGO GO annotation | `data/poaceaeHelixerOG_ancSeq_gapRemoved_v2_uppercase_output_CC_f.gmt` (tracked on GitHub) | 11 | T2 (regenerated via DeepGO) |
| Maize B73 v5 GO annotation | external — MaizeGDB (`data/Zm00001eb.1.fulldata_goList.txt` locally; not tracked in this repo, not planned for any external deposit by this project — just cite MaizeGDB as the source) | 08A, 11 | T1 |
| Maize v5 max-expression (FPKM) | external — MaizeGDB (`data/maxFPKM_v5.csv` locally; not tracked in this repo, not planned for any external deposit by this project — just cite MaizeGDB as the source) | expression-informed filtering | T1 |
| OrthoFinder protein MSAs | `output/orthofinderProteinMSAs_fullset_20250710/` (current, full OG set), `output/orthofinderProteinMSAs/` (older, pre-consolidation) | 07Ba, 07Cb, 08B, 11B | T2 (regenerated via OrthoFinder) |
| Per-OG MSA (raw, pre-gap-stripping) | `output/OrthofinderMAFFT/` — on this machine a renamed copy currently sits at `output/CDSMSAPerOG/` (see stage-04 README's note) | 05A (gap-stripping), 08B (GC composition) | T2 (regenerated by stage 04, `notebook/04_msaGeneration/README.md` Step 0 + Step 1) |
| Per-taxon (whole-genome) protein and CDS FASTAs (amino-acid/GC composition input) | `output/aminoAcidPerTaxa/`, `output/CDSperTaxa/` | 08B | T2 (regenerated via `gffread` from stage 03's miniprot GFFs, see stage-04 README's Step 2) |
| Filtered OG list (stage-03 OG universe, used as a null-sampling background) | `output/poaceaeHelixerOG_filtered_20250331.txt` (tracked on GitHub) | 10B | T2 (regenerated by 03B) |
| BUSCO gene ↔ OG mapping (QC) | `data/busco_gene_OG_20250828.txt` (tracked on GitHub) | QC / orthology cross-check | T2 |
| Maize BUSCO completeness (QC) | `data/maize_busco.tsv` (tracked on GitHub) | QC | T2 |

## Environmental raster data

| Data | Path | Used by | Tier |
|---|---|---|---|
| WorldClim 2.1 elevation | `p_evolBNI/data/GIS_env_data/WorldClim_raw_2.5m_files/wc2.1_2.5m_elev/` (external — Zenodo, link TBD) | 06A (`src/09_pulling_envData.r`, via `GIS_DATA_ROOT`) | T1 |
| WorldClim 2.1 bioclim variables | `p_evolBNI/data/GIS_env_data/WorldClim_raw_2.5m_files/wc2.1_2.5m_bio/` (external — Zenodo, link TBD) | 06B (envPC-pipeline map figure, via `GIS_DATA_ROOT`) | T1 |
| Global Hydrologic Soil Groups | `p_evolBNI/data/GIS_env_data/Global_Hydrologic_Soil_Group_1566/.../HYSOGs250m.tif` (external — Zenodo, link TBD) | 06A (`src/09_pulling_envData.r`, via `GIS_DATA_ROOT`) | T1 |
| GSDE soil dataset (33 NetCDF files) | `p_evolBNI/data/GIS_env_data/GSDE_raw_nc_files/` (external — Zenodo, link TBD) | 06A (`src/09_pulling_envData.r`, via `GIS_DATA_ROOT`) | T1 |

## Occurrence data (external DB)

| Data | Path | Used by | Tier |
|---|---|---|---|
| GBIF occurrence records | live query (`rgbif::occ_search`) → `output/metadataFormalOut/` | 06A | T1 |
| BIEN occurrence records | live query (`BIEN::BIEN_occurrence_species`) → `output/metadataFormalOut/` | 06A | T1 |
| Derived occurrence dataset (cleaned GBIF+BIEN) | Zenodo `14968186` | 06A, 06B | T2 (shortcut for the two rows above) |

## Model weights

| Data | Path | Used by | Tier |
|---|---|---|---|
| ESM2 protein language model | HuggingFace `facebook/esm2_t33_650M_UR50D` | 07 (zero-shot scoring) | T1 |
| PlantCAD nucleotide language model | HuggingFace `kuleshov-group/PlantCaduceus_l32` | 07 (zero-shot scoring) | T1 |

Both model weights are pulled directly from HuggingFace at runtime — no separate file to
gather. (Software tools — RAxML, HyPhy, miniprot, ASReml-R, etc. — are tracked in
`WORKFLOW.md`'s "Tools used across stages" table, not here; this manifest is data only.)
