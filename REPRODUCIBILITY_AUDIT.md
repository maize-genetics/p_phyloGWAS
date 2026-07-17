# p_phyloGWAS — Repository & Reproducibility Audit

**Scope:** read-only review of `notebook/`, `src/`, `data/`, `output/`, `log/`, and `README.md`. No files were edited or deleted.
**Manuscript:** "Convergent genome- and gene-level constraints shape repeated environmental adaptation in grasses" (Hsu, Schulz, Hale, et al.; Cornell IGD/USDA-ARS; bioRxiv 10.64898/2026.06.01.729361). **Update:** the full manuscript text (Cell-formatted PDF, including Methods, all main figures, and the Data Availability Statement) was supplied directly and has now been read in full — the earlier version of this audit only had the abstract (bioRxiv blocked automated fetching with a 403). Section 3 below is rewritten against the full Methods/Results text and resolves several previously-open questions. Per the author's note, **figure numbers/ordering in the manuscript have shifted since the code was last run** — several `output/figure/*.png` filenames (e.g. `Fig8_revised.png`) reflect an older numbering and no longer match the current Figure 1–6 layout; see §3.1 for the mapping.

---

## 1. What this repository is

A single Cornell (Buckler/Hsu) lab working directory serving as both the analysis codebase and the (partial) data/output store for a cross-species phylogenetic-GWAS study of environmental adaptation in Poaceae (707 genomes / 569 species). It is **not** a packaged, portable pipeline — it's a numbered sequence of notebooks/scripts (`02_...` through `10_...`) that were run interactively and iteratively on Cornell BioHPC (`cbsublfs1`, `/workdir/sh2246/...`) and USDA SCINET (`buckler_lab_panand` SLURM account), with results accumulated in `output/`.

Git only tracks `notebook/` (121 files), `src/` (254 files), `README.md`, and 6 files in `log/`. **`data/`, `output/`, and (mostly) `log/` are gitignored** — none of the actual genomic data or result tables are version-controlled. `src/interproscan-5.76-107.0/` is also gitignored (a vendored tool install). On disk, `data/` is 511 GB, `output/` is 1.2 TB, `src/` is 49 GB (dominated by vendored binaries/tool installs, see §5).

## 2. Pipeline stage map

Numbering follows `notebook/`'s own scheme. "Stage 00/01/04" notebooks do not exist (see §4.4).

| Stage | Directory | Purpose | Status |
|---|---|---|---|
| 02 | `02_metadataCuration` | Merge PanAnd/LIMS QC metadata + manual species-ID curation → filtered genome list (808→728 genomes) | Runs, but depends on a file from `p_phyloGWAS_archived` (see §4.2) |
| 03 | `03_orthogroup` | OrthoFinder (34 reps) → orthogroups → ancestral AA reconstruction → miniprot cross-mapping to Zm/At/Pv/Angiosperm353 → OG filtering (99,140→22,363 OGs) | Runs; several inputs sourced from `p_phyloGWAS_archived` |
| 05 | `05_phylotreeConstruction` | Gap-strip CDS MSAs → RAxML gene trees (all-OG + angiosperm353) → ASTRAL-pro species tree → phylogenetic K matrix | `05A_treeConstruction` is a plain shell script (no extension — confirmed by reading it, contra its file-manager appearance as an empty folder) with fully hardcoded `/workdir/sh2246/p_phyloGWAS/` paths |
| 06 | `06_envirotyping` | Species coordinates (BIEN/GBIF) → WorldClim environmental rasters → envPC1-3 → Fig1/env figures | `06D` supplementary-figure script depends on `p_phyloGWAS_archived` **and** the sibling repo `p_evolBNI` for raw WorldClim tif files |
| 07 | `07_summaryStats` | Premature stop codons (PMS), frameshifts, ESM2 & PlantCAD zero-shot scores per OG | `07B*`/`07C*` are SCINET-only SLURM pipelines (A100 GPU partition, `buckler_lab_panand` account, `~/miniconda3` user-specific env activation) |
| 08 | `08_linearModeling` | Master data table (genomic features × envPC × phyloK) → genome-wide feature association (MLM, → Fig3) → per-OG ASReml phylogenetic mixed models + permulation (→ Fig5) → power simulation (→ Fig2, see §3.1 on renumbering) | Requires licensed **ASReml-R v4.2.0**; `08B` pulls a mapping file from `p_phyloGWAS_archived` |
| 09 | `09_molEvolution` | MSA cleaning → RAxML gene trees → foreground/background node labeling → HyPhy RELAX per trait (temp/precip/soil/life-history/rhizome) → result tables | Depends on cluster-installed `/programs/hyphy-2.5.49/`; `09B` GO enrichment pulls a GO-term file from **`p_panAndOGASR`** (a third, different sibling repo) |
| 10 | `10_aprioriCandidate` | miniprot OG→gene-ID mapping (Helixer annotations) → integrate ASReml + RELAX + expression evidence → candidate OG lists → Fig6 Sankey of the filtering cascade | Same `p_panAndOGASR`/`p_phyloGWAS_archived` dependency as above |
| `slurm/`, `XX_archived/`, `*/archived/` | — | SLURM job templates; superseded/exploratory notebooks (scaffolding tests, haploidization QC, PAML, hackathon scripts, etc.) | Not part of the current pipeline; safe to ignore for reproducing the manuscript, but currently intermixed with active code with no clear "deprecated" signal besides folder name |

Vendored vs. authored code: `src/ASTER` (ASTRAL-pro, embedded as a **nested git repo with its own `.git/`**, not a registered submodule — this is why `git status` shows it as modified even though nothing in it changed), `src/standard-RAxML` (RAxML source, includes compiled `.o` files committed to git), `src/patch-scaffolds` (a compiled Kotlin/JVM tool with ~90 jars), and `src/interproscan-5.76-107.0` (gitignored) are all third-party tool installs living inside `src/`, not analysis code.

**Tool/package versions pinned in the Methods** (useful as the seed list for the `environment.yml`/`renv.lock` recommended in §7 — none of these are currently pinned anywhere in the repo itself):

| Tool/package | Version | Used for |
|---|---|---|
| Helixer | (Stiehler et al. 2021) | De novo annotation of the 32 representative genomes |
| OrthoFinder | v2.6.4 | Orthogroup construction |
| R/phangorn | v2.12.1 | Ancestral OG protein sequence reconstruction |
| miniProt | v0.13.0 | Querying ancestral OG sequences against all 727 genomes; angiosperm353/Zm/Pv/At cross-mapping |
| MAFFT | v7.520 | Nucleotide + amino acid MSA per orthogroup |
| RAxML | v8 (GAMMA+GTR) | Gene trees (all-OG and angiosperm353 subsets) |
| ASTRAL-Pro | v2 | Species tree reconciliation from gene trees |
| R/ape | v5.8-1 | Ancestral envPC state reconstruction / transition counting |
| R/canprot | v2.0.0 | Amino-acid physiochemical properties (N/C ratio, energetic cost, hydropathy, MW, density) |
| SeqKit | v2.13.0 | GC content of CDS |
| R/MSA2dist | v1.2.0 | Tip-to-outgroup dN/dS |
| PlantCAD (PlantCaduceus_l32) | kuleshov-group | Nucleotide LLM zero-shot score |
| ESM2 | facebook/esm2_t33_650M_UR50D | Protein LLM zero-shot score |
| ASReml-R | **v4.2.0** | All phylogenetic mixed models (licensed — see §4.5) |
| R/phylolm | v2.6.5 | Independent cross-check of ASReml-R results (unlicensed alternative; rho > 0.998 agreement reported) |
| R/asremlPlus | (unversioned in text) | Full-vs-reduced model likelihood ratio tests |
| Phylogenetic permulation | per Saputra et al. 2021 | Empirical p-values for all mixed-model and RELAX tests |
| HyPhy RELAX | (72) | Selection-intensity shift testing on labeled gene trees |
| R/topGO | v2.50.0 | GO enrichment ("elim"/"classic" algorithms) on maize v5 GO annotation |
| deepGO | — | De novo GO annotation of ancestral OG sequences (independent of maize-homolog-based topGO run) |
| targetP | 2.0 | Signal/transit peptide prediction for high-confidence candidates |
| InterProScan | (vendored in `src/interproscan-5.76-107.0`, gitignored) | Protein domain annotation |
| AlphaFold3 | — | Structure prediction for candidate mutational effects |

## 3. Cross-check against the full manuscript

The full manuscript text (not just the abstract) is now available and confirms most of the internal `output/` evidence gathered in the first pass, resolves the previously-unverified "330 genes" figure, and explains several filename/version numbers that looked arbitrary from the code alone.

### 3.1 — Figure/table renumbering (per the author's note)
The manuscript's current figure order (Fig. 1 phylogeny/climate zones & transition counts; Fig. 2 phylogenetic-mixed-model **power simulation**; Fig. 3 genome-wide molecular/physiochemical trait associations; Fig. 4 LLM score validation; Fig. 5 per-OG phylogenetic mixed model, envPC1-3; Fig. 6 Sankey + 17 high-confidence candidate table) **does not match** several filenames still sitting in `output/figure/`:
- The power-simulation figure is now **Figure 2**, but the repo file is named `Fig8_revised.png` (generated by `src/powerSimulation_XY_revised.R` / `08D_power_simulation.sh`) — a leftover from an earlier draft numbering where the power analysis was placed near the end.
- `08A_masterDataTableGeneration.ipynb` writes `Fig2.png`, `Fig3c_v2.png`, `Fig3d.png` — the genome-wide trait-association content now lives in **Figure 3**, so `Fig3c_v2.png`/`Fig3d.png` are plausibly still correct, but `Fig2.png` is stale (Figure 2 is now the power simulation, produced elsewhere).
- The manuscript's Results/Discussion only walk through Figures 1–6; there is no Figure 7 in the current text, and Figure 8 no longer exists as a main figure (its content — the power simulation — was pulled forward to Figure 2). `output/figure/Fig8.png` / `Fig8_revised.png` / `FigS6*.pdf` are therefore artifacts of a prior draft's numbering, not missing current figures.
- **Net effect for reproducibility**: the code that generates each figure's *content* still appears correct and traceable, but a reader trying to match `output/figure/*.png` files 1:1 against the current manuscript figure numbers will be misled by at least the power-simulation file. Recommend renaming figure output files to match the final manuscript numbering (or better, deriving the figure number from a config rather than a literal filename) as part of the cleanup pass in §7.

### 3.2 — Numbers that now resolve cleanly
- **727 assemblies / 589 species → 707 genomes / 569 species**: Methods states 727 total assemblies (217 public + 33 new long-read + 487 short-read) across 589 species; 20 species lacked habitat/occurrence data and were dropped, leaving the abstract's 707 genomes / 569 species used for envPC modeling. This lines up with `output/envData_707Poaceae_*.txt` and confirms that number's provenance (previously just observed as "the file is called 707," not derived).
- **707 → 687 genomes for OG/gene-activity testing**: Methods further reports 20 *additional* assemblies were dropped for the gene-activity analysis specifically (>15% frameshift/premature-stop rate, attributed to sequencing/assembly error), leaving 687 genomes for the 19,613-OG gene-activity tests. This almost certainly corresponds to `data/Poaceae_metadata_highErrorFiltered_2025.10.08.tsv`, which the first audit pass flagged as having an "unclear versioning strategy" — it's now explained as the QC-error-rate exclusion step, not ad hoc manual curation.
- **22,363 filtered OGs, "present in ≥8 taxa"**: matches `output/poaceaeHelixerOG_filtered_20250331.txt` and `03B_OGFilter.ipynb` exactly, including the ≥8-taxa threshold (Methods doesn't mention the "≤100 copies/taxon" upper bound the notebook also applies — worth a comment in the notebook clarifying that the extra filter is a QC step beyond what's stated in Methods).
- **19,613 "testable" OGs**: Methods defines these as present in the outgroup *Pharus latifolius* AND in ≥100 of the 687 genomes — this matches the `output/geneTree_allOGs_*` / RELAX pipeline's OG count referenced by the second audit pass.
- **32 vs. "34" representative genomes**: Methods states 32 representative long-read assemblies were used for OG construction, but `notebook/03_orthogroup/03A_buildHelixerOG.sh` has a comment reading "Charlie's 34 representative assemblies." This is a minor but real drift between the code comment and the published number — likely two long-read assemblies were dropped between when that script was written and the final run; worth a one-line note in the script reconciling the discrepancy.
- **envPC1 (27.0%)/envPC2 (17.8%)/envPC3 variance explained, from 282 environmental variables**: matches the envPC pipeline (`06B_spCoordEnvData.sh` → `src/08_pulling_geo_data.R`/`09_pulling_envData*.R`) and `output/envPC_load.txt`/`envPC_loading_20250804.png`.
- **330 candidate genes (previously "not independently verified")**: **now resolved.** Figure 6A's Sankey shows 160 (envPC1) + 109 (envPC2) + 68 (envPC3) candidate OGs with pairwise overlaps of 3 (envPC1∩2), 1 (envPC1∩3), 3 (envPC2∩3), and 0 (all three) — the inclusion-exclusion union is exactly 337 − 3 − 1 − 3 + 0 = **330**, confirming the abstract's number is the *union of all three envPC-specific candidate sets*, not a separately computed/stored figure. This union isn't materialized as its own `output/` file, but it's a trivial derived quantity from the three existing `candidateOG_envPC{1,2,3}.txt` files plus the pairwise-overlap counts hardcoded in `src/fig6a_sankey_v2.py` — reproducible, just not pre-computed anywhere.
- **17 high-confidence candidates**: unchanged from the first pass — `src/fig6a_sankey_v2.py` hardcodes envPC1:10 + envPC2:6 + envPC3:1 = 17, matching `output/figure/fig6a_sankey.png` and now also matching the manuscript's Figure 6B table of 17 named OGs (α-expansin/ZmEXPA5, RING-E3/OsDCA1, XTH orthologs, the PGAM ortholog, etc.). Still true that these counts are **hardcoded in the plotting script** rather than computed live at plot time from the upstream tables — so the figure is a manually transcribed snapshot of a specific run, not a live view.
- **Three-layer high-confidence filter**: manuscript text describes the filter as (1) phylogenetic-mixed-model significance (empirical p ≤ 0.001) → (2) RELAX selection-intensity shift on the same transition branches → (3) intersection with a curated list of *a priori* abiotic-stress-responsive genes (validated by expression evidence in ≥3 species-studies). This matches the repo's `08C → 09A/09B → 10B` chain identified in the first pass; Figure 6B's "Differential expression" column is that third (curated-list/expression-evidence) filter, not a separate DEG analysis.
- **Permulation-based empirical p-values**: Methods cites Saputra et al. 2021 "phylogenetic permulations" as the method underlying every reported empirical p-value (genome-wide trait models, per-OG models, and implicitly the RELAX/candidate calls) — consistent with `output/permulation/`, `log/perOGModeling_envPC*_permulation.log`, and `src/12_runPermulation_perOGModel.R`.

### 3.3 — New information from the Data Availability Statement (changes the framing of §4.1)
The manuscript's Data Availability Statement says:
> "Upon the publication of this manuscript, the associated raw sequence data will be available through [SRA] and the assemblies generated will be released through Ag Data Commons. Source scripts, analytical notebooks and derived datasets used will be available at https://github.com/maize-genetics/p_phyloGWAS. Derived occurrence datasets from BIEN and GBIF.org are available at https://zenodo.org/records/14968186."

Two things follow from this that change how §4.1 should be read:
1. **The raw sequence/assembly data embargo is intentional and time-boxed**, not a repo-hygiene failure — SRA/Ag Data Commons releases are explicitly tied to publication. Until then, no fresh clone (by anyone, including the authors on a clean machine) can pull the raw genomes; this is a *publication-timing* blocker, separate from the *code-hygiene* blockers in §4.2–4.4, and should be tracked/communicated separately (e.g. it will resolve itself at publication without any repo changes).
2. **The BIEN/GBIF occurrence data already has a public, citable location** (Zenodo record 14968186) that is **not referenced anywhere in this repo** — not in the README, not in `06_envirotyping/`, not in a DATA.md. This is a quick, concrete fix: add the Zenodo link next to the `06B_spCoordEnvData.sh`/`08_pulling_geo_data.R` step so a reader doesn't have to hunt for it in the manuscript's Data Availability section. Recommended addition to §7.

## 4. Reproducibility blockers, ranked by how much they hurt an outside reader

### 4.1 — Data acquisition step is broken (partly by design, partly not — see §3.3)
The README says: `bash notebook/00_gatherData.sh` to "download all data," pointing to a private Google Sheet on `cbsublfs1`. **`notebook/00_gatherData.sh` does not exist in the repo**, and no `environment.yml`/`requirements.txt`/conda spec exists either, despite the README's `conda create -n <blank>` instructions. A fresh clone has no way to obtain the 707 genome assemblies or any intermediate data.

Per the manuscript's Data Availability Statement (§3.3), part of this is **intentional and temporary**: raw sequence/assembly release to SRA/Ag Data Commons is explicitly tied to publication, so this specific gap should close on its own once the paper is out — it's not something to fix in the repo now. What *isn't* time-boxed, and should be fixed regardless: (a) the dead `00_gatherData.sh` reference and missing environment spec, and (b) the fact that the BIEN/GBIF occurrence data already has a public Zenodo location (record 14968186) that nothing in the repo points to.

### 4.2 — Silent dependency on repos that aren't this repo
Notebooks/scripts across stages 02, 03, 06, 08, 09, and 10 read files from three sibling directories that are **not part of `p_phyloGWAS`**:
- `p_phyloGWAS_archived` (metadata QC, miniprot GFFs, BIEN coordinates, an OG→Zm-gene mapping file used repeatedly downstream)
- `p_evolBNI` (raw WorldClim raster files)
- `p_panAndOGASR` (a GO-term database file used in RELAX/candidate GO enrichment)

All three **do exist on this machine**, right next to `p_phyloGWAS/` (`/local/workdir/sh2246/p_evolBNI`, `/local/workdir/sh2246/p_phyloGWAS_archived`), so the pipeline is reproducible *on this exact filesystem*. But `git clone git@github.com:maize-genetics/p_phyloGWAS` (the README's own instructions) would not bring any of these along, and nothing in the repo documents that they're needed, what they contain, or where to get them. This is the single biggest gap between "runs on my machine" and "reproducible from the public repo."

### 4.3 — Hardcoded absolute paths and machine-specific infrastructure
Virtually every script/notebook hardcodes `/workdir/sh2246/p_phyloGWAS/...`. Stage 07 (ESM2/PlantCAD scoring, premature-stop screening) and stage 09 (HyPhy) additionally hardcode:
- SCINET SLURM account `buckler_lab_panand`, `--partition=gpu-a100`
- Tool paths only present on that cluster: `/programs/hyphy-2.5.49/bin/hyphy`, `/programs/seqkit-0.15.0/`
- User-specific conda activation: `~/miniconda3/bin/activate test`, `.../activate cbsuenv`

None of this is parameterized (no config file, no env vars, no CLI flags for base path).

### 4.4 — Numbering gaps / orphaned scripts
There is no stage "00" (besides the missing `00_gatherData.sh`), "01", or "04" notebook. But `src/01_runGWAS.R`, `src/01_MSA2Dist.R`, `src/01B_extractDistToB73.R`, `src/04_ATGdiagnosis.py`, and `src/04_run_get4d.msa.R` exist with no notebook/shell script in the current (non-archived) tree that calls them. `01_runGWAS.R` in particular looks like it fits an early GWAS-fitting step but its caller is missing or was folded into `XX_archived/` predecessors. Worth clarifying whether these are truly dead code or missing orchestration.

### 4.5 — Proprietary/licensed dependency
The core per-OG statistical model (`src/12_runPermulation_perOGModel.R`, called by `08C_perOGmodel.sh`) requires **ASReml-R**, a commercial/licensed package (VSNi). Anyone without a license cannot run the central mixed-model analysis at all, regardless of data access.

### 4.6 — Undocumented manual curation
`02B_furtherFilter.ipynb` hardcodes a curator's exclusion list (`cohBadSamples <- c(...)`) and a TABASCO completeness cutoff (2796) with no comment explaining the rationale or a sensitivity check. OG filtering thresholds (≥8 taxa, ≤100 copies/taxon) in `03B_OGFilter.ipynb` are similarly bare constants.

### 4.7 — Notebook/version hygiene
`.ipynb_checkpoints/` are committed alongside the real notebooks throughout (adds noise, occasionally drifts from the real file — e.g. `06C_visualizationEnvAdapt.ipynb` currently shows as modified in git while its checkpoint is also modified, so it's mid-edit). Several stages have both an active notebook and 2-4 "archived" predecessors with overlapping scope (e.g. `08_linearModeling/archived/08B_perOGModeling.ipynb` vs. the current `08C_perOGmodel.sh`+`src/12_...R`), and there's no changelog or comment marking which superseded which or why.

## 5. What's actually in good shape

- The bulk genome assemblies (727 `.fa.gz` files in `data/assemblies/`) and major intermediate/result tables are physically present locally, not just referenced — this isn't a repo of dangling paths, it's a repo of *undeclared* paths.
- `log/` contains real run logs (e.g. `perOGModeling_envPC*.log`, up to 388 KB) that corroborate the pipeline was actually executed as described, not just scripted.
- Output file naming is fairly consistent and traceable to a specific figure/table (Fig1-Fig8, TableS5, `candidateOG_*`, `RELAX_resultTable_*`, `ASREML_res_*`), which made it possible to map most figures back to a generating notebook/script.
- The stage-numbering convention (02→10) is followed consistently enough that the overall pipeline order is easy to reconstruct even without a written DAG.

## 6. Bottom line: can this reproduce all the manuscript's results?

**Not from the public repo alone, no — though less of that is "broken code" than the first pass suggested, and more of it is either intentional pre-publication embargo or straightforward filename cleanup.** The code that *implements* essentially every stage does exist in `notebook/`+`src/`, and with the full manuscript now in hand, I could trace **every** headline number I checked (707/569, 687, 22,363/19,613, envPC variance %, 330-gene union, 17 final candidates) back to a specific `output/` file or a hardcoded value in a specific script — nothing is unaccounted for. But actually re-running the pipeline end-to-end — even on a machine with unlimited compute — is blocked by: raw genome/read data that's embargoed until publication (by design, per the Data Availability Statement — not a bug), silent reliance on two or three sibling repos that live outside `p_phyloGWAS` and aren't mentioned anywhere in the README, hardcoded single-user/single-cluster paths, a licensed statistics package (ASReml-R v4.2.0), and some drift between figure numbers in `output/figure/` and the manuscript's current Figure 1–6 layout (§3.1). On *this specific machine*, with the sibling directories present at their current paths, most stages would probably run today if someone patched the hardcoded paths and had ASReml — but that's "reproducible for the original author on this filesystem," not "reproducible by a reader of the paper" once the data embargo lifts.

## 7. Suggested cleanup priorities (not performed — for your review)

1. Either commit a real `00_gatherData.sh` / document the actual current data source, or write a short `DATA.md` stating what's expected in `data/`, where it currently lives, and — critically — **link the Zenodo BIEN/GBIF occurrence record (14968186) and note the SRA/Ag Data Commons publication-embargo status** so a reader isn't left guessing which gaps are permanent and which resolve at publication.
2. Inline or clearly vendor the specific files pulled from `p_phyloGWAS_archived`/`p_evolBNI`/`p_panAndOGASR` — right now the dependency is invisible until a script fails on a missing path.
3. Parameterize the base path (env var or config) instead of hardcoding `/workdir/sh2246/p_phyloGWAS/` in ~every script.
4. Add an `environment.yml`/`renv.lock` seeded from the tool/version table in §2 (flag the ASReml-R v4.2.0 license requirement explicitly since it can't be captured by a lockfile).
5. Register `src/ASTER` as a proper git submodule (it's currently an unregistered nested repo, which is why it shows as dirty in `git status`).
6. Add one line per archived notebook/script noting what superseded it, or move truly dead ones out of the numbered stage folders entirely.
7. Rename or regenerate `output/figure/*.png` so filenames match the manuscript's current Figure 1–6 numbering (§3.1) — at minimum, rename the power-simulation output from `Fig8_revised.png` to `Fig2*.png` and double check `Fig2.png`/`Fig3c`/`Fig3d` from `08A_masterDataTableGeneration.ipynb` against what Figure 3 actually shows now.
8. Reconcile the "34 representative assemblies" comment in `03A_buildHelixerOG.sh` with the Methods-reported 32 (§3.2) — either the script is stale or two assemblies were dropped later; a one-line comment would save a future reader the cross-check.
