# Handoff: where everything lives (as of 2026-10-13)

State of the p_phyloGWAS project and its sibling repos after the first revision round and the
archiving pass. Written so the work can be picked up on another machine without this
machine's session history. Details live in `revision-analyses.md` (analyses), `DATA.md`
(data provenance), `WORKFLOW.md` (pipeline) and, on the file server, `README_ARCHIVE.md`.

## Repositories

| Repo | Where | Notes |
|---|---|---|
| p_phyloGWAS | https://github.com/maize-genetics/p_phyloGWAS, `main` | Revision work merged via PR #2 (`revision-analyses`). Large data/outputs are gitignored and live in the file-server archive below. |
| p_phyloGWASGeneLookup | https://github.com/maize-genetics/p_phyloGWASGeneLookup, `main` | Gene lookup dashboard, live at https://maize-genetics.github.io/p_phyloGWASGeneLookup/. Self-contained (data snapshot committed, built 2026-08-18 with the 27 candidates). `build/buildData.py` rebuilds it from p_phyloGWAS outputs. |
| p_lifeHistoryPoaceae | https://github.com/maize-genetics/p_lifeHistoryPoaceae | Work is on stacked branches `stage0-scaffolding` ... `stage6-association` (tip = all work); `master` holds only the initial extraction commit. History was rewritten on 2026-10-13 (`git filter-repo`) to drop four files over GitHub's 100 MB limit; they stay on disk and are gitignored (see commit `e93bc9d`). All commit hashes changed then. |
| p_evolBNI | not a git repo | Backed up only to the file server (below). |

## File server: cbsublfs1:/data1/users/shengkaihsu/

| Path | Contents | Status |
|---|---|---|
| `p_phyloGWAS/` | Permanent archive at git commit `34f208c`: repo + `.git`, T1/T2 data (`DATA.md` tiers), master tables, HyPhy/RELAX runs, results, OrthoFinder results, figures. 869,116 files, 367.7 GB. `README_ARCHIVE.md`, `MANIFEST.tsv`, `MD5SUMS` at its top. | Uploaded and verified 2026-10-12: `md5sum -c` 0 failures. Mirrors the repo layout, so `data/` and `output/` drop straight into a fresh clone. |
| `GIS_env_data.tar.gz` | GIS rasters (11 GB; GSDE and WorldClim as their original .zip downloads, plus FAO-GAEZ, HYSOGs, PanAnd_src_data). Shared with p_evolBNI. | Checked file by file against the unpacked folder; covers all of it. Unpack recipe in `p_phyloGWAS/README_ARCHIVE.md`. |
| `p_phyloGWAS_old/` | Previous, partial archive (2023 - mid-2026). | Superseded; slated for deletion by the author. Its few unique small files were reviewed before deletion. |
| `p_lifeHistoryPoaceae/` | Working-directory copy (outputs, logs, figures git does not hold; scratch `tmp/` excluded). | Check it is present and current before relying on it. |
| `p_evolBNI/` | Backup of the local p_evolBNI. | Being brought up to date 2026-10-13 (`~/archive_check/evolBNI_sync.sh`): missing/newer files copied, outdated version folders removed on the server. Large regenerable items (RNA-seq SAM/trimmed reads, salmon indexes, gentromes, `promoter_mafft`) deliberately not backed up; raw RNA-seq reads are stored separately. |

Access: SSH with your password works; key-based login and `mount_server` do not for this
server. To avoid retyping the password across many commands, open one shared connection:
`ssh -M -S ~/.ssh/cm_cbsublfs1 -fN cbsublfs1`, use `rsync -e "ssh -S ~/.ssh/cm_cbsublfs1"`,
close with `ssh -S ~/.ssh/cm_cbsublfs1 -O exit cbsublfs1`. Inventories and comparison files
from the archiving pass are in `~/archive_check/` (shared home).

## Not in any archive (local to cbsumm26 only)

- `/workdir/sh2246/p_phyloGWAS_archived/` (1.5 TB): older pipeline generations. Holds a few
  files referenced elsewhere: `output/OGTo{Pv,At}_mapping_v2.txt` (also in the server's old
  archive; not used by the current pipeline) and `output/proteinMSAtoRescue/all_missingSeq.fa`
  (input to the `allMissingSeq_*` ESM scores).
- `p_phyloGWAS/data/assemblies/` (504 GB): to be deposited at agDataCommons (embargoed).
- `p_phyloGWAS/output/distPerMSA/` (350 GB), OrthoFinder `WorkingDirectory/` (22 GB), the
  InterProScan 5.76-107.0 and ASTER installs: excluded on purpose; see `README_ARCHIVE.md`.
- Safety copies from the lifeHistory history rewrite, deletable once satisfied:
  `/local/workdir/sh2246/p_lifeHistoryPoaceae_preFilterRepo_20261013.bundle` and
  `/local/workdir/sh2246/p_lifeHistoryPoaceae_bigfiles_aside/`.

## Things that are easy to get wrong

- **Taxon counts.** `envData_707Poaceae_20250804.txt` has 707 taxa; the per-OG models are fit on
  **686** (21 lack gene data: premature-stop genomes). The ASR and power simulation now use the
  **687-taxon tree** (686 + outgroup ASM1935983v1; `output/powerSimulation/tree.nwk`); 561 species.
- **Transitions on the 687 tree:** 43 / 43 / 54 for envPC1/2/3 (707 tree: 46 / 43 / 57).
  Revision analyses #2, #4, #5b and 12F still compare against the 707 baseline (`revision-analyses.md` §10).
- **High-confidence candidates: 27, not 17** (RELAX-layer bug fixed 2026-07-30; §9).
- **K matrix version.** `output/phyloK_728Poaceae_astral_20250407.txt` was rewritten on 2026-07-29,
  after the final per-OG models (`finalModels_20251002/`) were fit; the 2025-04-07 original was in
  the old server archive. Check they are numerically identical if exact reproduction matters.
- **Metadata with assembly stats:** `data/Poaceae_metadata_highErrorFiltered_2025.10.08_assemblyStats.tsv`
  (687 rows; adds scaffold/contig N50 and gene-bearing contig length). Not tracked in git.

## Restoring on a new machine

1. `git clone` the repo(s) from GitHub.
2. Copy untracked data/outputs from the file server into the clone
   (`rsync -a cbsublfs1:/data1/users/shengkaihsu/p_phyloGWAS/ ./ --exclude .git`).
3. Unpack the GIS tarball if environmental stages are rerun; set `GIS_DATA_ROOT`.
4. Fix machine-specific paths (`PHYLOGWAS_ROOT`; p_lifeHistoryPoaceae's `config.yaml`, checked by
   `Rscript src/S00_checkConfig.R config.yaml`; its masterDataTable symlink points into p_phyloGWAS).

## Open items

- Finish the p_evolBNI sync and confirm its checksum step reported 0 differences.
- Delete `p_phyloGWAS_old/` on the server once satisfied.
- Optionally copy `all_missingSeq.fa` into the archive (`output/proteinMSAtoRescue/`).
- p_lifeHistoryPoaceae: set the GitHub default branch (currently `master`, nearly empty) or merge
  the stages into it.
- At publication: Zenodo DOIs for the GitHub repos; fill the "link TBD" rows in `DATA.md`.
