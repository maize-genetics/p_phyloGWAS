#!/usr/bin/env python3
"""Build the file-server archive manifest for p_phyloGWAS (read-only; copies nothing).

Assigns every file under the repo (and the GIS rasters in the sibling repo p_evolBNI) to one
archive category, or to "excluded" / "unclassified":

  repo          git-tracked files + .git/
  T1-raw        GIS rasters + small public references (DATA.md T1 rows kept locally)
  T2            DATA.md T2 rows (regenerable intermediates)
  intermediate  master data tables + HyPhy/RELAX runs (not in DATA.md, kept by choice)
  results       model inputs/outputs, candidate tables, figures, revision outputs
  excluded      deliberately left out (raw assemblies, distPerMSA, superseded runs, strays)
  other         everything else (older dated runs, exploratory outputs, side threads); kept
                because it is small (~5 GB) - reviewed with the author 2026-10-12

Rules are path-prefix matches, first match wins, in the order of RULES below. Excluded
directories are not walked. Writes, under output/archive/:
  MANIFEST.tsv           relpath, category, size_bytes   (relpath relative to the archive root,
                         i.e. p_phyloGWAS/... or p_evolBNI/...; md5 is added at the checksum step)
  filelist_<cat>.txt     one list per category, paths relative to the archive root
  summary.tsv            files and bytes per category
  other_top.tsv          "other" files rolled up to their top-level entry, with sizes

Usage: python notebook/archive/buildArchiveManifest.py
"""
import os
import subprocess
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
WORK = os.path.dirname(REPO)
EVOL = os.path.join(WORK, "p_evolBNI")
OUT = os.path.join(REPO, "output", "archive")

# (category, prefixes relative to the repo root). A prefix ending in "/" matches a directory.
RULES = [
    ("excluded", [
        "data/assemblies/", "output/distPerMSA/", "output/OrthoFinder/", "output/finalModels/",
        "output/functionPredictionPerOG_20250331.txt", "output/functionPredictionPerOG_20250423.txt",
        "output/functionPredictionPerOG_20250504.txt", "output/distToB73_Charlie",
        "output/jonathanOG/", "output/archive/", "Rplots.pdf",
        # third-party tool installs, re-installable (versions recorded in README_ARCHIVE.md)
        "src/interproscan-5.76-107.0/", "src/ASTER/",
    ]),
    ("T2", [  # before "results" so the logits inside revisionAnalyses/ count as T2
        "output/revisionAnalyses/mutagenesisValidation/At_ESM2_650M_logits.npz",
        "output/revisionAnalyses/mutagenesisValidation/maize/maize_ESM2_650M_logits.npz",
        "data/fullSetOGs_240903.txt", "data/DEG_mappingFiles/",
        "output/orthofinderMiniProt/", "output/orthofinderMiniProt_additionalOGs/",
        "output/seqIDmapping.txt", "output/seqIDmapping_additionalOGs.txt",
        "output/orthofinderProteinMSAs/", "output/orthofinderProteinMSAs_fullset_20250710/",
        "output/CDSMSAPerOG/", "output/CDSMSAPerOG_gs/", "output/aminoAcidPerTaxa/",
        "output/CDSperTaxa/", "output/metadataFormalOut/",
    ]),
    ("T1-raw", [
        "data/Zmays_cds.fa", "data/Zea_mays_v5_mrna.fa", "data/Zm00001eb.1.fulldata_goList.txt",
        "data/maxFPKM_v5.csv", "data/Angiosperms353_orysaSequences.fasta", "data/DEG_CDS_fasta/",
        "data/mutagenesisValidation/B73v5/",
    ]),
    ("intermediate", [
        "output/masterDataTable_20251001.txt", "output/masterDataTable_PAVFill_20251001.txt",
        "output/HyPhyResult/", "output/hyphyRuns_warm_cold/", "output/CDSMSAPerOG_HyPhy_20250203/",
        "output/CDSMSAPerOG_HyPhy_rhizome_20260303/", "output/geneTree_allOGs_20250203/",
        "output/labeledGeneTree_warm/", "output/labeledGeneTree_cold/", "output/targetHyPhy_",
        "output/labeledGeneTree_drought/", "output/labeledGeneTree_wet/",
        "output/labeledGeneTree_sand/", "output/labeledGeneTree_clay/",
    ]),
    ("results", [
        "output/PoaceaeTree_", "output/phyloK_728Poaceae_astral_20250407.txt", "output/envData",
        "output/ePC_20250804.rds", "output/KG3_perSpecies_20250804.txt",
        "output/finalModels_20251002/", "output/permulation/", "output/candidateOG_envPC",
        "output/candidateGenes/", "output/finalCandidateOGs_summary.txt", "output/RELAX_resultTable_",
        "output/combined_ESM_zeroShotScores_20250414.txt",
        "output/combined_plantCAD_zeroShotScores_20250421.txt", "output/combined_PMS_20250421.txt",
        "output/functionPredictionPerOG_20251023.txt", "output/TableS", "output/figure/",
        "output/powerSimulation/", "output/assemblyStats/", "output/hitContigLength/",
        "output/revisionAnalyses/",
    ]),
]
EXCLUDED_DIRS = tuple(p for p in RULES[0][1] if p.endswith("/"))


def classify(rel, tracked):
    if "/.ipynb_checkpoints/" in "/" + rel:
        return "excluded"
    if rel in tracked or rel.startswith(".git/"):
        return "repo"
    for cat, prefixes in RULES:
        if any(rel.startswith(p) for p in prefixes):
            return cat
    return "other"


def walk(root, skip_dirs=()):
    """Yield (relpath, size) for every regular file under root, not entering skip_dirs."""
    stack = [""]
    while stack:
        d = stack.pop()
        with os.scandir(os.path.join(root, d)) as it:
            for e in it:
                rel = os.path.join(d, e.name) if d else e.name
                if e.is_dir(follow_symlinks=False):
                    if (rel + "/") not in skip_dirs:
                        stack.append(rel)
                elif e.is_file(follow_symlinks=False) or e.is_symlink():
                    yield rel, e.stat(follow_symlinks=False).st_size


def main():
    os.makedirs(OUT, exist_ok=True)
    tracked = set(subprocess.run(["git", "-C", REPO, "ls-files"], capture_output=True, text=True,
                                 check=True).stdout.split("\n"))
    rows = []
    for rel, size in walk(REPO, EXCLUDED_DIRS):
        rows.append(("p_phyloGWAS/" + rel, classify(rel, tracked), size))
    for d in EXCLUDED_DIRS:  # record excluded directories as single entries
        p = os.path.join(REPO, d)
        if os.path.isdir(p):
            rows.append(("p_phyloGWAS/" + d, "excluded", -1))
    gis = os.path.join(EVOL, "data", "GIS_env_data")
    for rel, size in walk(gis):
        rows.append(("p_evolBNI/data/GIS_env_data/" + rel, "T1-raw", size))

    cats = defaultdict(lambda: [0, 0])
    unc = defaultdict(lambda: [0, 0])
    with open(os.path.join(OUT, "MANIFEST.tsv"), "w") as fh:
        fh.write("relpath\tcategory\tsize_bytes\n")
        for rel, cat, size in sorted(rows):
            fh.write(f"{rel}\t{cat}\t{size}\n")
            cats[cat][0] += 1
            cats[cat][1] += max(size, 0)
            if cat == "other":
                parts = rel.split("/")
                top = "/".join(parts[:3]) if parts[1] in ("data", "output") else "/".join(parts[:2])
                unc[top][0] += 1
                unc[top][1] += size
    for cat in cats:
        if cat == "excluded":
            continue
        with open(os.path.join(OUT, f"filelist_{cat}.txt"), "w") as fh:
            fh.writelines(r[0] + "\n" for r in sorted(rows) if r[1] == cat)
    with open(os.path.join(OUT, "summary.tsv"), "w") as fh:
        fh.write("category\tfiles\tGB\n")
        for cat, (n, b) in sorted(cats.items()):
            fh.write(f"{cat}\t{n}\t{b / 1e9:.2f}\n")
    with open(os.path.join(OUT, "other_top.tsv"), "w") as fh:
        fh.write("entry\tfiles\tGB\n")
        for top, (n, b) in sorted(unc.items(), key=lambda kv: -kv[1][1]):
            fh.write(f"{top}\t{n}\t{b / 1e9:.3f}\n")
    print(open(os.path.join(OUT, "summary.tsv")).read())


if __name__ == "__main__":
    main()
