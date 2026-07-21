# The Genetic basis of adaptive evolution in Poaceae using Phylogenetic mixed model
This is the repository to store the notebooks and source scripts for cross-species GWAS projects (Life history & env. adaptation).

## Environment setup and data loading
- Fork this repo
```
git clone git@github.com:maize-genetics/p_phyloGWAS
```

- Conda environments (see `envs/`) — split by domain because pinning CLI tools, R, and the
  Python ML stack into one solve is slow/fragile. Activate whichever one a given stage needs
  (see `WORKFLOW.md`'s "Tools used across stages" table):
```
conda env create -f envs/environment-tools.yml   # mafft, orthofinder, hyphy, seqkit, miniprot, gffread, parallel
conda env create -f envs/environment-r.yml       # R 4.2 + envirotyping/modeling/tree packages
conda env create -f envs/environment-py.yml      # Python 3.11 + torch/transformers (ESM2/PlantCAD scoring)
```
  RAxML (8.2.12) and ASTER (astral-pro3) are **not** in any environment file — they're vendored
  as source in `src/standard-RAxML` and `src/ASTER` and need to be compiled once (see each
  directory's own build instructions).

- Licensed software (not distributed via conda, not covered by any `environment*.yml`):
  **ASReml-R**, used in stage 08 (`08C_perOGmodel.sh`, power simulation). The installer is
  vendored at `src/asreml_4.2.0.302_linux-intel64_R4.2.0.tar.gz` (built for R 4.2.0, matching
  `envs/environment-r.yml`'s `r-base` pin) but requires your own VSN International license to
  activate — see https://vsni.co.uk/software/asreml-r for licensing. Install it into the
  `phyloGWAS-r` conda env after creating it:
  ```
  conda activate phyloGWAS-r
  R CMD INSTALL src/asreml_4.2.0.302_linux-intel64_R4.2.0.tar.gz
  R -e 'install.packages("asremlPlus", repos="https://cloud.r-project.org")'
  ```

- Key data can be found on cbsublfs1: https://docs.google.com/spreadsheets/d/1XzyAfph4fuamt4h22uD1OMzkXGQXlFWmwgw7zDncKKY/edit?usp=sharing
  — see `DATA.md` for the full input-data manifest and `WORKFLOW.md` for how each stage consumes it.

## Key steps in these projects include
1. assembly and QC of the genomes from Ag Data Common
2. manual curation of the species identity of the genomes and life history habits
3. collection of species occurrence coordinates and habitat environment characterization (rbien and rgbif -> env GIS data -> quantiles summary -> envPC variables)
4. OG construction and miniProt search (orthoFinder among 32 high-quality representative genomes -> ancestral protein sequence reconstruction -> miniProt query against all assemblies)
5. MSA generation (mafft)
6. Species tree construction -> phyloK estimates (angiosperm 353 loci -> 353 gene trees (pruned to 1 tips per taxa) -> species tree
7. per OG summary stats:
   (1) premature stop codon calling 
   (2) tip-to-outgroup dN/dS calculation
   (3) avg. residue/nucleotide confidence score per sequence (ESM2 & PlantCAD)
   
8. Linear modeling: ASREML & ASREMLPlus
   (1) genome-wide features
   (2) per OG
   (3) power simulation

9. HyPhy - RELAX
10. Candidate gene investigation


# Contact
Sheng-Kai Hsu (sh2246@cornell.edu), Aimee Schulz

