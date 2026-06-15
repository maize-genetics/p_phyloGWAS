# The Genetic basis of adaptive evolution in Poaceae using Phylogenetic mixed model
This is the repository to store the notebooks and source scripts for cross-species GWAS projects (Life history & env. adaptation).

## Environment setup and data loading
- Environment loading
```
# load conda env. from yaml
conda create -n 
# Fork this repo
git clone git@github.com:maize-genetics/p_phyloGWAS
```

- Key data can be found on cbsublfs1: https://docs.google.com/spreadsheets/d/1XzyAfph4fuamt4h22uD1OMzkXGQXlFWmwgw7zDncKKY/edit?usp=sharing
```
#download all data
bash notebook/00_gatherData.sh
```

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

