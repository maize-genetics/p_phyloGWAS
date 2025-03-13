# The Genetic basis of adaptive and life history evolution in Poaceae using Phylogenetic mixed model
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
1. assembly and QC of the genomes (Aimee, Charlie)
2. manual curation of the species identity of the genomes and life history habits
3. collection of species occurrence coordinates and habitat environment characterization (Sheng-Kai: rbien and rgbif -> env GIS data -> quantiles summary)
4. OG construction and miniProt search (Sheng-Kai: orthoFinder among 32 high-quality representative genomes -> ancestral protein sequence reconstruction -> Aimee: query against all assemblies)
5. MSA generation (Aimee)
6. Species tree construction -> phyloK estimates (Sheng-Kai: angiosperm 353 loci -> 353 gene trees (pruned to 1 tips per taxa) -> species tree
7. per OG summary stats:
 
 (1) premature stop codon calling (Sheng-Kai) 

 (2) tip-to-outgroup dN/dS calculation (Sheng-Kai) 

 (3) avg. residue/nucleotide confidence score per sequence (ESM2: Sheng-Kai; PlantCaduceus: Aimee)

9. Linear modeling: ASREML & ASREMLPlus (Sheng-Kai, Aimee and Charlie: still working on)

 (1) association (genome-wide features; per OG)

 (2) prediction

10. HyPhy - RELAX (Sheng-Kai: under construction; need SCINET codes)
11. Literature curation for stress-inducible genes (Sheng-Kai & teams)


# Contact
Sheng-Kai Hsu (sh2246@cornell.edu), Aimee Schulz, Charlie Hale

