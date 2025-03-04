# Introduction
This is the repository to store the scripts for cross-species GWAS projects (Life history & env. adaptation).

## Environment configuration and dataloading
- Environment loading
```
# load conda env. from yaml
conda create -n 
# Fork this repo
git clone git@github.com:maize-genetics/dna_language_model.git
```

- Key data can be found on cbsublfs1: https://docs.google.com/spreadsheets/d/1XzyAfph4fuamt4h22uD1OMzkXGQXlFWmwgw7zDncKKY/edit?usp=sharing
```
#download all data
bash 00_gatherData.sh
```

## Key steps in these projects include
1. assembly and QC of the genomes (Aimee, Charlie)
2. manual curation of the species identity of the genomes and life history habits
3. collection of species occurrence coordinates and habitat environment characterization (Sheng-Kai: rbien and rgbif -> env GIS data -> quantiles summary)
4. OG construction and miniProt search (Sheng-Kai: orthoFinder among 32 high-quality representative genomes -> ancestral protein sequence reconstruction -> query against all assemblies)
5. MSA generation (Aimee)
6. Species tree construction -> phyloK estimates (Sheng-Kai: angiosperm 353 loci -> 353 gene trees (pruned to 1 tips per taxa) -> species tree
7. per OG summary stats: 
(1) premature stop codon calling (Sheng-Kai) 
(2) tip-to-outgroup dN/dS calculation (Sheng-Kai) 
(3) avg. residue/nucleotide confidence score per sequence (ESM2: Sheng-Kai; PlantCaduceus: Aimee)
8. Linear modeling: ASREML & ASREMLPlus (Sheng-Kai, Aimee and Charlie)
9. HyPhy - RELAX (under construction)
10. Literature curation for stress-inducible genes (Sheng-Kai & teams)
11. 
