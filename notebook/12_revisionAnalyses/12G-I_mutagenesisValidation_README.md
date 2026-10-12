> **Archived README of the standalone `protein_LLMscoring_pkg`** (2026-09-17), kept for its run
> commands and environment notes. The package has been deleted; every file it contained is in this
> repo under different paths:
>
> | package path | repo path |
> |---|---|
> | `src/4_ESM_logits.py` | `src/4_ESM_logits.py` |
> | `src/10_logit2zeroShot.R` | `src/10_logit2zeroShot_MutagenesisValidation.R` |
> | `src/visualization.R` | `src/mutagenesisValidation_visualization.R` |
> | `src/compare_datasets.R` | `src/mutagenesisValidation_compareDatasets.R` |
> | `src/maize_01_build_deleterious.py` | `notebook/12_revisionAnalyses/12G_buildDeleteriousSet.py` |
> | `src/maize_02a_fetch_hapmap_regions.sh` | `notebook/12_revisionAnalyses/12H_fetchHapmapRegions.sh` |
> | `src/maize_02b_hapmap_neutral.py` | `notebook/12_revisionAnalyses/12I_buildNeutralSet.py` |
> | `data/` | `data/mutagenesisValidation/` |
> | `output/` | `output/revisionAnalyses/mutagenesisValidation/` |
> | `environment/` | `envs/*mutagenesisValidation*` |
> | `MD5SUMS` | `data/mutagenesisValidation/MD5SUMS_sourcePackage` |
>
> The commands below use the package paths; substitute the repo paths above. See also `DATA.md`
> and `revision-analyses.md` §6.

---

# protein_LLMscoring

Zero-shot scoring of missense mutations with a protein language model
(ESM2, `facebook/esm2_t33_650M_UR50D`), benchmarked on whether the score
separates deleterious from neutral missense mutants relative to their wild-type
protein. Two datasets:

| dataset | genes | deleterious | neutral | PR AUC | ROC AUC |
| --- | --- | --- | --- | --- | --- |
| arabidopsis (TAIR10 subset) | 996 | 2,906 | 1,574 | 0.923 | 0.862 |
| maize (B73v5, experimentally validated) | 35 | 54 | 53 | 0.956 | 0.965 |

Overlay: `output/figure/At_vs_maize_PR_ROC.png`.

## Layout

```
src/
  4_ESM_logits.py                  step 1: per-residue ESM2 probabilities -> .npz
  10_logit2zeroShot.R              step 2: .npz -> per-sequence zero-shot score (TSV)
  visualization.R                  step 3: gene-level contrast, PR/ROC, figures, tables
  compare_datasets.R               overlay PR/ROC of several datasets on one figure
  maize_01_build_deleterious.py    maize: WT + deleterious FASTA from the mutation table
  maize_02a_fetch_hapmap_regions.sh maize: download HapMap3.2.1 282-panel VCFs, extract gene regions
  maize_02b_hapmap_neutral.py      maize: pick common missense SNPs as the neutral set, assemble FASTA
data/
  At.TAIR10.subset.trainingData.fa            arabidopsis WT + _neutralMut_N + _mut_N records
  At_gene2mutate.txt                          gene list from the original analysis (not read by the pipeline)
  deleteriousMissenseMutationsMaize.txt       maize input: protein isoform + missense (e.g. L129F)
  B73v5/                                      peptide FASTA, GFF3, canonical-transcript list, genome .fai
                                              (genome FASTA itself is NOT shipped - see "Maize" below)
  hapmap282/hapmap282_targetGenes.vcf.gz      35-gene extract of the 282-panel VCF (+ .tbi)
output/
  At_ESM2_650M_logits.npz                     step-1 output (270 MB; shipped so no GPU is needed)
  At_ESM2_650M_zeroShotScores.txt             step-2 output
  At_ESM2_650M_summary.txt / _perGeneMeanByLabel.txt / figure/   step-3 output
  maize/                                      same set for maize, plus the prep tables (below)
  logs/
environment/                                  python + R pins, full conda export
MD5SUMS
```

## Environment

- Python: conda env recipe in `environment/pytorch-2-1.conda-export.yml`; minimal
  pins in `environment/requirements-python.txt` (torch 2.1.1, transformers 4.40.0,
  biopython 1.83, numpy 1.26.4). The ESM2 weights are pulled from the Hugging Face
  hub on first use (~2.6 GB, cached in `~/.cache/huggingface`).
- R 4.0.5 with reticulate, Biostrings, limma, PRROC (`environment/requirements-R.txt`).
- `10_logit2zeroShot.R` reads the `.npz` through reticulate. It defaults to
  `/usr/bin/python`; if that python has no numpy, run
  `export RETICULATE_PYTHON=/path/to/python-with-numpy` first.
- samtools / bcftools / tabix / bgzip / curl for the maize prep scripts only.
- GPU: only step 1 benefits. Both `.npz` files are shipped, so steps 2-3 and all
  figures reproduce on CPU. The 142-sequence maize set also runs on CPU with
  `-device cpu` in reasonable time; the 5,476-sequence arabidopsis set takes
  ~7 min on an H100.

## Reproduce

All commands from the repo root.

### Arabidopsis

```bash
# step 1 (only if regenerating the shipped npz)
python src/4_ESM_logits.py -input data/At.TAIR10.subset.trainingData.fa \
    -output output/At_ESM2_650M_logits -model facebook/esm2_t33_650M_UR50D -device cuda:0

# step 2
Rscript src/10_logit2zeroShot.R output/At_ESM2_650M_logits.npz \
    data/At.TAIR10.subset.trainingData.fa output/At_ESM2_650M_zeroShotScores.txt

# step 3  (args: fasta scores outDir [outPrefix] [groupRegex])
Rscript src/visualization.R data/At.TAIR10.subset.trainingData.fa \
    output/At_ESM2_650M_zeroShotScores.txt output At_ESM2_650M
```

### Maize

The genome FASTA (`Zm-B73-REFERENCE-NAM-5.0.fa`, 2.2 GB, MaizeGDB) is needed
only by `maize_02b` and is not in the package. Link your copy in first; the
shipped `.fai` matches it:

```bash
ln -s /path/to/Zm-B73-REFERENCE-NAM-5.0.fa data/B73v5/Zm-B73-REFERENCE-NAM-5.0.fa
```

```bash
# WT + deleterious FASTA, with per-mutation validation table
python src/maize_01_build_deleterious.py -table data/deleteriousMissenseMutationsMaize.txt \
    -peptide data/B73v5/Zm-B73-REFERENCE-NAM-5.0_peptide.fa \
    -out output/maize/maize_WT_deleterious.fa -report output/maize/maize_deleterious_validation.txt

# HapMap region extract (already shipped; rerun only to regenerate). Downloads the
# 10 per-chromosome VCFs (5.75 GB) into data/hapmap282/full/ unless present.
bash src/maize_02a_fetch_hapmap_regions.sh output/maize/target_genes.bed data/hapmap282/hapmap282_targetGenes.vcf.gz

# neutral set + final FASTA (WT, neutralMut_N, mut_N per gene)
python src/maize_02b_hapmap_neutral.py -vcf data/hapmap282/hapmap282_targetGenes.vcf.gz \
    -gff data/B73v5/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz \
    -genome data/B73v5/Zm-B73-REFERENCE-NAM-5.0.fa \
    -peptide data/B73v5/Zm-B73-REFERENCE-NAM-5.0_peptide.fa \
    -deleterious output/maize/maize_deleterious_validation.txt -outdir output/maize

# steps 1-3
python src/4_ESM_logits.py -input output/maize/maize_trainingData.fa \
    -output output/maize/maize_ESM2_650M_logits -model facebook/esm2_t33_650M_UR50D -device cuda:0
Rscript src/10_logit2zeroShot.R output/maize/maize_ESM2_650M_logits.npz \
    output/maize/maize_trainingData.fa output/maize/maize_ESM2_650M_zeroShotScores.txt
Rscript src/visualization.R output/maize/maize_trainingData.fa \
    output/maize/maize_ESM2_650M_zeroShotScores.txt output/maize maize_ESM2_650M

# overlay
Rscript src/compare_datasets.R output/figure/At_vs_maize_PR_ROC.png \
    arabidopsis:data/At.TAIR10.subset.trainingData.fa:output/At_ESM2_650M_zeroShotScores.txt \
    maize:output/maize/maize_trainingData.fa:output/maize/maize_ESM2_650M_zeroShotScores.txt
```

`target_genes.bed` was derived from the GFF3 gene records of the 35 genes in
the mutation table (`chrom start end strand geneID`).

## Data provenance

- `At.TAIR10.subset.trainingData.fa`: from the 2024-06 hackathon analysis. WT
  headers are Ensembl-style; mutants are `<transcript>_neutralMut_N` / `<transcript>_mut_N`.
- `deleteriousMissenseMutationsMaize.txt`: experimentally validated deleterious
  missense mutations, B73v5 protein IDs with isoform (`_P001` etc.). 53/54 match
  the reference protein; `Zm00001eb252380_P001 A139T` is applied as **A140T**
  (position 139 is G, 140 is A - treated as an off-by-one in the source) and is
  flagged `OK_CORRECTED_FROM_A139T` in `output/maize/maize_deleterious_validation.txt`.
- B73v5 references (MaizeGDB, `Zm-B73-REFERENCE-NAM-5.0`, annotation `Zm00001eb.1`):
  GFF3 from https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz;
  peptide FASTA and genome from the same release. Protein `_P###` corresponds to
  transcript `_T###`.
- HapMap3.2.1, 282 association panel (Bukowski et al. 2018, GigaScience), lifted
  to B73v5 and SnpEff-annotated by MaizeGDB (Panzea_HapMap3.2.1_282, 277 samples):
  https://iastate.app.box.com/s/9edkm5f51q88ltzrb0hxaz4p36fdsr7q/folder/364949363629
  (`panzea_hapmap_282_chr{1..10}.vcf.gz` + `.tbi`). Box file ids are recorded in
  `maize_02a_fetch_hapmap_regions.sh`.

## How the score and the benchmark are defined

- **Step 1** runs one unmasked forward pass per sequence (sliding 1022-residue
  windows for longer proteins) and stores, per residue, the softmax over the 20
  canonical amino acids only, in the order `LAGVSERTIDPKQNFYMHWC`.
- **Step 2** score per sequence = mean over standard-AA sites of
  `log( P(observed) / max_{a != observed} P(a) )`. `scaledZeroShotScore` is the
  z-score within gene (gene = first token of the seqID split on `_` or space;
  optional 4th argument to change the regex).
- **Step 3** z-scales `scaledZeroShotScore` within gene again, subtracts the
  gene's WT value, negates so that "more damaging" is larger, and computes PR/ROC
  with deleterious mutants as positives and neutral mutants as negatives, over all
  genes. There is no train/test split: the score is unsupervised, nothing is fit.
- **Maize neutral set** (`maize_02b`): for every biallelic SNP inside the CDS of
  the isoform named in the mutation table, the amino-acid change is re-derived
  from GFF3 CDS + genome on that isoform (the VCF's SnpEff `SUB` field refers to
  whichever isoform SnpEff chose; all 35 CDS translations were checked to be
  identical to the official peptides). Kept if missense, ALT-allele frequency
  >= 0.1 and call rate >= 0.5 computed from the 277 genotypes, and not at a
  position already used by a deleterious mutation. Two per gene are drawn at
  random (`seed 123`); genes with fewer eligible sites keep what exists. All
  candidates with their frequencies, reasons and the selected flag are in
  `output/maize/maize_hapmap_missense_candidates.txt`; per-gene counts in
  `maize_neutral_summary.txt`.

Known caveats for maize: 6 of 35 genes have no eligible neutral site
(`Zm00001eb153340, 190120, 194210, 215710, 234730, 373430`), so their deleterious
mutants enter the PR/ROC without a within-gene neutral comparator; and
`Zm00001eb336930` drew two SNPs in the same codon (N77S, N77D). With 107 mutants
the maize AUCs have wide confidence intervals.

## Integrity

`md5sum -c MD5SUMS` from the repo root checks every shipped file.
