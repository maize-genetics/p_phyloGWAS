# 04 — MSA Generation (mafft)

Per-orthogroup multiple sequence alignment of the CDS sequences assembled in stage
`03_orthogroup`, ahead of gap-stripping and gene-tree construction in stage
`05_phylotreeConstruction`.

Run once per orthogroup (looped across all OG FASTA files coming out of stage 03):

```
mafft --ep 0 --genafpair --maxiterate 1000 <input> > <output>
```

- `--ep 0` — offset value 0, allows longer gaps in the alignment
- `--genafpair` — E-INS-i strategy, suited for sequences with multiple conserved domains
  separated by long, gappy regions (appropriate for CDS sequences across diverged taxa)
- `--maxiterate 1000` — up to 1000 rounds of iterative refinement

`<input>` is a per-OG CDS FASTA; `<output>` follows the `{OG}_mafft.fa` naming convention
consumed downstream (e.g. `output/OrthofinderMAFFT/{OG}_mafft.fa`, read by
`05A_treeConstruction`).
