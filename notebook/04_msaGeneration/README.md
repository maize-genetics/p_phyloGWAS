# 04 — MSA Generation (mafft)

Per-orthogroup multiple sequence alignment of the CDS/protein sequences extracted from
stage `03_orthogroup`'s OG set, ahead of gap-stripping and gene-tree construction in stage
`05_phylotreeConstruction`.

## Step 0: extract orthologous sequences per assembly

Before alignment, the orthologous CDS/protein sequence for each OG needs to be pulled out
of every assembly, using the miniprot GFF cross-mapping produced in stage `03_orthogroup`
(`03C_ancestralSeqReconstruction.sh`'s miniprot ID-matching block). This is done with the
vendored `patch-scaffolds` tool (`src/patch-scaffolds`):

```
src/patch-scaffolds/bin/patch-scaffolds build-msa \
  --input-fasta-dir <assemblies dir> \
  --input-gff-dir <miniprot GFF dir> \
  --output-dir <output dir> \
  --delimiter ":" \
  [--single-fasta] [--protein-mode] \
  --mapping-file-name <MPID-to-OGID mapping output>
```

- `--input-fasta-dir` / `--input-gff-dir` / `--output-dir` — self-explanatory; GFFs are the
  per-assembly miniprot results from stage 03
- `--delimiter ":"` — delimiter used in the output sequence IDs
- `--single-fasta` — combine everything into one FASTA instead of one file per OG
- `--protein-mode` — output protein sequences rather than CDS
- `--mapping-file-name` — writes an MPID→OGID mapping table alongside the extracted sequences

`--single-fasta` and `--protein-mode` are independent flags, run in all 4 combinations for
different downstream needs:

| `--protein-mode` | `--single-fasta` | Output | Feeds |
|---|---|---|---|
| off (CDS) | off | per-OG CDS FASTAs | MAFFT alignment (below) |
| off (CDS) | on | one combined nucleotide FASTA | PlantCAD zero-shot scoring (stage 07) |
| on | off | per-OG protein FASTAs | MAFFT alignment (below) |
| on | on | one combined protein FASTA | ESM2 zero-shot scoring (stage 07) |

**Requires Java 17+** — the vendored jar is compiled for class file version 61 (Java 17); this
machine's default `java` resolves to OpenJDK 13 and fails with `UnsupportedClassVersionError`.
Point `JAVA_HOME`/`PATH` at a Java 17+ install before running (e.g. `/programs/jdk-17.0.10` on
this machine).

## Step 1: per-OG alignment

Run once per orthogroup (looped across all OG FASTA files coming out of stage 03):

```
mafft --ep 0 --genafpair --maxiterate 1000 <input> > <output>
```

- `--ep 0` — offset value 0, allows longer gaps in the alignment
- `--genafpair` — E-INS-i strategy, suited for sequences with multiple conserved domains
  separated by long, gappy regions (appropriate for CDS sequences across diverged taxa)
- `--maxiterate 1000` — up to 1000 rounds of iterative refinement

`<input>` is a per-OG FASTA from step 0 (CDS or protein, per the table above); `<output>`
follows the `{OG}_mafft.fa` naming convention consumed downstream (e.g.
`output/OrthofinderMAFFT/{OG}_mafft.fa`, read by `05A_treeConstruction` for the CDS case).
