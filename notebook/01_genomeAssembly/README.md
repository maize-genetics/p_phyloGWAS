# 01 — Genome Assembly (megahit)

Short-read genome assembly of the newly sequenced accessions, producing the raw
`*.fa.gz` assemblies later QC'd and curated in stage `02_metadataCuration`.

This step is **not run from within `p_phyloGWAS`** — the assembly pipeline (read QC,
megahit assembly, post-assembly filtering) is more involved than a single command and is
maintained in a separate repository:

**https://bitbucket.org/bucklerlab/p_reelgene/src/master/short_read_assembly/**

Refer to that repo for the actual megahit invocation and surrounding pipeline steps. The
output of that pipeline — one assembly FASTA per accession — is what lands in this repo's
`data/assemblies/` (see `DATA.md`) and feeds stage `02`.
