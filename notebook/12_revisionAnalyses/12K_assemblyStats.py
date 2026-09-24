#!/usr/bin/env python3
"""Standard assembly statistics for one gzipped FASTA; prints one TSV row.

Called per assembly by 12K_assemblyStats.sh.

  A. as-deposited sequences: count, total/min/mean/max length, N50/L50,
     N90/L90, GC% over A/C/G/T only (Ns and IUPAC codes excluded from the
     denominator), N count and N%.
  B. contigs: sequences split at runs of >= MIN_GAP Ns (10, the QUAST
     convention), with the gap runs removed. Shorter N runs stay inside the
     contig. For ungapped assemblies B equals A.
  C. count and % of total length in sequences >= 1 kb and >= 10 kb.

Pass 1 (every assembly): `seqkit fx2tab` emits per-sequence length and G+C /
A+C+G+T / N counts (case-insensitive) in one Go pass. Pass 2 (only assemblies
with Ns): `seqkit seq -s -w 0` streams one sequence per line and gaps are
located with numpy. Doing the base counting in Python instead cost ~6x the
runtime (repeated bytes.count scans), so keep it in seqkit.

Note seqkit's own `stats` column N50_num is NOT L50; L50 here is the number
of sequences needed to reach half the total length.

Usage: 12K_assemblyStats.py <assembly.fa.gz>
       12K_assemblyStats.py --header
"""
import re
import subprocess
import sys

import numpy as np
import pandas as pd

SEQKIT = "/home/sh2246/miniconda3/envs/daily/bin/seqkit"
MIN_GAP = 10

COLUMNS = [
    "assemblyID",
    "num_seqs", "total_length", "min_length", "mean_length", "max_length",
    "N50", "L50", "N90", "L90",
    "GC_pct", "N_count", "N_pct",
    "contig_count", "contig_N50", "contig_L50", "contig_N90", "contig_L90",
    "num_seqs_ge1kb", "pct_length_ge1kb", "num_seqs_ge10kb", "pct_length_ge10kb",
]


def nx(lengths, frac):
    """Nx and Lx: the length at which sorted (desc) sequences reach frac of total."""
    s = np.sort(lengths)[::-1]
    idx = int(np.searchsorted(np.cumsum(s), frac * s.sum()))
    return int(s[idx]), idx + 1


def base_counts(path):
    """Pass 1: per-sequence length, G+C, A+C+G+T and N counts via seqkit."""
    proc = subprocess.Popen(
        [SEQKIT, "fx2tab", "-j", "2", "-n", "-i", "-l", "-C", "GC", "-C", "ACGT", "-C", "N", path],
        stdout=subprocess.PIPE)
    tab = pd.read_csv(proc.stdout, sep="\t", header=None, usecols=[1, 2, 3, 4],
                      names=["len", "gc", "acgt", "n"], dtype="int64")
    if proc.wait() != 0:
        raise RuntimeError(f"seqkit fx2tab failed on {path}")
    return tab


def contig_lengths(path):
    """Pass 2: split each sequence at runs of >= MIN_GAP N/n, drop the runs."""
    out = []
    proc = subprocess.Popen([SEQKIT, "seq", "-s", "-w", "0", "-j", "2", path],
                            stdout=subprocess.PIPE, bufsize=1 << 24)
    for line in proc.stdout:
        seq = line.rstrip(b"\n")
        if b"N" not in seq and b"n" not in seq:
            out.append(len(seq))
            continue
        # Work from the N positions only (typically a tiny fraction of the
        # sequence) rather than full-length boolean/diff arrays, which peaked
        # at ~14 GB on a 556 Mb chromosome. (b | 0x20) == 'n' matches N and n.
        pos = np.flatnonzero((np.frombuffer(seq, dtype=np.uint8) | 0x20) == 110)
        brk = np.flatnonzero(np.diff(pos) > 1)
        rs = pos[np.concatenate(([0], brk + 1))]
        re_ = pos[np.concatenate((brk, [pos.size - 1]))] + 1
        keep = (re_ - rs) >= MIN_GAP
        gap_s, gap_e = rs[keep], re_[keep]
        seg_s = np.concatenate(([0], gap_e))
        seg_e = np.concatenate((gap_s, [len(seq)]))
        seg = seg_e - seg_s
        out.extend(seg[seg > 0].tolist())
    proc.stdout.close()
    if proc.wait() != 0:
        raise RuntimeError(f"seqkit seq failed on {path}")
    return np.asarray(out, dtype=np.int64)


def main(path):
    tab = base_counts(path)
    L = tab["len"].to_numpy()
    total, n_count = int(L.sum()), int(tab["n"].sum())
    contigs = contig_lengths(path) if n_count > 0 else L

    n50, l50 = nx(L, 0.5)
    n90, l90 = nx(L, 0.9)
    cn50, cl50 = nx(contigs, 0.5)
    cn90, cl90 = nx(contigs, 0.9)
    ge1, ge10 = L >= 1_000, L >= 10_000

    asm_id = re.sub(r"\.(fa|fasta|fna)\.gz$", "", path.rsplit("/", 1)[-1])
    row = [
        asm_id,
        len(L), total, int(L.min()), round(total / len(L), 1), int(L.max()),
        n50, l50, n90, l90,
        round(100 * tab["gc"].sum() / tab["acgt"].sum(), 3), n_count, round(100 * n_count / total, 4),
        len(contigs), cn50, cl50, cn90, cl90,
        int(ge1.sum()), round(100 * L[ge1].sum() / total, 3),
        int(ge10.sum()), round(100 * L[ge10].sum() / total, 3),
    ]
    print("\t".join(map(str, row)))


if __name__ == "__main__":
    if sys.argv[1] == "--header":
        print("\t".join(COLUMNS))
    else:
        main(sys.argv[1])
