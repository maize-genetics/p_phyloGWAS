#!/usr/bin/env python3
"""Host-contig length of every miniprot hit in one assembly; writes a per-hit TSV.gz.

Called per assembly by 12L_hitContigLength.sh.

For each mRNA record in the assembly's miniprot GFFs (the main OG run plus the
additional-OG run), find the contig that carries it. Contigs are sequences split
at runs of >= MIN_GAP Ns, the same rule as 12K, so for scaffolded assemblies a
hit's host is the gap-free stretch around it rather than the whole scaffold. A
hit that crosses a gap is flagged spans_gap; its contig_len then covers the
first-to-last contig it touches (the scaffold-joined stretch).

Only hit-bearing sequences are read at base level (pass 2), and only if they
contain Ns; sequence lengths for everything else come from seqkit fx2tab.

Output columns:
  assemblyID OG rank seqid hit_start hit_end hit_span q_len q_cov identity
  seq_len contig_len dist_left dist_right spans_gap

hit_span is the genomic span (introns included), q_cov is (target end - start
+ 1) / query length, dist_left/right are bases between the hit and the
contig ends.

Usage: 12L_hitContigLength.py <assembly.fa.gz> <queryLengths.tsv> <outdir> <gff> [<gff> ...]
"""
import gzip
import os
import re
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

SEQKIT = "/home/sh2246/miniconda3/envs/daily/bin/seqkit"
MIN_GAP = 10
ATTR = re.compile(r"Rank=(\d+);Identity=([\d.]+).*Target=(\S+) (\d+) (\d+)")


def read_hits(gffs):
    rows = []
    for gff in gffs:
        with open(gff) as fh:
            for line in fh:
                f = line.split("\t", 8)
                if len(f) < 9 or f[2] != "mRNA":
                    continue
                m = ATTR.search(f[8])
                rows.append((m.group(3), int(m.group(1)), f[0], int(f[3]), int(f[4]),
                             float(m.group(2)), int(m.group(5)) - int(m.group(4)) + 1))
    return pd.DataFrame(rows, columns=["OG", "rank", "seqid", "hit_start", "hit_end",
                                       "identity", "q_aln"])


def seq_table(path):
    """Pass 1: length and N count of every sequence."""
    proc = subprocess.Popen([SEQKIT, "fx2tab", "-j", "2", "-n", "-i", "-l", "-C", "N", path],
                            stdout=subprocess.PIPE)
    tab = pd.read_csv(proc.stdout, sep="\t", header=None, usecols=[0, 1, 2],
                      names=["seqid", "seq_len", "n"], dtype={"seqid": str})
    if proc.wait() != 0:
        raise RuntimeError(f"seqkit fx2tab failed on {path}")
    return tab


def gaps(path, ids):
    """Pass 2: {seqid: (gap_starts, gap_ends)} (0-based, half-open) for the given ids."""
    out = {}
    with tempfile.NamedTemporaryFile("w", suffix=".ids") as tf:
        tf.write("\n".join(ids) + "\n")
        tf.flush()
        grep = subprocess.Popen([SEQKIT, "grep", "-j", "2", "-f", tf.name, path],
                                stdout=subprocess.PIPE)
        tab = subprocess.Popen([SEQKIT, "fx2tab", "-i"], stdin=grep.stdout,
                               stdout=subprocess.PIPE, bufsize=1 << 24)
        grep.stdout.close()
        for line in tab.stdout:
            sid, seq = line.rstrip(b"\n").split(b"\t")[:2]
            # Same N-position-only approach as 12K (see contig_lengths there).
            pos = np.flatnonzero((np.frombuffer(seq, dtype=np.uint8) | 0x20) == 110)
            brk = np.flatnonzero(np.diff(pos) > 1)
            rs = pos[np.concatenate(([0], brk + 1))]
            re_ = pos[np.concatenate((brk, [pos.size - 1]))] + 1
            keep = (re_ - rs) >= MIN_GAP
            out[sid.decode()] = (rs[keep], re_[keep])
        if tab.wait() != 0 or grep.wait() != 0:
            raise RuntimeError(f"seqkit grep/fx2tab failed on {path}")
    return out


def place(hits, seq_len, gap):
    """Contig coordinates for hits on one sequence."""
    s0 = hits["hit_start"].to_numpy() - 1          # 0-based, half-open
    e0 = hits["hit_end"].to_numpy()
    if gap is None or gap[0].size == 0:
        z = np.zeros(len(hits), dtype=bool)
        return np.full(len(hits), seq_len), s0, seq_len - e0, z
    gs, ge = gap
    seg_s = np.concatenate(([0], ge))
    seg_e = np.concatenate((gs, [seq_len]))
    i = np.searchsorted(seg_s, s0, side="right") - 1   # segment holding the start
    j = np.searchsorted(seg_s, e0 - 1, side="right") - 1  # segment holding the end
    # A start inside a gap belongs to the next contig; an end inside a gap
    # already resolves to the previous one via j.
    i = np.where(s0 >= seg_e[i], i + 1, i)
    cs, ce = seg_s[i], seg_e[j]
    return ce - cs, np.maximum(s0 - cs, 0), np.maximum(ce - e0, 0), j > i


def main(fa, qlen_path, outdir, gffs):
    asm = re.sub(r"\.(fa|fasta|fna)\.gz$", "", os.path.basename(fa))
    hits = read_hits(gffs)
    qlen = pd.read_csv(qlen_path, sep="\t", header=None, names=["OG", "q_len"])
    hits = hits.merge(qlen, on="OG", how="left")
    hits["q_cov"] = (hits["q_aln"] / hits["q_len"]).round(4)

    seqs = seq_table(fa)
    seqs = seqs[seqs["seqid"].isin(set(hits["seqid"]))].set_index("seqid")
    gapped = seqs.index[seqs["n"] >= MIN_GAP].tolist()
    gap = gaps(fa, gapped) if gapped else {}

    parts = []
    for sid, h in hits.groupby("seqid", sort=False):
        L = int(seqs.at[sid, "seq_len"])
        cl, dl, dr, sg = place(h, L, gap.get(sid))
        parts.append(h.assign(seq_len=L, contig_len=cl, dist_left=dl, dist_right=dr,
                              spans_gap=sg.astype(int)))
    out = pd.concat(parts)
    out.insert(0, "assemblyID", asm)
    out["hit_span"] = out["hit_end"] - out["hit_start"] + 1
    cols = ["assemblyID", "OG", "rank", "seqid", "hit_start", "hit_end", "hit_span", "q_len",
            "q_cov", "identity", "seq_len", "contig_len", "dist_left", "dist_right", "spans_gap"]
    with gzip.open(os.path.join(outdir, f"{asm}.tsv.gz"), "wt") as fh:
        out[cols].to_csv(fh, sep="\t", index=False)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
