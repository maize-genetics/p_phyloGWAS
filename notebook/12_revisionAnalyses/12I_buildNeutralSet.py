"""Pick common missense variants from the HapMap3.2.1 282 panel as the neutral
set for the maize genes, and assemble the final training FASTA (WT +
neutralMut + mut per gene, mirroring the At layout).

Usage:
  python src/maize_02b_hapmap_neutral.py \
      -vcf data/hapmap282/hapmap282_targetGenes.vcf.gz \
      -gff data/B73v5/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz \
      -genome data/B73v5/Zm-B73-REFERENCE-NAM-5.0.fa \
      -peptide data/B73v5/Zm-B73-REFERENCE-NAM-5.0_peptide.fa \
      -deleterious output/maize/maize_deleterious_validation.txt \
      -outdir output/maize \
      [-minAltFreq 0.1] [-minCallRate 0.5] [-nPerGene 2] [-seed 123]

The protein change of every SNP is re-derived on the *target isoform* (the one
named in the deleterious table) from GFF3 CDS + genome, rather than trusting the
VCF's SnpEff SUB field, which is relative to whichever isoform SnpEff picked.
ALT-allele frequency and call rate are recomputed from the 277 genotypes.
"""
import argparse, gzip, os, random, subprocess
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-vcf", required=True)
    p.add_argument("-gff", required=True)
    p.add_argument("-genome", required=True)
    p.add_argument("-peptide", required=True)
    p.add_argument("-deleterious", required=True, help="validation table from maize_01")
    p.add_argument("-outdir", required=True)
    p.add_argument("-minAltFreq", type=float, default=0.1)
    p.add_argument("-minCallRate", type=float, default=0.5)
    p.add_argument("-nPerGene", type=int, default=2)
    p.add_argument("-seed", type=int, default=123)
    return p.parse_args()

def faidx(genome, chrom, start, end):
    out = subprocess.run(["samtools", "faidx", genome, f"{chrom}:{start}-{end}"],
                         capture_output=True, text=True, check=True).stdout
    return "".join(out.split("\n")[1:]).upper()

def load_cds(gff, protein_ids):
    """protein_id -> dict(chrom, strand, exons=[(start,end,phase)] in genomic order)"""
    cds = defaultdict(lambda: {"exons": []})
    with gzip.open(gff, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "CDS": continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            pid = attrs.get("protein_id")
            if pid not in protein_ids: continue
            cds[pid]["chrom"] = f[0]; cds[pid]["strand"] = f[6]
            cds[pid]["exons"].append((int(f[3]), int(f[4]), int(f[7])))
    for pid in cds:
        cds[pid]["exons"].sort()
    return cds

def build_cds_seq(genome, rec):
    """Return spliced CDS in transcript orientation plus a map: genomic pos -> CDS index (0-based)."""
    exons = rec["exons"]
    seq_parts, pos_map = [], {}
    idx = 0
    order = exons if rec["strand"] == "+" else exons[::-1]
    for (s, e, _) in order:
        chunk = faidx(genome, rec["chrom"], s, e)
        if rec["strand"] == "-":
            chunk = str(Seq(chunk).reverse_complement())
            positions = range(e, s - 1, -1)
        else:
            positions = range(s, e + 1)
        for gpos in positions:
            pos_map[gpos] = idx; idx += 1
        seq_parts.append(chunk)
    return "".join(seq_parts), pos_map

def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    random.seed(a.seed)

    # ---- deleterious table -> target isoforms and occupied protein positions ----
    targets = {}        # protein_id -> set(protein positions already used by deleterious muts)
    with open(a.deleterious) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f[6].startswith("OK"): continue
            targets.setdefault(f[0], set()).add(int(f[3][1:-1]))
    protein_ids = set(targets)
    gene_of = {pid: pid.split("_")[0] for pid in protein_ids}

    pep = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(a.peptide, "fasta")}

    # ---- CDS models + consistency check against peptide FASTA ----
    cds = load_cds(a.gff, protein_ids)
    models = {}
    for pid in sorted(protein_ids):
        seq, pos_map = build_cds_seq(a.genome, cds[pid])
        prot = str(Seq(seq).translate()).rstrip("*")
        ref = pep[pid.replace("_P", "_T")]
        assert prot == ref, f"CDS translation != peptide FASTA for {pid}"
        models[pid] = {"cds": seq, "pos_map": pos_map, "chrom": cds[pid]["chrom"], "prot": prot}
    print(f"CDS models built and verified against peptide FASTA for {len(models)} isoforms")

    # ---- walk the VCF ----
    q = subprocess.run(["bcftools", "query", "-f",
                        "%CHROM\t%POS\t%REF\t%ALT\t%INFO/TYPE\t%INFO/GENEMODEL\t%INFO/SUB\t%INFO/MAF[\t%GT]\n", a.vcf],
                       capture_output=True, text=True, check=True).stdout
    cand_path = os.path.join(a.outdir, "maize_hapmap_missense_candidates.txt")
    cand = open(cand_path, "w")
    cand.write("protein\tgene\tchrom\tpos\tref\talt\tprot_pos\taa_change\tsnpeff_SUB\tsnpeff_TYPE\tvcf_MAF\talt_freq\tcall_rate\tn_called\treason\tselected\tmut_id\n")

    by_protein = defaultdict(list)
    for line in q.splitlines():
        f = line.split("\t")
        chrom, pos, ref, alt, vtype, gmodel, sub, vmaf = f[:8]
        gts = f[8:]
        pos = int(pos)
        if len(ref) != 1 or len(alt) != 1 or alt == ".":
            continue  # biallelic SNPs only
        for pid, m in models.items():
            if m["chrom"] != chrom or pos not in m["pos_map"]:
                continue
            # genotype-derived frequency
            n_alt = n_all = 0
            for g in gts:
                if g in ("./.", ".", "./"): continue
                al = g.replace("|", "/").split("/")
                for x in al:
                    if x == ".": continue
                    n_all += 1; n_alt += (x != "0")
            n_called = sum(1 for g in gts if not g.startswith("."))
            call_rate = n_called / len(gts)
            alt_freq = n_alt / n_all if n_all else float("nan")

            ci = m["pos_map"][pos]
            codon_i = ci // 3; off = ci % 3
            codon = m["cds"][codon_i*3:(codon_i+1)*3]
            strand_alt = alt if cds[pid]["strand"] == "+" else str(Seq(alt).complement())
            strand_ref = ref if cds[pid]["strand"] == "+" else str(Seq(ref).complement())
            assert codon[off] == strand_ref, f"REF mismatch at {chrom}:{pos} for {pid}"
            mut_codon = codon[:off] + strand_alt + codon[off+1:]
            aa_ref = str(Seq(codon).translate()); aa_alt = str(Seq(mut_codon).translate())
            prot_pos = codon_i + 1
            aa_change = f"{aa_ref}{prot_pos}{aa_alt}"

            if aa_ref == aa_alt: reason = "synonymous"
            elif aa_ref == "*": reason = "stop_lost"
            elif aa_alt == "*": reason = "nonsense"
            elif prot_pos in targets[pid]: reason = "overlaps_deleterious_site"
            elif call_rate < a.minCallRate: reason = "low_call_rate"
            elif alt_freq < a.minAltFreq: reason = "alt_freq_below_threshold"
            else: reason = "eligible"
            row = dict(protein=pid, gene=gene_of[pid], chrom=chrom, pos=pos, ref=ref, alt=alt, prot_pos=prot_pos,
                       aa_change=aa_change, sub=sub, vtype=vtype, vmaf=vmaf, alt_freq=alt_freq,
                       call_rate=call_rate, n_called=n_called, reason=reason)
            by_protein[pid].append(row)

    # ---- select neutral set ----
    selected = {}
    for pid in sorted(models):
        elig = [r for r in by_protein[pid] if r["reason"] == "eligible"]
        elig.sort(key=lambda r: (r["chrom"], r["pos"]))
        pick = random.sample(elig, min(a.nPerGene, len(elig)))
        pick.sort(key=lambda r: r["pos"])
        for i, r in enumerate(pick, 1):
            r["selected"] = True; r["mut_id"] = f"{pid}_neutralMut_{i}"
        selected[pid] = pick
        for r in by_protein[pid]:
            r.setdefault("selected", False); r.setdefault("mut_id", ".")
            cand.write("\t".join(str(r[k]) for k in ["protein","gene","chrom","pos","ref","alt","prot_pos","aa_change","sub","vtype","vmaf"])
                       + f"\t{r['alt_freq']:.4f}\t{r['call_rate']:.3f}\t{r['n_called']}\t{r['reason']}\t{int(r['selected'])}\t{r['mut_id']}\n")
    cand.close()

    # ---- final FASTA: WT, neutralMut_*, mut_* per gene ----
    del_fa = {r.id: str(r.seq) for r in SeqIO.parse(os.path.join(a.outdir, "maize_WT_deleterious.fa"), "fasta")}
    out_fa = os.path.join(a.outdir, "maize_trainingData.fa")
    n_neut = 0
    with open(out_fa, "w") as fa:
        for pid in sorted(models):
            wt = del_fa[pid]
            fa.write(f">{pid}\n{wt}\n")
            for r in selected[pid]:
                p = r["prot_pos"]; aa_ref, aa_alt = r["aa_change"][0], r["aa_change"][-1]
                assert wt[p-1] == aa_ref
                fa.write(f">{r['mut_id']}\n{wt[:p-1] + aa_alt + wt[p:]}\n"); n_neut += 1
            k = 1
            while f"{pid}_mut_{k}" in del_fa:
                fa.write(f">{pid}_mut_{k}\n{del_fa[f'{pid}_mut_{k}']}\n"); k += 1

    # ---- summary ----
    summ = os.path.join(a.outdir, "maize_neutral_summary.txt")
    with open(summ, "w") as s:
        s.write("protein\tn_snps_in_cds\tn_missense\tn_eligible\tn_selected\tn_deleterious\n")
        for pid in sorted(models):
            rows = by_protein[pid]
            n_mis = sum(1 for r in rows if r["reason"] not in ("synonymous", "nonsense", "stop_lost"))
            n_el = sum(1 for r in rows if r["reason"] == "eligible")
            s.write(f"{pid}\t{len(rows)}\t{n_mis}\t{n_el}\t{len(selected[pid])}\t{len(targets[pid])}\n")
    dist = defaultdict(int)
    for pid in models: dist[len(selected[pid])] += 1
    print(f"neutral mutants selected: {n_neut} across {len(models)} genes; genes by #selected: {dict(sorted(dist.items()))}")
    print(f"wrote {out_fa}\nwrote {cand_path}\nwrote {summ}")

if __name__ == "__main__":
    main()
