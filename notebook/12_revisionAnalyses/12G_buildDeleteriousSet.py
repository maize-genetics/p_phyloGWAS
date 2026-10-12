"""Build WT + deleterious-mutant protein FASTA for the maize set, mirroring the
At layout (>ID for WT, >ID_mut_N for each deleterious missense).

Usage:
  python src/maize_01_build_deleterious.py \
      -table data/deleteriousMissenseMutationsMaize.txt \
      -peptide data/B73v5/Zm-B73-REFERENCE-NAM-5.0_peptide.fa \
      -out output/maize/maize_WT_deleterious.fa \
      -report output/maize/maize_deleterious_validation.txt

The table lists protein IDs (Zm..._P001); the peptide FASTA is keyed by
transcript ID (Zm..._T001). P### <-> T### is the B73v5 convention.
"""
import argparse, os
from Bio import SeqIO

# manual corrections for entries that do not match the reference protein
# (position 139 of Zm00001eb252380_P001 is G; the A is at 140 - treated as an
# off-by-one in the source table, per discussion)
CORRECTIONS = {("Zm00001eb252380_P001", "A139T"): "A140T"}

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-table", required=True)
    p.add_argument("-peptide", required=True)
    p.add_argument("-out", required=True)
    p.add_argument("-report", required=True)
    return p.parse_args()

def main():
    a = parse_args()
    pep = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(a.peptide, "fasta")}

    muts = []  # (protein, original_mut, used_mut)
    with open(a.table) as fh:
        for line in fh:
            if line.startswith("Protein") or not line.strip():
                continue
            prot, mut = line.rstrip("\n").split("\t")
            used = CORRECTIONS.get((prot, mut), mut)
            muts.append((prot, mut, used))

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fa = open(a.out, "w")
    rep = open(a.report, "w")
    rep.write("protein\ttranscript\tmutation_in_table\tmutation_used\tprotein_len\tref_in_fasta\tstatus\tmut_id\n")

    written_wt = set()
    counter = {}
    n_ok = n_bad = 0
    for prot, orig, used in muts:
        tid = prot.replace("_P", "_T")
        seq = pep[tid]
        ref, pos, alt = used[0], int(used[1:-1]), used[-1]
        if pos > len(seq) or seq[pos - 1] != ref:
            found = seq[pos - 1] if pos <= len(seq) else "OUT_OF_RANGE"
            rep.write(f"{prot}\t{tid}\t{orig}\t{used}\t{len(seq)}\t{found}\tMISMATCH_SKIPPED\t.\n")
            n_bad += 1
            continue
        if prot not in written_wt:
            fa.write(f">{prot}\n{seq}\n")
            written_wt.add(prot)
        counter[prot] = counter.get(prot, 0) + 1
        mut_id = f"{prot}_mut_{counter[prot]}"
        mutseq = seq[:pos - 1] + alt + seq[pos:]
        fa.write(f">{mut_id}\n{mutseq}\n")
        status = "OK" if orig == used else f"OK_CORRECTED_FROM_{orig}"
        rep.write(f"{prot}\t{tid}\t{orig}\t{used}\t{len(seq)}\t{ref}\t{status}\t{mut_id}\n")
        n_ok += 1

    fa.close(); rep.close()
    print(f"genes (WT written): {len(written_wt)}; deleterious mutants written: {n_ok}; skipped: {n_bad}")
    print(f"wrote {a.out}\nwrote {a.report}")

if __name__ == "__main__":
    main()
