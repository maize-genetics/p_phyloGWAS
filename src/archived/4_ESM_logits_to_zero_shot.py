# Archived: earlier draft of ESM2 zero-shot score conversion, superseded by the complete,
# actually-wired implementation, src/10_logit2zeroShot.R (called from
# notebook/07_summaryStats/07Ca_ESM_SCINET.sh). Has a real bug - `resDic['name'] = df` uses
# the literal string 'name' instead of the loop variable, so every iteration overwrites the
# same dict key and only the last sequence's scores survive np.savez_compressed(). Kept for
# reference only; do not use.
import numpy as np
import pandas as pd
from Bio import SeqIO
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-inputFASTA", dest="inputFASTA", type=str, default=None,
                        help="The directory of input fasta")
    parser.add_argument("-inputLogit", dest="inputLogit", type=str, default=None, help="The directory of input logit")
    parser.add_argument("-output", dest = "output", help = "The directory of output")
    args = parser.parse_args()
    return args


def read_fasta(inputFASTA):
    sequences = []
    names = []
    for record in SeqIO.parse(inputFASTA, "fasta"):
        sequences.append(str(record.seq))
        names.append(record.id)
    return sequences, names


def calculate_log_ratio(row):
    """
    Calculate the log ratio of the reference amino acid probability to the largest probability of the remaining amino acids
    """
    ref_col = row['ref']
    ref_prob = row[ref_col]
    
    remaining_probs = row.drop(['ref', ref_col])
    largest_prob = remaining_probs.max()
    
    log_ratio = np.log(ref_prob / largest_prob)
    return log_ratio

args = parse_args()

sequences, names = read_fasta(args.inputFASTA)
logits = np.load(args.inputLogit)

amino_acids = list('LAGVSERTIDPKQNFYMHWC')
resDic = {}
for sequence, name in zip(sequences, names):
    logit = logits[name]
    df = pd.DataFrame(logit, columns = amino_acids)
    df['ref'] = list(sequence)
    df['zero_shot'] = df.apply(calculate_log_ratio, axis=1)
    df = df.drop('ref', axis=1)
    resDic['name'] = df

np.savez_compressed(args.output, **resDic)
    
