import argparse
from Bio import SeqIO

def find_premature_stops(input_file, output_file):
    with open(output_file, 'w') as out:
        out.write("Sequence_ID\tPremature_Stop\tStop_Position\n")  # Header for the output file
        for record in SeqIO.parse(input_file, "fasta"):
            sequence = str(record.seq)
            stop_pos = sequence.find('*')
            if stop_pos != -1 and stop_pos != len(sequence) - 1:
                out.write(f"{record.id}\t1\t{stop_pos}\n")
            else:
                out.write(f"{record.id}\t0\t{len(sequence) + 1}\n")

def main():
    parser = argparse.ArgumentParser(description="Identify premature stop codons in protein sequences.")
    parser.add_argument("input_file", type=str, help="Path to the input FASTA file containing protein sequences.")
    parser.add_argument("output_file", type=str, help="Path to the output file where results will be saved.")
    args = parser.parse_args()

    find_premature_stops(args.input_file, args.output_file)

if __name__ == "__main__":
    main()

