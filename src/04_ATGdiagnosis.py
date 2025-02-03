import pandas as pd
import sys

def parse_fasta(file_path):
    """Parse a FASTA file and return a list of sequences."""
    sequences = []
    with open(file_path, 'r') as file:
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                sequence = ''
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence)
    return sequences

def find_atg_sites_case_insensitive(sequences):
    """Find 'ATG' sites in the reference sequence and count matches in other sequences."""
    reference = sequences[-1].lower()
    atg_sites = {}
    for i in range(len(reference) - 2):
        if reference[i:i+3] == 'atg':
            atg_sites[i] = sum(seq[i:i+3].lower() == 'atg' for seq in sequences[:-1])
    return atg_sites

def correct_ungapped_coordinates(reference, gapped_coordinates):
    """Calculate the ungapped coordinates in the reference sequence."""
    return [pos - reference[:pos].count('-') for pos in gapped_coordinates]

def main(input_file, output_file):
    sequences = parse_fasta(input_file)
    
    # Analysis
    atg_sites = find_atg_sites_case_insensitive(sequences)
    ungapped_coords = correct_ungapped_coordinates(sequences[-1], list(atg_sites.keys()))

    # Preparing data for output
    data = {
        'Reference Coordinate': list(atg_sites.keys()),
        'Number of Sequences': list(atg_sites.values()),
        'Ungapped Coordinate': ungapped_coords
    }

    df = pd.DataFrame(data)
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Output saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_fasta_file> <output_file>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        main(input_file, output_file)
