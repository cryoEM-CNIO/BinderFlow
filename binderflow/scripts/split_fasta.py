import argparse
from Bio import SeqIO
import os

def create_output_directory(output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

def extract_sequence_name(sequence):
    return sequence.id.split()[0]  # Take the first part of the sequence ID as the name

def add_prefix_to_sequence(sequence, prefix):
    return f"{prefix}{sequence.seq}"

def main(input_file, output_dir, prefix):
    create_output_directory(output_dir)

    # Read the input file
    sequences = list(SeqIO.parse(input_file, "fasta"))

    # Create a new .fasta file for each sequence
    for sequence in sequences:
        sequence_name = extract_sequence_name(sequence)
        output_file = os.path.join(output_dir, f"{sequence_name}.fasta")

        # Add the prefix to the sequence
        modified_sequence = add_prefix_to_sequence(sequence, prefix)

        # Write the modified sequence to the output file
        with open(output_file, "w") as output_handle:
            output_handle.write(f">{sequence_name}\n{modified_sequence}\n")

        print(f"Sequence '{sequence_name}' written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split a multi-sequence FASTA file into individual files and add a prefix.")
    parser.add_argument("input_file", help="Path to the input FASTA file")
    parser.add_argument("output_dir", help="Path to the output directory")
    parser.add_argument("--prefix", default="ASDF", help="Prefix to add to the beginning of each sequence")
    args = parser.parse_args()

    main(args.input_file, args.output_dir, args.prefix)
