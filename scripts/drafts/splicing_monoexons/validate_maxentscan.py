import sys
import os
from Bio import SeqIO

def validate_maxentscan_sequences(input_fasta, expected_length):
    """
    Validates sequences in a FASTA file for MaxEntScan.
    
    INPUT:
        input_fasta: str, path to the input FASTA file
        expected_length: int, the expected length of sequences
    
    OUTPUT:
        A cleaned FASTA file with only valid sequences.
    """
    # Define valid bases
    valid_bases = {"A", "C", "G", "T"}
    
    # Create the output file path
    base, ext = os.path.splitext(input_fasta)
    output_fasta = f"{base}.clean{ext}"

    # Open the input and output files
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            # Check length and valid characters
            seq = str(record.seq).upper()
            if len(seq) == expected_length and set(seq).issubset(valid_bases):
                # Write valid sequences to the output file
                SeqIO.write(record, outfile, "fasta")

    print(f"Cleaned sequences written to: {output_fasta}")

# Example usage
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_fasta> <expected_length>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    try:
        expected_length = int(sys.argv[2])
    except ValueError:
        print("Error: expected_length must be an integer.")
        sys.exit(1)
    
    validate_maxentscan_sequences(input_fasta, expected_length)

