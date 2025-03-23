from Bio import SeqIO
import os
import sys

def extract_splicing_sites(fasta_file, exon_bed_file, output_file=None):
    """
    Extract acceptor and donor splice sites from an input BED file of exons.

    Parameters:
    fasta_file (str): Path to the genome FASTA file.
    exon_bed_file (str): Path to the input BED file with exons.
    output_file (str): Optional path to save the output TSV file.
                       If not provided, writes to <input_base_name>.splicing_sites.tsv in the current directory.
    """
    # Determine output file name if not provided
    if output_file is None:
        base_name = os.path.basename(exon_bed_file)
        output_file = f"{os.getcwd()}/{base_name.rsplit('.', 1)[0]}.splicing_sites.tsv"

    # Load genome sequences
    print("Loading genome sequences...")
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # Process the BED file
    print(f"Processing exons from {exon_bed_file}...")
    with open(exon_bed_file) as bed, open(output_file, "w") as out:
        # Write header
        out.write("chrom\tstart\tend\texonID\tstrand\tdonor_site\tacceptor_site\n")
        
        for line in bed:
            if not line.strip():  # Skip empty lines
                continue
            
            # Parse BED columns
            chrom, start, end, exonID, _, strand = line.strip().split("\t")
            start, end = int(start), int(end)

            # Initialize default values
            donor_site, acceptor_site = "NA", "NA"

            # Ensure chromosome exists in the FASTA file
            if chrom in genome:
                seq = genome[chrom].seq
                seq_len = len(seq)

                if strand == "+":
                    # Donor site: 2 bases after the exon end
                    if end + 2 <= seq_len:  # Check that slicing is within bounds
                        donor_site = str(seq[end:end + 2])
                    # Acceptor site: 2 bases before the exon start
                    if start - 2 >= 0:  # Check that slicing is within bounds
                        acceptor_site = str(seq[start - 2:start])

                elif strand == "-":
                    # Donor site: 2 bases before the exon start (reverse complement)
                    if start - 2 >= 0:  # Check that slicing is within bounds
                        donor_site = str(seq[start - 2:start].reverse_complement())
                    # Acceptor site: 2 bases after the exon end (reverse complement)
                    if end + 2 <= seq_len:  # Check that slicing is within bounds
                        acceptor_site = str(seq[end:end + 2].reverse_complement())

            # Write results to the output file
            out.write(f"{chrom}\t{start}\t{end}\t{exonID}\t{strand}\t{donor_site}\t{acceptor_site}\n")

    print(f"Splicing sites saved to {output_file}")


if __name__ == "__main__":
    # Check arguments
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python extract_splicing_sites.py <genome.fasta> <input_exons.bed> [output_file.tsv]")
        sys.exit(1)

    # Parse arguments
    fasta_path = sys.argv[1]
    bed_path = sys.argv[2]
    output_path = sys.argv[3] if len(sys.argv) == 4 else None

    # Run the extraction
    extract_splicing_sites(fasta_path, bed_path, output_path)

