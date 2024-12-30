import re
import os
import sys
from Bio import SeqIO
import pandas as pd

# Define YTNAY regex for branch point scanning
branch_point_pattern = re.compile(r"[CT]T[ACGT]A[CT]")  # YTNAY

def extract_introns_and_analyze(fasta_file, gtf_file, output_dir):
    """
    Extract introns from GTF and FASTA, detect donor/acceptor sites and branchpoints.
    """
    # Containers for results
    introns_data = []
    sequences_5ss = []
    sequences_3ss = []

    # Read genome sequences
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # Group exons by transcript 
    exons_by_transcript = {}
    with open(gtf_file) as gtf:
    	for line in gtf:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            #print(fields)
            if fields[2] != "exon":
                continue

            #attributes = {kv.split(" ")[0]: kv.split(" ")[1].replace('"', '') for kv in fields[8].split(";") if kv.strip()}
            attributes = fields[8]
            
            gene_id=None
            transcript_id=None
            for attr in attributes.split(";"):
                if "gene_id" in attr:
                    gene_id = attr.split('"')[1]
                if "transcript_id" in attr:
                    transcript_id = attr.split('"')[1]
            
            print(gene_id)
            print(transcript_id)
            # Exon coordinates
            chrom = fields[0]
            strand = fields[6]
            exon_start = int(fields[3]) - 1  # Convert to 0-based
            exon_end = int(fields[4])  # End is 1-based

            # Add exon information
            if transcript_id not in exons_by_transcript:
                exons_by_transcript[transcript_id] = {
                    "gene_id": gene_id,
                    "chrom": chrom,
                    "strand": strand,
                    "exons": []
                }
            exons_by_transcript[transcript_id]["exons"].append((exon_start, exon_end))

    # Process each transcript
    for transcript_id, data in exons_by_transcript.items():
        gene_id = data["gene_id"]
        chrom = data["chrom"]
        strand = data["strand"]
        exons = sorted(data["exons"], key=lambda x: x[0])  # Keep genomic order

        # Assign intron numbers
        num_introns = len(exons) - 1
        if strand == "+":
            intron_number = list(range(1, num_introns + 1))  # [1, 2, ...]
        else:
            intron_number = list(range(num_introns, 0, -1))  # [n, n-1, ...]

        # Calculate introns
        for i in range(num_introns):
            prev_exon = exons[i]
            next_exon = exons[i + 1]

            intron_start = prev_exon[1]  # End of upstream exon
            intron_end = next_exon[0]   # Start of downstream exon

            # Extract intron sequence
            intron_seq = genome[chrom].seq[intron_start:intron_end] if strand == "+" else genome[chrom].seq[intron_start:intron_end].reverse_complement()
            intron_length = len(intron_seq)

            # Splice site sequences
            donor_site = intron_seq[:2]  # First 2 bases of intron
            acceptor_site_2 = intron_seq[-2:]  # Last 2 bases of intron
            acceptor_site_3 = intron_seq[-3:]  # Last 3 bases of intron

            # Branch point scanning (YTNAY)
            # Adjust intron sequence for branchpoint scanning
            branchpoint_seq = str(intron_seq[-600:]) if intron_length > 600 else str(intron_seq)

            # Branch point scanning (YTNAY)
            branchpoints = []
            distances_to_3ss = []
            for match in branch_point_pattern.finditer(branchpoint_seq):
                branchpoints.append(match.group())
                distances_to_3ss.append(len(branchpoint_seq) - match.end())
                
            # Add intron data to results
            introns_data.append({
                "transcript_id": transcript_id,
                "gene_id": gene_id,
                "chrom": chrom,
                "start": prev_exon[1],
                "end": next_exon[0],
                "strand": strand,
                "intron_number": intron_number[i],
                "intron_length": intron_length,
                "donor_site": str(donor_site),
                "acceptor_site_2": str(acceptor_site_2),
                "acceptor_site_3": str(acceptor_site_3),
                "branchpoints": ",".join(branchpoints),
                "distances_to_3ss": ",".join(map(str, distances_to_3ss))
            })

            # Prepare FASTA sequences for splice sites
            # 5'ss: donor site - last 3 bases of upstream exon, first 6 bases of intron
            if strand == "+":
                five_ss_seq = genome[chrom].seq[prev_exon[1] - 3:prev_exon[1]] + genome[chrom].seq[prev_exon[1]:prev_exon[1] + 6]
            else:
                # Next exon and reverse complement for negative strand
                five_ss_seq = genome[chrom].seq[next_exon[0] - 6:next_exon[0]] + genome[chrom].seq[next_exon[0]:next_exon[0] + 3]
                five_ss_seq = five_ss_seq.reverse_complement()

            sequences_5ss.append(f">{transcript_id}_intron_{intron_number[i]}\n{five_ss_seq}")

            # 3'ss: acceptor site - last 20 bases of intron, first 3 bases of downstream exon
            if strand == "+":
                three_ss_seq = genome[chrom].seq[next_exon[0] - 20:next_exon[0]] + genome[chrom].seq[next_exon[0]:next_exon[0] + 3]
            else:
                # Prev exon and reverse complement for negative strand
                three_ss_seq = genome[chrom].seq[prev_exon[1] - 3:prev_exon[1]] + genome[chrom].seq[prev_exon[1]:prev_exon[1] + 20]
                three_ss_seq = three_ss_seq.reverse_complement()

            sequences_3ss.append(f">{transcript_id}_intron_{intron_number[i]}\n{three_ss_seq}")

    # Save results to output directory
    os.makedirs(output_dir, exist_ok=True)
    print(introns_data)
    # Save table
    pd.DataFrame(introns_data).to_csv(os.path.join(output_dir, "introns_analysis.csv"), index=False)

    # Save FASTA files
    with open(os.path.join(output_dir, "5ss_sequences.fasta"), "w") as f:
        f.write("\n".join(sequences_5ss))

    with open(os.path.join(output_dir, "3ss_sequences.fasta"), "w") as f:
        f.write("\n".join(sequences_3ss))

# Parse command-line arguments
if len(sys.argv) != 4:
    print("Usage: python script.py <fasta_file> <gtf_file> <output_dir>")
    sys.exit(1)

fasta_file = sys.argv[1]
gtf_file = sys.argv[2]
output_dir = sys.argv[3]

# Run the script
extract_introns_and_analyze(fasta_file, gtf_file, output_dir)

