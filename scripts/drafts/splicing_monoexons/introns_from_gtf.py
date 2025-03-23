import sys  # 1
import os  # 2
from collections import defaultdict  # 3

def parse_gtf(gtf_file):  # 5
    """
    Parse a GTF file and extract exon coordinates grouped by transcript ID.
    """  # 7
    exons = defaultdict(list)  # 8
    with open(gtf_file) as gtf:  # 9
        for line in gtf:  # 10
            if line.startswith("#") or not line.strip():  # 11
                continue  # Skip comments and empty lines # 12
            
            fields = line.strip().split("\t")  # 13
            if fields[2] != "exon":  # 14
                continue  # Only process exon entries # 15
            
            chrom = fields[0]  # 16
            start = int(fields[3]) - 1  # Convert to 0-based # 17
            end = int(fields[4])  # GTF end is already exclusive # 18
            strand = fields[6]  # 19
            
            attributes = fields[8]  # 20
            transcript_id = None  # 21
            for attr in attributes.split(";"):  # 22
                if "transcript_id" in attr:  # 23
                    transcript_id = attr.split('"')[1]  # 24
                    break  # 25
            
            if transcript_id:  # 26
                exons[transcript_id].append((chrom, start, end, strand))  # 27
    return exons  # 28


def find_introns(exons_by_transcript):  # 30
    """
    Find intronic regions from exon coordinates and track their transcript IDs.
    """  # 32
    introns = []  # 33
    monoexonic_count = 0  # Counter for monoexonic transcripts # 34
    for transcript_id, exons in exons_by_transcript.items():  # 35
        if len(exons) == 1:  # Check if the transcript is monoexonic # 36
            print(f"Skipping monoexonic transcript: {transcript_id}")  # Log monoexonic transcript # 37
            monoexonic_count += 1  # Increment counter # 38
            continue  # Skip processing for monoexonic transcripts # 39
	# Assign intron numbers
        num_introns = len(exons) - 1
        strand = exons[0][3]
        if strand == "+":
            intron_numbers = list(range(1, num_introns + 1))  # [1, 2, ...]
        else:
            intron_numbers = list(range(num_introns, 0, -1))  # [n, n-1, ...]
            
        exons = sorted(exons, key=lambda x: x[1])  # Sort exons by start # 40
        for i in range(num_introns):  # 41
            chrom, end, next_start, strand,  intron_number = exons[i][0], exons[i][2], exons[i + 1][1], exons[i][3], intron_numbers[i]  # 42
            combined_id = f"{transcript_id}_{intron_number}" 
            introns.append((chrom, end, next_start, strand, combined_id))  # Add combined_id # 43

    print(f"Total monoexonic transcripts skipped: {monoexonic_count}")  # Final log # 44
    return introns  # 45


def write_introns_to_bed(introns, output_file):  # 47
    """
    Write intronic regions to a BED file, including transcript IDs.
    """  # 49
    with open(output_file, "w") as out:  # 50
        for chrom, start, end, strand, transcript_id in introns:  # Include transcript_id # 51
            out.write(f"{chrom}\t{start}\t{end}\t{transcript_id}\t1000\t{strand}\n")  # Include transcript_id in BED # 52


if __name__ == "__main__":  # 54
    # Check arguments
    if len(sys.argv) < 2 or len(sys.argv) > 3:  # 56
        print("Usage: python extract_introns_from_gtf.py <input.gtf> [output.bed]")  # 57
        sys.exit(1)  # 58

    input_gtf = sys.argv[1]  # 60

    # Determine the output file name
    if len(sys.argv) == 3:  # 62
        output_bed = sys.argv[2]  # User-specified output path
    else:  # Auto-generate the output name
        base_name = os.path.basename(input_gtf)  # Extract basename of input GTF
        output_bed = f"{os.getcwd()}/{base_name.rsplit('.', 1)[0]}.introns.bed"  # Replace extension

    exons_by_transcript = parse_gtf(input_gtf)  # 66
    introns = find_introns(exons_by_transcript)  # 67
    write_introns_to_bed(introns, output_bed)  # 68
    print(f"Intronic regions saved to {output_bed}")  # 69

# cd ~/Rprojects/StemLinc_analyses
# python scripts/introns_from_gtf.py data/raw/LSK_StemLinc.combined.sorted.gtf outputs/bed_files/LSK_StemLinc.combined.sorted.introns.bed 

