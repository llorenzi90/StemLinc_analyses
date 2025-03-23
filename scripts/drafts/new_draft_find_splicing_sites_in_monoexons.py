import re
import csv
from Bio import SeqIO
import os
from pybedtools import BedTool
from random import choice
import logging
from collections import Counter
import math

logging.basicConfig(level=logging.INFO)

# Splice site regex
DONOR_REGEX = re.compile(r"GT")
ACCEPTOR_REGEX = re.compile(r"AG")

# Define constants
MIN_INTRON_LENGTH = 20
GENOME_EXTENSION_FORMULA = lambda mono_length: max(1000, 2 * (1000 - mono_length))

# Define functions
def ensure_output_directory(output_prefix): # Ensure the output directory exists
    output_dir = os.path.dirname(output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

def extract_sequences(genome, chrom, start, end, strand):
    """Extract sequence from genome with coordinates."""
    if not (0 <= start < len(genome[chrom].seq)) or not (0 < end <= len(genome[chrom].seq)):
        raise ValueError(f"Invalid coordinates for {chrom}: start={start}, end={end}")

    sequence = genome[chrom].seq[start:end]
    
    if strand == "-":
        sequence = sequence.reverse_complement()
    return str(sequence)
   
def calculate_gc_content(sequence):
    """Calculate GC content of a sequence."""
    sequence = sequence.upper() 
    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count / len(sequence) if len(sequence) > 0 else 0

def calculate_shannon_entropy(sequence, k=1):
    """Calculate Shannon entropy for a sequence."""
    # Count k-mers in the sequence
    sequence = sequence.upper() 
    kmer_counts = Counter(sequence[i:i+k] for i in range(len(sequence) - k + 1))
    total_kmers = sum(kmer_counts.values())

    # Calculate frequencies
    frequencies = [count / total_kmers for count in kmer_counts.values()]

    # Calculate entropy
    entropy = -sum(p * math.log2(p) for p in frequencies) if frequencies else 0
    return entropy

def search_donor_sites(chrom, start, end, strand, genome):
    """Search for all donor (GT) sites in a sequence."""
    donors = []
    sequence = extract_sequences(genome, chrom, start, end, strand)
    
    match_number = 0
     
    for match in DONOR_REGEX.finditer(sequence):
        match_number += 1
        donor_pos = match.start()
        abs_donor_start = start + donor_pos if strand == "+" else end - donor_pos - 2
        abs_donor_end = abs_donor_start + 2
        donor_id = f"{chrom}:{abs_donor_start}-{abs_donor_end + 2}:{strand}"
        ext_donor_start = abs_donor_start - 2 if strand == "+" else abs_donor_start - 4 
        ext_donor_end = abs_donor_end + 4 if strand == "+" else abs_donor_end + 2 
        donor_seq = extract_sequences(genome, chrom, ext_donor_start, ext_donor_end, strand)
        donors.append((donor_id, donor_pos, donor_seq, match_number, abs_donor_start))
    return donors
    
def search_acceptor_sites(chrom, start, end, strand, genome):
    """Search for all acceptor (GT) sites in a sequence."""
    acceptors = []
    sequence = extract_sequences(genome, chrom, start, end, strand)
    
    match_number = 0
    
    for match in ACCEPTOR_REGEX.finditer(sequence):
        match_number += 1
        acceptor_pos = match.start()
        abs_acceptor_start = start + acceptor_pos if strand == "+" else end - acceptor_pos - 2
        abs_acceptor_end = abs_acceptor_start + 2
        acceptor_id = f"{chrom}:{abs_acceptor_start}-{abs_acceptor_end + 2}:{strand}"
        ext_acceptor_start = abs_acceptor_start - 12 if strand == "+" else abs_acceptor_start - 1 
        ext_acceptor_end = abs_acceptor_end + 1 if strand == "+" else abs_acceptor_end + 12 
        acceptor_seq = extract_sequences(genome, chrom, ext_acceptor_start, ext_acceptor_end, strand)
        acceptors.append((acceptor_id, acceptor_pos, acceptor_seq, match_number, abs_acceptor_start))
    return acceptors

def find_potential_introns(donors, acceptors):
    """Find valid introns based on donor and acceptor positions."""
    potential_introns = []
    
    match_number = 0
    
    for donor_id, donor_pos, _, _, _ in donors:
        for acceptor_id, acceptor_pos, _, _, _ in acceptors:
            if acceptor_pos > donor_pos + MIN_INTRON_LENGTH:  # Minimum intron length
                match_number += 1
                intron_length = acceptor_pos - donor_pos
                potential_introns.append((donor_id, acceptor_id, donor_pos, acceptor_pos, match_number, intron_length))
    return potential_introns

def extract_5_splice_site_sequences(genome, chrom, donor_start, strand):
    """Extract sequences for 5' splice sites specifically for FASTA."""
    # 5'ss: donor site - last 3 bases of upstream exon, first 6 bases of intron
    if strand == "+":
        donor_seq = genome[chrom].seq[donor_start - 3:donor_start + 6]  # 5' splice site: 3 bases upstream + 6 bases downstream
        
    else:
        donor_seq = genome[chrom].seq[donor_start - 4:donor_start + 5].reverse_complement()  # Reverse complement for negative strand
        
    return str(donor_seq)
    
def extract_3_splice_site_sequences(genome, chrom, acceptor_start, strand):
    """Extract sequences for 3' splice sites specifically for FASTA."""
    # 3'ss: acceptor site - last 20 bases of intron, first 3 bases of downstream exon
    if strand == "+":
        acceptor_seq = genome[chrom].seq[acceptor_start - 18:acceptor_start + 5]  # 3' splice site: 20 bases upstream + 3 bases downstream
    
    else: 
        acceptor_seq = genome[chrom].seq[acceptor_start - 3:acceptor_start + 20].reverse_complement() # Reverse complement for negative strand
    
    return str(acceptor_seq)

def save_fasta(sequences, output_file):
    """Save sequences to a fasta file."""
    with open(output_file, "w") as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n{seq}\n")

def find_random_region(gtf_file, seq_id, chrom, monoexon_start, monoexon_end, strand, region_size, chromosome_sizes, MIN_DISTANCE=10000, DISTANCE_INCREMENT=500):
    """Find a random region avoiding overlaps with GTF annotations."""
    annotated_regions = BedTool(gtf_file).filter(lambda x: x[0] == chrom).saveas()
    distance = MIN_DISTANCE
    chrom_size = chromosome_sizes[chrom]
    
    # Define direction-specific calculations once, outside the loop
    my_directions = {
        "upstream": lambda dist: (
            max(0, monoexon_start - dist - region_size),
            max(0, monoexon_start - dist - region_size) + region_size,
        ),
        "downstream": lambda dist: (
            min(chrom_size, monoexon_end + dist + region_size) - region_size,
            min(chrom_size, monoexon_end + dist + region_size),
        ),
    }

    while True:
        # Randomly select a direction and compute the other
        rand_dir = choice(list(my_directions.keys()))
        other_dir = "downstream" if rand_dir == "upstream" else "upstream"

        # Attempt with the randomly selected direction
        random_start, random_end = my_directions[rand_dir](distance)
        random_region= BedTool(f"{chrom}\t{random_start}\t{random_end}\t{seq_id}_random\t.\t{strand}\n", from_string=True)

        # Check if the region is valid
        if random_region.intersect(annotated_regions, u=True).count() == 0:
            return random_region

        # If fallback needed, try the other direction
        random_start, random_end = my_directions[other_dir](distance)
        random_region= BedTool(f"{chrom}\t{random_start}\t{random_end}\t{seq_id}_random\t.\t{strand}\n", from_string=True)

        if random_region.intersect(annotated_regions, u=True).count() == 0:
            return random_region

        # Increment distance if neither works
        distance += DISTANCE_INCREMENT


def process_sequences(sequences, genome, output_prefix):
    """Process sequences and generate outputs for donors, acceptors, and introns."""
    fasta_5ss = {}
    fasta_3ss = {}
    
    donor_file = open(f"{output_prefix}_donor_sites.tsv", "w", newline="")
    acceptor_file = open(f"{output_prefix}_acceptor_sites.tsv", "w", newline="")
    intron_file = open(f"{output_prefix}_potential_introns.tsv", "w", newline="")
    sequence_composition_file = open(f"{output_prefix}_sequence_composition.tsv", "w", newline="")

    donor_writer = csv.writer(donor_file, delimiter="\t")
    acceptor_writer = csv.writer(acceptor_file, delimiter="\t")
    intron_writer = csv.writer(intron_file, delimiter="\t")
    sequence_composition_writer = csv.writer(sequence_composition_file, delimiter="\t")
    
    donor_writer.writerow(["seq_ID", "chr","start","end", "extended_start","extended_end", "strand", "donor_site_id", "donor_number", "donor_pos", "extended_donor_seq"])
    acceptor_writer.writerow(["seq_ID", "chr","start","end", "extended_start","extended_end", "strand", "acceptor_site_id", "acceptor_number", "acceptor_pos", "extended_acceptor_seq"])
    intron_writer.writerow(["seq_ID", "donor_site_id", "acceptor_site_id", "donor_pos", "acceptor_pos", "potential_intron_number", "potential_intron_length"])
    sequence_composition_writer.writerow(["seq_ID", "GC_content","GC_content_extended", "Shannon_entropy_1", "Shannon_entropy_1_extended", "Shannon_entropy_3","Shannon_entropy_3_extended"])
    
    for seq_id, chrom, start, end, extended_start, extended_end, strand in sequences:
        # Search for donor and acceptor sites
        donors = search_donor_sites(chrom, extended_start, extended_end, strand, genome)
        acceptors = search_acceptor_sites(chrom, extended_start, extended_end, strand, genome)
        introns = find_potential_introns(donors, acceptors)
        sequence = extract_sequences(genome, chrom, start, end, strand)
        sequence_extended = extract_sequences(genome, chrom, extended_start, extended_end, strand)
        
        # Calculate Shannon entropy and GC content
        GC_content = calculate_gc_content(sequence)
        Shannon_entropy_1 = calculate_shannon_entropy(sequence, k=1)
        Shannon_entropy_3 = calculate_shannon_entropy(sequence, k=3)
        GC_content_extended = calculate_gc_content(sequence_extended)
        Shannon_entropy_1_extended = calculate_shannon_entropy(sequence_extended, k=1)
        Shannon_entropy_3_extended = calculate_shannon_entropy(sequence_extended, k=3)
        
        # Write sequence composition info
        logging.info("Writing sequence composition info to %s", f"{output_prefix}_sequence_composition.tsv")
        sequence_composition_writer.writerow([seq_id, GC_content,GC_content_extended, Shannon_entropy_1, Shannon_entropy_1_extended, Shannon_entropy_3,Shannon_entropy_3_extended])
        
        # Write donor sites
        for donor_id, donor_pos, donor_seq, donor_number, abs_donor_start in donors:
            logging.info("Writing donor sites to %s", f"{output_prefix}_donor_sites.tsv")
            donor_writer.writerow([seq_id, chrom, start, end, extended_start, extended_end, strand, donor_id, donor_number, donor_pos, donor_seq])
            if donor_id not in fasta_5ss:
                fasta_5ss[donor_id] = extract_5_splice_site_sequences(genome, chrom, abs_donor_start, strand)  # Specific 5' splice site sequence

        # Write acceptor sites
        for acceptor_id, acceptor_pos, acceptor_seq, acceptor_number, abs_acceptor_start in acceptors:
            logging.info("Writing acceptor sites to %s", f"{output_prefix}_acceptor_sites.tsv")
            acceptor_writer.writerow([seq_id, chrom, start, end, extended_start, extended_end, strand, acceptor_id, acceptor_number, acceptor_pos, acceptor_seq])
            if acceptor_id not in fasta_3ss:
                fasta_3ss[acceptor_id] = extract_3_splice_site_sequences(genome, chrom, abs_acceptor_start, strand)  # Specific 3' splice site sequence

        # Write potential introns
        for donor_id, acceptor_id, donor_pos, acceptor_pos, match_number, intron_length in introns:
            intron_writer.writerow([seq_id, donor_id, acceptor_id, donor_pos, acceptor_pos, match_number, intron_length])

    donor_file.close()
    acceptor_file.close()
    intron_file.close()

    # Save FASTA files
    save_fasta(fasta_5ss, f"{output_prefix}_5ss.fasta")
    save_fasta(fasta_3ss, f"{output_prefix}_3ss.fasta")

# Main processing function
def analyze_monoexons_and_random(monoexons, genome, output_prefix, gtf_file):
    """Analyze monoexons and random regions."""
    random_sequences = []
    monoexon_sequences = []
    chromosome_sizes = {chrom: len(seq_record.seq) for chrom, seq_record in genome.items()}
    fasta_random = {}
    
    for monoexon in monoexons:
        chrom, start, end, seq_id, _, strand = monoexon
        
        logging.info("Processing sequence: %s", seq_id)  # Log the current sequence being processed

        start, end = int(start), int(end)
        mono_length = end - start

        # Extend downstream region for monoexons
        extension = GENOME_EXTENSION_FORMULA(mono_length)
        extended_start = start if strand == "+" else start - extension
        extended_end = end + extension if strand == "+" else end
        
        # Calculate extended regions size
        extended_length = extended_end - extended_start

        # Find random region
        logging.info("Finding random region for sequence: %s", seq_id)
        random_region = find_random_region(gtf_file, seq_id, chrom, extended_start, extended_end, strand, extended_length, chromosome_sizes)
        _, ran_start_extended, ran_end_extended, ran_seq_id, _, _ = random_region[0].fields
        ran_start_extended, ran_end_extended = int(ran_start_extended), int(ran_end_extended)
        ran_start = ran_start_extended if strand == "+" else ran_start_extended + extension
        ran_end = ran_end_extended - extension if strand == "+" else ran_end_extended
        
        # Add random sequence to fasta
        fasta_random[ran_seq_id] = extract_sequences(genome, chrom, ran_start_extended, ran_end_extended, strand)
        
        # Add monoexon and random sequence to lists
        monoexon_sequences.append((seq_id, chrom, start, end, extended_start, extended_end, strand))
        random_sequences.append((ran_seq_id, chrom, ran_start, ran_end, ran_start_extended, ran_end_extended, strand ))
    
    # Write fasta of random regions
    save_fasta(fasta_random, f"{output_prefix}_random.fasta")
    # Process both monoexons and random sequences
    process_sequences(monoexon_sequences + random_sequences, genome, output_prefix)

#

if __name__ == "__main__":
    # Input files
    genome_fasta = "/home/llorenzi/Documentos/test_data/fasta/chr1.mm39.fasta"  # Replace with actual genome FASTA file
    gtf_file = "/home/llorenzi/Documentos/test_data/gtf/LSK_StemLinc.combined.sorted.chr1.gtf"  # Replace with actual GTF file
    #monoexon_file = "/home/llorenzi/Rprojects/StemLinc_analyses/outputs/bed_files/tr_level/LSK_StemLinc.combined.sorted.chr1_monoexons.20_random.bed"  # Replace with actual monoexon BED file
    monoexon_file = "/home/llorenzi/Rprojects/StemLinc_analyses/outputs/bed_files/tr_level/LSK_StemLinc.combined.sorted.chr1_monoexons.100_random.bed"  # Replace with actual monoexon BED file
    output_prefix = "/home/llorenzi/Documentos/test_outputs/monoexon_analysis_test"  # Replace with your desired output file prefix
    
    
    if not os.path.exists(genome_fasta):
        raise FileNotFoundError(f"Genome file not found: {genome_fasta}")

    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file not found: {gtf_file}")

    if not os.path.exists(monoexon_file):
        raise FileNotFoundError(f"Monoexon BED file not found: {monoexon_file}")

    # Ensure output directory exists
    ensure_output_directory(output_prefix)
    
    # Load genome
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    
    # Load monoexons
    monoexons = BedTool(monoexon_file)
    
    # Run the analysis
    analyze_monoexons_and_random(monoexons, genome, output_prefix, gtf_file)


