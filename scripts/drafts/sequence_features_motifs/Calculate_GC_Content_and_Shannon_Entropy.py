import csv
from collections import Counter
import math

def calculate_gc_content(sequence):
    """Calculate GC content of a sequence."""
    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count / len(sequence) if len(sequence) > 0 else 0

def calculate_shannon_entropy(sequence, k=1):
    """Calculate Shannon entropy for a sequence."""
    # Count k-mers in the sequence
    kmer_counts = Counter(sequence[i:i+k] for i in range(len(sequence) - k + 1))
    total_kmers = sum(kmer_counts.values())

    # Calculate probabilities
    probabilities = [count / total_kmers for count in kmer_counts.values()]

    # Calculate entropy
    entropy = -sum(p * math.log2(p) for p in probabilities) if probabilities else 0
    return entropy

def analyze_sequences(sequences, genome, output_file):
    """
    Analyze sequences for GC content and Shannon entropy.
    Save results to a CSV file.
    """
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Sequence_ID", "Type", "GC_Content", "Shannon_Entropy"])

        for seq_id, chrom, start, end, _, _, strand in sequences:
            # Extract sequence
            sequence = extract_sequences(genome, chrom, start, end, strand)

            # Calculate metrics
            gc_content = calculate_gc_content(sequence)
            complexity = calculate_shannon_entropy(sequence)

            # Write to CSV
            writer.writerow([seq_id, "monoexon" if "mono" in seq_id else "random", gc_content, complexity])

    print(f"Analysis complete. Results saved to {output_file}")

# Example Usage
if __name__ == "__main__":
    # Define your input sequences and genome
    genome_fasta = "/path/to/genome.fasta"
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    # Example input sequences (replace with real monoexons and random sequences)
    monoexons = [
        ("mono_1", "chr1", 1000, 2000, 1000, 2000, "+"),
        ("mono_2", "chr1", 3000, 4000, 3000, 4000, "+"),
    ]
    random_sequences = [
        ("random_1", "chr1", 5000, 6000, 5000, 6000, "+"),
        ("random_2", "chr1", 7000, 8000, 7000, 8000, "+"),
    ]

    # Combine monoexons and random sequences
    all_sequences = monoexons + random_sequences

    # Output file
    output_file = "sequence_analysis_results.csv"

    # Run analysis
    analyze_sequences(all_sequences, genome, output_file)

