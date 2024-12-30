#!/bin/bash

# Script to extract unique exon coordinates from a GTF file and output in BED format

# Check for correct number of arguments
if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 <input.gtf> [output.bed]"
    exit 1
fi

# Input GTF file
input_gtf=$1

# Determine output file name
if [ "$#" -eq 2 ]; then
    output_bed=$2
else
    # Replace .gtf with .exons.bed in the input filename, and write to the current directory
    base_name=$(basename "$input_gtf")            # Extract filename from path
    output_bed="${PWD}/${base_name%.gtf}.exons.bed"  # Replace .gtf with .exons.bed
fi

# Extract exons and format as BED with unique exon IDs
awk '
    $3 == "exon" {
        exonID = $1 ":" ($4-1) "-" $5 ":" $7;
        print $1 "\t" ($4-1) "\t" $5 "\t" exonID "\t.\t" $7;
    }
' "$input_gtf" | sort -k1,1 -k2,2n -k3,3n -u > "$output_bed"

# Completion message
echo "Unique exons saved to $output_bed"

