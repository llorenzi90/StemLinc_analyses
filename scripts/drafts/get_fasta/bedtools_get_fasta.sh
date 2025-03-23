#!/bin/bash

genome_fasta="data/references/genomes/mm39.fa"

#in_bed="outputs/bed_files/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.promoters.bed"
in_bed=$1

echo "Input bed:"
echo "$in_bed"

out_dir="outputs/fasta"
out_file=$(basename $in_bed)
out_file=${out_file/bed/fa}
out_fasta=$out_dir/$out_file

bedtools getfasta -fi $genome_fasta -fo $out_fasta -bed $in_bed -s -name

