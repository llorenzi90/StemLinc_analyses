in_bed="outputs/bed_files/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.promoters.bed"

bedtools_get_fasta.sh $in_bed


in_bed="outputs/bed_files/LSK_StemLinc.combined.introns.bed"

bedtools_get_fasta.sh $in_bed
