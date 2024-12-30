awk '{for(i=1;i<=NF;i++) if ($i == ".") $i = "1000"; print}' outputs/bed_files/LSK_StemLinc.combined.sorted.exons.bed > tmp
mv tmp outputs/bed_files/LSK_StemLinc.combined.sorted.exons.bed
Rscript scripts/bigWigAverageOverBed_phastCons35way.R outputs/bed_files/LSK_StemLinc.combined.sorted.exons.bed



Rscript scripts/bigWigAverageOverBed_phastCons35way.R outputs/bed_files/LSK_StemLinc.combined.sorted.introns.bed

