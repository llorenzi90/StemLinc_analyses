#!/bin/bash

#fasta=outputs/fasta/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.promoters.fa
fasta=$1
echo "Input fasta:"
echo "$fasta"
outname=$(basename $fasta)
outname=${outname/.fa/}


/home/llorenzi/software/meme-5.5.6/src/fasta-get-markov -m 5 $fasta > outputs/bg_models/$outname.5th_order_background_model.txt 
