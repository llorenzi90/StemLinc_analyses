#!/bin/bash

# run from project home

#gtf1=$1
#gtf2=$2
#out_pref=$3

mkdir outputs/gffcompare
gtf1="data/raw/LSK_StemLinc.combined.gtf"
gtf2="data/raw/Kli_MarkedDups_guided.combined.gtf"
out_pref="outputs/gffcompare/SLvsKli_LSK"

./scripts/run_gffCompare_compare_two.sh $gtf1 $gtf2 $out_pref # check where output is

#input_R1 = ""
#input_R2 = ""
# move some files to new folders if needed in this line
#./Rscript compare_two_transcriptomes.R $input_R1 $input_R2
