#!/bin/bash
# run gffCompare 

gtf1=$1
gtf2=$2
out_pref=$3

ref_gtf=${4:-'data/references/merged_refs_annotation/merged_refs.combined.annot_trnames.gtf'}


# check hpc runs
gffcompare -r $ref_gtf -o $out_pref $gtf1 $gtf2

