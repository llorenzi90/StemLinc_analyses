#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --nodes=1
#SBATCH --time=06:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=gffcompare_combined_LSK

module load GffCompare/0.12.6-GCC-11.2.0

refgtf="/mnt/beegfs/llorenzi/jobs/StemLinc/gffcompare_merge_ref_annotations/merged_refs.combined.gtf"
smaskedgenome="/mnt/beegfs/public/references/genome/mouse/soft_masked/mm39.fa"

outpref=$1

gffcompare -r $refgtf -s $smaskedgenome -V -j gffcompare/combined_$outpref.new_junction.tsv -o gffcompare/combined_$outpref gffcompare/$outpref.combined.gtf


