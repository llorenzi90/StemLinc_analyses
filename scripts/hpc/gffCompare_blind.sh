#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --nodes=1
#SBATCH --time=06:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=gffcompare_LSK

module load GffCompare/0.12.6-GCC-11.2.0

refgtf="/mnt/beegfs/llorenzi/jobs/StemLinc/gffcompare_merge_ref_annotations/merged_refs.combined.gtf"
smaskedgenome="/mnt/beegfs/public/references/genome/mouse/soft_masked/mm39.fa"

outpref=$1 #LSK_StemLinc.blind

mkdir gffcompare_blind

gffcompare -i strtie_blind_gtfs.txt -r $refgtf -s $smaskedgenome -V -j gffcompare_blind/$outpref.new_junction.tsv -o gffcompare_blind/$outpref


gffcompare -r $refgtf -s $smaskedgenome -V -j gffcompare_blind/combined_$outpref.new_junction.tsv -o gffcompare_blind/combined_$outpref gffcompare_blind/$outpref.combined.gtf

