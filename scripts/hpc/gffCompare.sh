#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --nodes=1
#SBATCH --time=06:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=gffcompare_LSK

module load GffCompare/0.12.6-GCC-11.2.0

refgtf="/mnt/beegfs/llorenzi/jobs/StemLinc/gffcompare_merge_ref_annotations/merged_refs.combined.gtf"
smaskedgenome="/mnt/beegfs/public/references/genome/mouse/soft_masked/mm39.fa"

query_gtfs_file=$1 # strtie_guided_gtfs.txt
outpref=$2 # LSK_StemLinc
outdir=$3 # gffcompare

mkdir -p $outdir

# gffCompare with various GTF queries
gffcompare -i $query_gtfs_file -r $refgtf -s $smaskedgenome -V -j $outdir/$outpref.new_junction.tsv -o $outdir/$outpref

# gffCompare with the combined GTF against ref GTF
gffcompare -r $refgtf -s $smaskedgenome -V -j $outdir/combined_$outpref.new_junction.tsv -o $outdir/combined_$outpref $outdir/$outpref.combined.gtf
