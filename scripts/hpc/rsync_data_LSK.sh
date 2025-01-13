#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --job-name=rsync_data_to_jobs


datadir="/ijc/LABS/CUARTERO/RAW/StemLinc/20240219_StemLinc_LSK_RNA_seq/01.RawData/"

mkdir analyses

samps=$(ls $datadir)
for samp in $samps
do
	mkdir analyses/$samp
	rsync -rlv $datadir/$samp/*fq.gz analyses/$samp
done


md5sum analyses/*/*.fq.gz > md5sum_jobs.txt



