#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --nodes=1
#SBATCH --time=05:00:00
##SBATCH --cpus-per-task=10
#SBATCH --mem=50GB
#SBATCH --job-name=read_distribution

module load RSeQC/4.0.0-foss-2021b

bam=$1

bed=$2

read_distribution.py -r $bed -i $bam > $bam.readDist
