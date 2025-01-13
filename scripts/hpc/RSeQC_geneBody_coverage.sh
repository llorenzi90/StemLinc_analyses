#!/bin/bash
#SBATCH -o %x.%j.geneBody_cov.out
#SBATCH --nodes=1
#SBATCH --time=25:00:00
##SBATCH --cpus-per-task=10
#SBATCH --mem=50GB
#SBATCH --job-name=gene_body_cov


module load RSeQC/4.0.0-foss-2021b



bamsIn_gbc=$1
bed=$2
out=$3

geneBody_coverage.py -r $bed -i $bamsIn_gbc -o $out
