#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --nodes=1
#SBATCH --time=1:00:00
##SBATCH --cpus-per-task=10
#SBATCH --mem=24GB
#SBATCH --job-name=infer_experiment


module load RSeQC/4.0.0-foss-2021b

bam=$1
bed=$2


infer_experiment.py -r $bed -i $bam > $bam.infer_experiment.txt

