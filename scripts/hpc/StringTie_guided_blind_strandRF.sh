#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --nodes=1
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=50GB
#SBATCH --job-name=StringTie

module load StringTie/2.1.7-GCC-11.2.0

PPN=$SLURM_CPUS_PER_TASK
bam=$1

ref_gtf=$2

#option -M:
# -M fraction of bundle allowed to be covered by multi-hit reads (default:1)
#I will set it at 0.3 


#guided 
stringtie -p $PPN -v -G $ref_gtf -M 0.3 --rf -C $bam.cov_refs.gtf -o $bam.strtie.guided.gtf $bam
	
#blind
stringtie -p $PPN -v -M 0.3 --rf -o $bam.strtie.blind.gtf $bam  


