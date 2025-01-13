#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --nodes=1
#SBATCH --time=25:00:00
#SBATCH --mem=50GB
#SBATCH --job-name=get_bigWigs

sdir=$1 #each sample has its own directory, for example /mnt/beegfs/llorenzi/jobs/StemLinc/Experiments2-3/analyses/LT-HSC_R2  


cd $sdir

# load required module
module load deepTools/3.3.1-foss-2021b-Python-3.8.5

##Generate bigWig files

for input in *minq1.primarychrs.bam
do

	bamCoverage --binSize 50 --filterRNAstrand reverse --normalizeUsing None --bam $input -o $input.minus.bw
	bamCoverage --binSize 50 --filterRNAstrand forward --normalizeUsing None --bam $input -o $input.plus.bw

	bamCoverage --samFlagExclude 256 --binSize 50 --filterRNAstrand reverse --normalizeUsing None --bam $input -o $input.minus.primaryalns.bw
	bamCoverage --samFlagExclude 256 --binSize 50 --filterRNAstrand forward --normalizeUsing None --bam $input -o $input.plus.primaryalns.bw

done
#256 SAM flag: not primary alignment
