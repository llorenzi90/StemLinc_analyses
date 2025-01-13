#!/bin/bash

bed="/mnt/beegfs/public/references/annotation/mouse/GRCm39_GENCODE_VM31.bed"


wdir=$1 #/mnt/beegfs/llorenzi/jobs/LSK_RNA_seq
gbcov_out=$2 # EXP4_geneBodycov

cd $wdir

mdbams=$(ls analyses/*/*.MarkedDups.minq1.primarychrs.bam)
bams=$(ls analyses/*/*.minq1.primarychrs.bam)

# run infer experiments in markedDups bams
for bam in $mdbams

do	
	bn=${bam/Aligned.sortedByCoord.out.MarkedDups.minq1.primarychrs.bam/}
	#infer experiment
	sbatch -J $bn.infer_exp RSeQC_infer_experiment.sh $bam $bed
done	

# run read distribution in both markedDups and NoDups bams
for bam in $bams
do
	bn=${bam/Aligned.sortedByCoord.out./}
	bn=${bn/.minq1.primarychrs.bam/}
	#read distribution
	sbatch -J $bn.readDist RSeQC_read_distribution.sh $bam $bed
done


#gene body coverage
bamsIn_gbc=$(echo $bams|tr ' ' ,)
sbatch RSeQC_geneBody_coverage.sh $bamsIn_gbc $bed $gbcov_out


