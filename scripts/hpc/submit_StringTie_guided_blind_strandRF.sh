#!/bin/bash
gtf=/mnt/beegfs/public/references/annotation/mouse/gencode.vM31.primary_assembly.annotation.gtf

wdir=$1 #/mnt/beegfs/llorenzi/jobs/LSK_RNA_seq

cd $wdir

bams=$(ls analyses/*/*.MarkedDups.minq1.primarychrs.bam)

for bam in $bams

do 
	bn=$(basename $bam)
	bn=${bn/Aligned.sortedByCoord.out.MarkedDups.minq1.primarychrs.bam/}
	sbatch -J $bn.strtie StringTie_guided_blind_strandRF.sh $bam $gtf
done

