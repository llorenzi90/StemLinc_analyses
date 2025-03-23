#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --job-name=featureCounts_oct24

PPN=$SLURM_CPUS_PER_TASK

##load required modules
module load Subread/2.0.3-GCC-11.2.0

#go to working dir

wdir=/mnt/beegfs/llorenzi/jobs/LSK_RNA_seq

cd $wdir


bams=$(cat all_public_and_SL_bams_oct2024.txt|tr "\n" " ")
strands=$(cat all_public_and_SL_bams.featureCount_strand.txt|tr "\n" ",")

## Quantify LSK potential novel genes with featureCounts
gtf=/ijc/USERS/llorenzi/gtfs/LSK_StemLinc.combined.20240930_173216.filtered.gene_name.gtf
outname=LSK_StemLinc.combined.20240930_173216.filtered.gene_name.all_blood_cells.primary.txt

featureCounts -T $PPN -s $strands --primary --verbose -O --fracOverlap 0.5 -p -M -B --countReadPairs -a $gtf -o $outname $bams


outname=LSK_StemLinc.combined.20240930_173216.filtered.gene_name.all_blood_cells.txt

featureCounts -T $PPN -s $strands --verbose -O --fracOverlap 0.5 -p -M -B --countReadPairs -a $gtf -o $outname $bams



#####Description of input parameters:

# -s <int or string>  Perform strand-specific read counting. A single integer
# 	value (applied to all input files) or a string of comma-
#   	separated values (applied to each corresponding input
#                     file) should be provided. Possible values include:
#   	0 (unstranded), 1 (stranded) and 2 (reversely stranded).
# 	Default value is 0 (ie. unstranded read counting carried
#                     out for all input files).


# --primary           Count primary alignments only. Primary alignments are 
# identified using bit 0x100 in SAM/BAM FLAG field.
#


# -O                  Assign reads to all their overlapping meta-features (or 
#                      features if -f is specified).


# --fracOverlap <float> Minimum fraction of overlapping bases in a read that is
# required for read assignment. Value should be within range
# [0,1]. 0 by default. Number of overlapping bases is
# counted from both reads if paired end. Both this option
# and '--minOverlap' option need to be satisfied for read
# assignment.


# -p                  Specify that input data contain paired-end reads. To
# perform fragment counting (ie. counting read pairs), the
# '--countReadPairs' parameter should also be specified in
# addition to this parameter.
 

# -M                  Multi-mapping reads will also be counted. For a multi-
#                     mapping read, all its reported alignments will be 
#                     counted. The 'NH' tag in BAM/SAM input is used to detect 
#                     multi-mapping reads.



# -B                  Only count read pairs that have both ends aligned.
  

# --countReadPairs    Count read pairs (fragments) instead of reads. This option
# is only applicable for paired-end reads.





