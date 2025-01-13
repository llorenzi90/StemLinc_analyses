#!/bin/bash
#SBATCH -o %x.%j.filterbam_rmDups.out
#SBATCH --nodes=1
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=100GB
##SBATCH --job-name=filterbam_rmDups

PPN=$SLURM_CPUS_PER_TASK
primary_chrs_file=/mnt/beegfs/llorenzi/jobs/StemLinc/GRCm39_primary_chrs.bed

#cat $primary_chrs_file
#chr1	0	195154279
#chr2	0	181755017
#chr3	0	159745316
#chr4	0	156860686
#chr5	0	151758149
#chr6	0	149588044
#chr7	0	144995196
#chr8	0	130127694
#chr9	0	124359700
#chr10	0	130530862
#chr11	0	121973369
#chr12	0	120092757
#chr13	0	120883175
#chr14	0	125139656
#chr15	0	104073951
#chr16	0	98008968
#chr17	0	95294699
#chr18	0	90720763
#chr19	0	61420004
#chrX	0	169476592
#chrY	0	91455967
#chrM	0	16299



sdir=$1 #each sample has its own directory, for example /mnt/beegfs/llorenzi/jobs/StemLinc/Experiments2-3/analyses/LT-HSC_R2  
sample=$(basename $sdir)

module load SAMtools/1.13-foss-2021b
module load picard/2.26.3-Java-11

cd $sdir

inbam=${sample}Aligned.sortedByCoord.out.bam 

# Step 1: mark duplicates with picard mark duplicates
	
input_mark=$inbam
output_mark=${input_mark/.bam/.MarkedDups.bam}

java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$input_mark O=$output_mark M=$input_mark.dupmatrix CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false TAGGING_POLICY=All

# Step 2: filter out non-primary chrs and MAPQ smaller than 1

input_st=$output_mark
output_st=${input_st/.bam/.minq1.primarychrs.bam}

samtools view -@ $PPN -L $primary_chrs_file -q 1 -b -h -o $output_st $input_st


# Step 3: remove duplicates with samtools

input_rd=$output_st
output_rd=${input_rd/MarkedDups/NoDups}

samtools view -@ $PPN -F 1024 -b -h -o $output_rd $input_rd



# Step 4: run samtools index stats for all bam files

for bam in *.bam
do

	samtools index -@ $PPN $bam 	
	samtools stats -@ $PPN $bam > $bam.stats
done

