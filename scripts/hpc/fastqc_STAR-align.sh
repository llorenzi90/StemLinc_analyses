#!/bin/bash
#SBATCH --job-name=generate_STAR_index
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=10
#SBATCH --time=15:35:30
#SBATCH --mem-per-cpu=10G

PPN=$SLURM_CPUS_PER_TASK

STARidx=/mnt/beegfs/public/references/index/GRCm39_STAR

load_STAR="module load STAR/2.7.6a-GCC-11.2.0"
load_fastqc="module load FastQC/0.11.9"

sdir=$1 #each sample has its own directory, for example /mnt/beegfs/llorenzi/analyses/StemLinc/Experiments2-3/analyses/LT-HSCs_R2  

sample=$(basename $sdir) 

cd $sdir

read1=$(ls *_1.fq.gz)
read2=$(ls *_2.fq.gz)
        
        
if ! [ -f "$read1" ]; then
                die "input file $read1 not found"
fi


if ! [ -f "$read2" ]; then
                die "input file $read2 not found"
fi
echo "read1: $read1"
echo "read2: $read2"

echo `date`
echo -e "################################\n\tFASTQC started\n################################\n"
$load_fastqc
mkdir fastqc
fastqc $read1 -o fastqc
fastqc $read2 -o fastqc
echo -e "################################\n\tFASTQC done\n################################\n"
#echo `date  +'%r'`

echo `date`
echo -e "################################\n\tSTAR started\n################################\n"
$load_STAR
STAR --runThreadN $PPN --readFilesIn $read1 $read2 --readFilesCommand zcat --genomeDir $STARidx --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $sample --outSAMattributes NH HI AS nM NM MD jM jI XS MC ch
echo -e "################################\n\tSTAR done\n################################\n"
echo `date  +'%r'`


