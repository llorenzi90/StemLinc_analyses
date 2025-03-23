############################################
##
## Purpose of script: Remove repetitive regions from XLOC_079903
##
## Author: Lucia Lorenzi
##
## Date Created: 2025-01-29
##
## Email: lucialorenzi90@gmail.com
##
#  Notes ---------------------------
##
##
##
#  Setup ---------------------------

options(scipen = 999)
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
library(rtracklayer)
library(trastools)
require(bedtoolsr)
library(DBI)
library(RSQLite)
library(AnnotationHub)
#repeats=read.table("/home/llorenzi/work_local/UCSC_RepeatMasker_mm39.bed")
db <- dbConnect(RSQLite::SQLite(), "outputs/dbs/StemLinc.db")
GTF <- readGFF("outputs/gtf/LSK_StemLinc.combined.filtered.20250127_112131.gene_name.gtf")
# annotation of the reference transcriptome used, a merge between GENCODE and RefSeq

transcript_level_info <- read.delim("outputs/transcriptome_characterization/jan25/LSK_StemLinc.combined_annotated_tracking.filtered.20250127_112131.tsv")

gene_level_info <-read.delim("outputs/transcriptome_characterization/jan25/LSK_StemLinc.combined_annotated_tracking.filtered.20250127_112131.gene_level.tsv")




# Initialize AnnotationHub
ah <- AnnotationHub()

# Search for RepeatMasker data for the mm39 genome
query_mm39 <- query(ah, c("mm39", "RepeatMasker"))

# List available resources
query_mm39

# Retrieve the RepeatMasker data
rmsk_mm39 <- query_mm39[["AH99013"]]

# Inspect the GRanges object
rmsk_mm39

# Convert RepeatMasker GRanges to a BED-like data frame
rmsk_bed <- data.frame(
  chrom = as.character(seqnames(rmsk_mm39)),     # Chromosome
  start = start(rmsk_mm39) - 1,                 # 0-based start
  end = end(rmsk_mm39),                         # End
  name = rmsk_mm39$repName,                     # Repeat name
  score = 0,                                    # Default score
  strand = as.character(strand(rmsk_mm39)),     # Strand
  repClass = rmsk_mm39$repClass,                # Repeat class
  repFamily = rmsk_mm39$repFamily               # Repeat family
)

# Replace NA strands with "."
rmsk_bed$strand[is.na(rmsk_bed$strand)] <- "."


# generate exon bed
exon_bed <- GTF %>% filter(type=="exon") %>% mutate(start=start-1,score=1000) %>% select(seqid,start,end,transcript_id,score,strand)

subtracted_reps <- bt.subtract(exon_bed,rmsk_bed)

write.table(subtracted_reps,"outputs/bed_files/LSK_StemLinc.combined_annotated_tracking.filtered.20250127_112131.subtractedReps.bed",
            quote = F,col.names = F,row.names = F,sep = "\t")

# write in GTF format
colnames(GTF)
subtracted_reps$type="exon"
subtracted_reps$gene_id=GTF$gene_name[match(subtracted_reps$V4,GTF$transcript_id)]
subtracted_reps$source="StringTie_subtracted_reps"
colnames(subtracted_reps) <- c("seqid","start","end","transcript_id","score","strand","type","gene_id","source")

subtracted_reps$start=subtracted_reps$start +1

subtracted_reps <- subtracted_reps[order(subtracted_reps$seqid,
                                         subtracted_reps$start),]
export(subtracted_reps,"outputs/gtf/LSK_StemLinc.combined.filtered.20250127_112131.gene_name.subtracted_reps.gtf",format = "gtf")
