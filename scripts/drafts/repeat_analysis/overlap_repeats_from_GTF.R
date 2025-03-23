############################################
##
## Purpose of script: Compute overlap of genes with repeats from RepeatMasker
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-12-19
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
args <- commandArgs(trailingOnly = TRUE)

# Check if an input file is provided
if (length(args) == 0) {
  stop("Please provide the input files as a command-line argument.")
}

gtf_path <- args[1] # "data/raw/LSK_StemLinc.combined.sorted.gtf"
GTF <- readGFF(gtf_path)
# annotation of the reference transcriptome used, a merge between GENCODE and RefSeq

GTF$gene_id <- GTF$gene_name

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



get_overlap_with_repeats_from_GTF=function(gtf){

  pgme=trastools::get_per_gene_merged_exons(gtf)
  pgme=pgme[,c(3:5,2,6,7)]
  pgme=pgme[order(pgme$seqnames,pgme$start),]
  exonic_length=pgme%>%group_by(group_name)%>%
    summarise(exonic_length=sum(width))
  pgme$start <- pgme$start - 1 # 0-based for bed format
  ol=bt.intersect(a = pgme,b=rmsk_bed,wo = T)
  ol=ol%>% rowwise()%>% mutate(ol_st=max(c(V2,V8)),ol_en=min(c(V3,V9)))
  colnames(ol) <- c("chr", "start","end","gene_name","merged_exon_width","strand",
                    "rep_chr","rep_start","rep_end","rep_name","score","rep_strand",
                    "rep_class","rep_family","ol_len","ol_start","ol_end")
  olsm=ol%>%group_by(gene_name)%>%
    summarise(repeats=paste(unique(rep_name[order(-ol_len)]),
                            collapse = ","),
              repeats_classes=paste(unique(rep_class[order(-ol_len)]),
                                    collapse = ","),
              t.overlap_length=sum(IRanges::width(IRanges::reduce(IRanges(
                start = ol_start,end = ol_end)))))
  olsm=left_join(olsm,exonic_length,by=c(gene_name="group_name"))
  olsm$repeat.fraction=olsm$t.overlap_length/olsm$exonic_length

  # strand-specific
  olsm_stranded <- ol%>%group_by(gene_name,rep_strand)%>%
    summarise(repeats=paste(unique(rep_name[order(-ol_len)]),
                            collapse = ","),
              repeats_classes=paste(unique(rep_class[order(-ol_len)]),
                                    collapse = ","),
              t.overlap_length=sum(IRanges::width(IRanges::reduce(IRanges(
                start = ol_start,end = ol_end)))))

  olsm_stranded=left_join(olsm_stranded,exonic_length,by=c(gene_name="group_name"))
  olsm_stranded$repeat.fraction=olsm_stranded$t.overlap_length/olsm_stranded$exonic_length


  return(list(olsm,olsm_stranded,ol, exonic_length))
}


#

# calculate repeat overlap ----
ol_repeats=get_overlap_with_repeats_from_GTF(GTF)

# calculate the total fraction of bases covered by repeats ----
# per biotype


exonic_length <- ol_repeats[[4]]
exonic_length=left_join(exonic_length,
                        gene_level_info%>%dplyr::select(gene_name,
                                                        biotype),
                        by=c(group_name="gene_name"))

olsm=ol_repeats[[1]]


# stranded
olsm_stranded=ol_repeats[[2]]

olsm_stranded$strand <- GTF$strand[match(olsm_stranded$gene_name,GTF$gene_name)]
olsm_stranded <- mutate(olsm_stranded, orientation = ifelse(rep_strand==strand,"sense","antisense"))

olsm_stranded_summ <- pivot_wider(olsm_stranded, id_cols = gene_name,names_from = orientation,
                                  values_from = repeat.fraction, names_prefix = "repeat.fraction_")
## write tables in db ----

### all orientations ----
head(olsm)

gene_overlap_RepeatMasker_any_orientation <- olsm %>% dplyr::select(gene_name,repeats,
                                                                    repeats_classes,total_overlap_length = t.overlap_length,
                                                                    total_exonic_length = exonic_length,
                                                                    fraction_covered_by_repeats = repeat.fraction)



### sense orientation ----

gene_overlap_RepeatMasker_sense_orientation <-  olsm_stranded %>% filter(orientation=="sense") %>%
  dplyr::select(gene_name,gene_strand = strand,  sense_repeats = repeats,sense_repeats_classes = repeats_classes,
                total_sense_overlap_length = t.overlap_length,
                total_exonic_length = exonic_length, fraction_covered_by_sense_repeats = repeat.fraction)


### antisense orientation ----

gene_overlap_RepeatMasker_antisense_orientation <-  olsm_stranded %>% filter(orientation=="antisense") %>%
  dplyr::select(gene_name,gene_strand = strand,  antisense_repeats = repeats,antisense_repeats_classes = repeats_classes,
                total_antisense_overlap_length = t.overlap_length,
                total_exonic_length = exonic_length, fraction_covered_by_antisense_repeats = repeat.fraction)





