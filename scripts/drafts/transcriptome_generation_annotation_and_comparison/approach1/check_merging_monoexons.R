options(scipen = 999)
## Load packages---------------------------
require(tidyverse)
require(data.table)
require(rtracklayer)
library(trastools)
args <- commandArgs(trailingOnly = TRUE)

# Check if an input dir is provided
if (length(args) == 0) {
   stop("Please provide the working directory as a command-line argument.")
}

# Setup ----
dir_path <- args[1] # path to wdir
setwd(dir_path) # dir_path="outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation_1"

tracking_path <- "~/Documentos/all_samples_analyses/gffcompare_replicates.vs.VM36/StemLinc_replicates.vs.VM36.tracking"
gtf_path <- "~/Documentos/all_samples_analyses/gffcompare_replicates.vs.VM36/StemLinc_replicates.vs.VM36.combined.gtf"
tracking <- read.table(tracking_path)
gtf <- readGFF(gtf_path)

gtf_new <- gtf %>%filter(transcript_id%in%tracking$V1[!tracking$V4%in%trastools::overlapping_class_codes])
export(gtf_new,gsub(".gtf",".NOTolapRef.gtf",gtf_path),format = "gtf")
# Extract transcript features:
N_exons <- as.integer(trastools::get_n_exons_per_transcript(tracking))
N_libraries <- trastools::get_n_samples_per_transcript_from_tracking(tracking)
TPMs <- trastools::get_expression_values(tracking)
table(N_exons,N_libraries)

gtf_XLOC_527626 <- gtf%>%filter(gene_id=="XLOC_527626",type=="transcript")

library(bedtoolsr)
source("scripts/functions/get_promoter_regions_gtf2bed_sort_bed_merge_beds_intersects_etc.R")
XLOC_bed=extract_bed(gtf_XLOC_527626,name = "transcript_id")
merge=bt.merge(XLOC_bed,s=T,c = "6",o="distinct")

bt_cov=bedtoolsr::bt.coverage(a =merge,b=XLOC_bed,d = T)
barplot(bt_cov$V6)
write.table(merge,"~/Documentos/test_bedtools/merged_XLOC_527626.bed",col.names = F,row.names = F,quote = F,sep = "\t")
write.table(XLOC_bed,"~/Documentos/test_bedtools/XLOC_527626.bed",col.names = F,row.names = F,quote = F,sep = "\t")

bt_cov_filtered <- bt_cov %>% filter(V6>=4)
bt_cov_filtered_coords <- bt_cov_filtered %>%mutate(start=V2+V5,end=start+1,score=V6) %>%
  dplyr::select(V1,start,end,V5,score,V4)
filtered_merge=bt.merge(bt_cov_filtered_coords,c = "6",o="distinct")
write.table(filtered_merge,"~/Documentos/test_bedtools/filtered_merged_XLOC_527626.bed",col.names = F,row.names = F,quote = F,sep = "\t")

# filtering by TPMs
XLOC_bed_tpm=TPMs$mean_tpm[match(XLOC_bed$name,tracking$V1)]
XLOC_bed_filtered <- XLOC_bed[XLOC_bed_tpm>=0.2,]
write.table(XLOC_bed_filtered,"~/Documentos/test_bedtools/XLOC_527626_filtered.bed",col.names = F,row.names = F,quote = F,sep = "\t")

XLOC_filtered_merged <- bt.merge(XLOC_bed_filtered,s=T,c = "6",o="distinct")
write.table(XLOC_filtered_merged,"~/Documentos/test_bedtools/XLOC_527626_filtered_merged.bed",col.names = F,row.names = F,quote = F,sep = "\t")
bt_cov2=bedtoolsr::bt.coverage(a =XLOC_filtered_merged,b=XLOC_bed_filtered,d = T)
bt_cov2_filtered <- bt_cov2 %>% filter(V6>=4)
bt_cov2_filtered_coords <- bt_cov2_filtered %>%mutate(start=V2+V5,end=start+1,score=V6) %>%
  dplyr::select(V1,start,end,V5,score,V4)
bt_cov2_filtered_merge=bt.merge(bt_cov2_filtered_coords,c = "6",o="distinct")
write.table(bt_cov2_filtered_merge,"~/Documentos/test_bedtools/test2filtered_merged_XLOC_527626.bed",col.names = F,row.names = F,quote = F,sep = "\t")

DGEA <- read.delim("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation_1/DGEA/filtered_genes.CPAT.CAGE.DGEA.tsv")
