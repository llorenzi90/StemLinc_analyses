options(scipen = 999)
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
require(rtracklayer)
require(bedtoolsr)
require(trastools)

chr_order=paste0("chr",c(seq(1:22),"M","X","Y"))

## Load data---------------------------
source("scripts/source_all_functions.R")
#dir.create(outdir)

gene_level_data_path <- "outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.biotypes_of_interest.last.tsv"

check_file_header(gene_level_data_path)
gene_level_info <- read.table(gene_level_data_path,header = T)
#gene_level_bed <- extract_bed(gene_level_info)

# transcript_level_data_path <- "outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.gene_classif.tsv"
# transcript_level_data <- read.table(transcript_level_data_path,header = T)
# transcript_level_bed <- extract_bed(transcript_level_data,name = "V1",score = "gene_name")
# TSS_bed <- transcript_level_bed %>% mutate(end=start+1)
# TSS_bed <- TSS_bed %>% dplyr::filter(!duplicated(TSS_bed%>%dplyr::select(seqid,start,end,strand)))

gl_out_pref=gsub(".tsv","",basename(gene_level_data_path))
#tl_out_pref=gsub(".tsv","",basename(transcript_level_data_path))


biots_of_interest=c("potNovel","lncRNA","TEC","pseudogene","protein_coding")
ncbiots=c("lncRNA","TEC","potNovel","pseudogene")
lncRNA_potNovel_TEC=c("lncRNA","TEC","potNovel")
potNovel_TEC=c("TEC","potNovel")
biotypes2select=list(biots_of_interest,
                     ncbiots,
                     lncRNA_potNovel_TEC,
                     potNovel_TEC)

classes=sort(unique(gene_level_info$best_classif_to_PCG))

simpl_classes=c("antisense","convergent",
                "convergent","divergent","divergent",
                "intergenic","intronic_antisense",
                "intronic_sense","sense","sense","sense","sense",
                "sense")
