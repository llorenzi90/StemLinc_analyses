
## Setup ---------------------------

options(scipen = 999)
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
require(rtracklayer)
library(trastools)
## Load data---------------------------
#setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/")
#setwd("public_datasets/metadata_RNA_seq_studies/Delas_et_al/")
#source("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/scripts/StemLinc_code/get_genomic_regions_function.R")
delas_genes=fread("data/public_data/lncRNA_studies/Delas_et_al/Delas_all_lncRNAs_2017_from_website.txt")
delas_2017gtfmm39=readGFF("data/public_data/lncRNA_studies/Delas_et_al/elife-25607-supp1-v2_Delas_lncRNA_catalog.liftovermm39.gtf")
delas_genes_mm39=get_genomic_range_by_gene(gtf_df = delas_2017gtfmm39)
table(unlist(strsplit(delas_genes$lnc,split = ";"))%in%delas_genes_mm39$gene_id)

delas_genes_mm39$coord=paste0(delas_genes_mm39$seqnames,
                              ":",
                              delas_genes_mm39$start,
                              "-",
                              delas_genes_mm39$end)
delas_genes=delas_genes[rep(1:nrow(delas_genes),sapply(strsplit(delas_genes$lnc,split = ";"),length)),]
delas_genes$lnc=unlist(strsplit(delas_genes$lnc[!duplicated(delas_genes$lnc)],split = ";"))
colnames(delas_genes)[2]="mm10coord"

delas_genes_mm39 <- left_join(delas_genes_mm39,delas_genes,by=c(gene_id="lnc"))
delas_genes_mm39$lnc_name=NA
delas_genes_mm39$lnc_name[delas_genes_mm39$gene_id=="XLOC_166788"]="Sphed"
delas_genes_mm39$lnc_name[delas_genes_mm39$gene_id=="XLOC_104449"]="Pilna"
delas_genes_mm39$lnc_name[delas_genes_mm39$gene_id=="XLOC_177417"]="Lilam"

#perform overlap with gffcompare

#gffcompare -r '/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/metadata_RNA_seq_studies/Delas_et_al/elife-25607-supp1-v2_Delas_lncRNA_catalog.liftovermm39.gtf' -V -o StemLincpotNovLSK_vs_Delas '/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare/LSK_StemLinc.combined.gtf'

tracking=read.table("data/raw/LSK_StemLinc.tracking")
refmap=read.table('/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare/StemLincpotNovLSK_vs_Delas/StemLincpotNovLSK_vs_Delas.LSK_StemLinc.combined.gtf.refmap',
                  header = T)

length(unique(refmap$ref_id))
length(unique(delas_2017gtfmm39$transcript_id))

length(unique(refmap$ref_gene))/length(unique(delas_2017gtfmm39$gene_id))*100
# 53 % of the Delas genes have partial or full match with StemLinc assemblies

delas_genes_mm39_not_in_SL=delas_genes_mm39 %>% filter(!gene_id%in%refmap$ref_gene)

# check overlap of delas with other public data and perform venn between the 3 datasets

colnames(refmap)[1]="gene_id"

matched_Delas_in_StemLinc <- left_join(refmap, delas_genes_mm39 )

source("scripts/load_gene_level_data_biotypes_of_interest.R")
transcript_level_data <- read.table("outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.gene_classif.tsv",header = T)

length(unique(matched_Delas_in_StemLinc$gene_id))

library(trastools)
matched_Delas_in_StemLinc <- split_samples_info(tracking = matched_Delas_in_StemLinc,cols =4,qnames = "StemLinc",remove = F )

matched_Delas_in_StemLinc$StemLinc_transcript=sapply(strsplit(matched_Delas_in_StemLinc$StemLinc_transcript,split = ","),function(x)x[1])
matched_Delas_in_StemLinc$StemLinc_gene_name=transcript_level_data$gene_name[match(
  matched_Delas_in_StemLinc$StemLinc_transcript, transcript_level_data$V1
)]

#tracking_delas=read.table('/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare/StemLincpotNovLSK_vs_Delas/StemLincpotNovLSK_vs_Delas.tracking')
tmap_delas=read.table('/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare/StemLincpotNovLSK_vs_Delas/StemLincpotNovLSK_vs_Delas.LSK_StemLinc.combined.gtf.tmap',header = T)

table(delas_genes_mm39$gene_id%in%tmap_delas$ref_gene_id)
tmap_delas$gene_name=transcript_level_data$gene_name[match(tmap_delas$qry_id,
                                                           transcript_level_data$V1)]
tmap_delas <- tmap_delas%>%filter(!is.na(gene_name))
tmap_delas_overlap=tmap_delas%>%filter(class_code%in%overlapping_class_codes)
length(unique(tmap_delas_overlap$ref_gene_id))

tmap_delas_overlap <- left_join(tmap_delas_overlap,delas_genes_mm39,by=c(ref_gene_id="gene_id"))

tmap_delas_overlap_summary <- tmap_delas_overlap %>% group_by(gene_name) %>%
  arrange(match(tmap_delas_overlap$class_code,overlapping_class_codes)) %>%
  summarise(Delas_best_gene=ref_gene_id[1],
            Delas_best_cc=class_code[1],
            Delas_genes=paste(unique(ref_gene_id),collapse = ","),
            Delas_ccs=paste(unique(class_code),collapse = ","))

tmap_delas_overlap_summary <- left_join(tmap_delas_overlap_summary,gene_level_info %>% select(gene_name, biotype))

tmap_delas_overlap_summary <- left_join(tmap_delas_overlap_summary,
                                        delas_genes_mm39 %>%select(gene_id,lnc_name,orientation,mm10coord),by=c(Delas_best_gene="gene_id"))

dir.create("outputs/overlap_lncRNA_catalogs")

write.table(tmap_delas_overlap_summary,"outputs/overlap_lncRNA_catalogs/overlap_StemLinc_genes_Delas_genes.txt")
