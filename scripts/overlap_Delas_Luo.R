## Setup ---------------------------

options(scipen = 999)
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
require(rtracklayer)
library(trastools)
source("scripts/source_all_functions.R")
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

# read Stem Linc transcript level data ----
transcript_level_data <- read.table("outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.gene_classif.tsv",header = T)

# read tmap from gffcomparison StemLinc vs Delas (Delas used as reference) ----
#tracking_delas=read.table('/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare/StemLincpotNovLSK_vs_Delas/StemLincpotNovLSK_vs_Delas.tracking')
tmap_delas=read.table("data/public_data/lncRNA_studies/StemLincpotNovLSK_vs_Delas.LSK_StemLinc.combined.gtf.tmap")

table(delas_genes_mm39$gene_id%in%tmap_delas$ref_gene_id)

# add gene name from the transcript level data ----
tmap_delas$gene_name=transcript_level_data$gene_name[match(tmap_delas$qry_id,
                                                           transcript_level_data$V1)]


tmap_delas <- tmap_delas%>%filter(!is.na(gene_name))
tmap_delas_overlap=tmap_delas%>%filter(class_code%in%overlapping_class_codes)
length(unique(tmap_delas_overlap$ref_gene_id))

# add info on delas genes
tmap_delas_overlap <- left_join(tmap_delas_overlap,delas_genes_mm39,by=c(ref_gene_id="gene_id"))

# make summary per StemLinc gene name
tmap_delas_overlap_summary <- tmap_delas_overlap %>% group_by(gene_name) %>%
  arrange(match(tmap_delas_overlap$class_code,overlapping_class_codes)) %>%
  summarise(Delas_best_gene=ref_gene_id[1],
            Delas_best_cc=class_code[1],
            Delas_genes=paste(unique(ref_gene_id),collapse = ","),
            Delas_ccs=paste(unique(class_code),collapse = ","))

# merge with gene level info
source("scripts/load_gene_level_data_biotypes_of_interest.R")
tmap_delas_overlap_summary <- left_join(tmap_delas_overlap_summary,gene_level_info %>% select(gene_name, biotype))

tmap_delas_overlap_summary <- left_join(tmap_delas_overlap_summary,
                                        delas_genes_mm39 %>%select(gene_id,lnc_name,orientation,mm10coord),by=c(Delas_best_gene="gene_id"))



# add annotation of Delas genes relative to enrichment
Delas_enrichment_sheet_id=list("AML_vs_rest"=1,
    "HSC_enriched"=2,
    "lymphoid_enriched"=3,
    "progenitor_vs_differentiated"=4)

Delas_enrichment <- lapply(names(Delas_enrichment_sheet_id),function(da){
  sheetnumber=Delas_enrichment_sheet_id[[da]]
  readxl::read_xls("data/public_data/lncRNA_studies/Delas_et_al/elife-25607-supp2-v2.xls",
                   sheet = sheetnumber,col_names = F)

})


names(Delas_enrichment) <- names(Delas_enrichment_sheet_id)

Delas_enrichment <- sapply(Delas_enrichment,function(da)
  unlist(strsplit(as.data.frame(da)[,1],split = ",")))

for (na in names(Delas_enrichment)) {
  tmap_delas_overlap_summary[,na] <- tmap_delas_overlap_summary$Delas_best_gene %in%
    Delas_enrichment[[na]]

}
summ_enrichment <- apply(tmap_delas_overlap_summary[,names(Delas_enrichment)],
                         1,function(x)paste0(names(Delas_enrichment)[x],
                                             collapse = ","))


tmap_delas_overlap_summary$Delas_enrichment <- summ_enrichment


for (na in names(Delas_enrichment)) {
  delas_genes[,na] <- delas_genes$lnc %in%
    Delas_enrichment[[na]]

}
delas_genes=as.data.frame(delas_genes)
summ_enrichment <- apply(delas_genes[,names(Delas_enrichment)],
                         1,function(x)paste0(names(Delas_enrichment)[x],
                                             collapse = ","))


delas_genes$Delas_enrichment <- summ_enrichment


# write summary table Delas ----
dir.create("outputs/overlap_lncRNA_catalogs")

write.table(tmap_delas_overlap_summary,"outputs/overlap_lncRNA_catalogs/overlap_StemLinc_genes_Delas_genes.txt")




# overlap Luo genes
Luo_503_novel <- read.table("data/public_data/lncRNA_studies/Luo_et_al/503_novel_lncRNAs.liftOver.mm39.bed")
Luo_503_novel_id <- Luo_503_novel %>% mutate(id = paste0(V1,":",V2,"-",V3,":",V6))
Luo_159_LncHSCs <- read.table("data/public_data/lncRNA_studies/Luo_et_al/Goodell_159_LncHSCs.liftOver.mm39.bed")
Luo_159_LncHSCs_id <- Luo_159_LncHSCs %>% mutate(id = paste0(V1,":",V2,"-",V3,":",V6))
table(Luo_159_LncHSCs_id$id%in%Luo_503_novel_id$id)

Luo_503_novel_id <- Luo_503_novel_id %>% mutate(Luo_enrichment=ifelse(id%in%Luo_159_LncHSCs_id$id,
                                                       "HSC_enriched",""))

Luo_503_novel$V4[Luo_503_novel_id$id%in%Luo_159_LncHSCs_id$id] <- Luo_159_LncHSCs_id$V4[match(Luo_503_novel_id$id[Luo_503_novel_id$id%in%Luo_159_LncHSCs_id$id],
                                                                                         Luo_159_LncHSCs_id$id)]

Luo_503_novel_id$V4=Luo_503_novel$V4
# compute overlap with StemLinc genes
gl_bed=sort_bed(extract_bed(gene_level_info))

Luo_intersect <- bt.intersect(gl_bed,sort_bed(Luo_503_novel),wo = T,s = T)
colnames(Luo_intersect)[4] <- "gene_name"

Luo_intersect$Luo_len <- Luo_intersect$V9 - Luo_intersect$V8
Luo_intersect <- Luo_intersect %>% mutate(frac_ol1=round(V13/V5,2),
                                          frac_ol2=round(V13/Luo_len,2))
Luo_intersect <- left_join(Luo_intersect,Luo_503_novel_id %>% select(V4,id,Luo_enrichment),by=c(V10="V4"))

Luo_intersect$sum_ol=Luo_intersect$frac_ol1+Luo_intersect$frac_ol2

Luo_intersect_summ <- Luo_intersect %>% group_by(gene_name) %>% arrange(
                                                                         match(Luo_enrichment,
                                                                               c("HSC_enriched","")),
                                                                         desc(sum_ol)) %>%
  summarise(best_Luo=id[1],
            best_Luo_enrichment=Luo_enrichment[1],
            best_Luo_name=V10[1],
            Luo_ids=paste0(id,collapse = ","),
            Luo_enrichment=paste0(Luo_enrichment,collapse = ","),
            Luo_names=paste0(V10,collapse = ","))

##
write.table(Luo_intersect_summ,"outputs/overlap_lncRNA_catalogs/overlap_StemLinc_genes_Luo_genes.txt",quote = F,row.names = F,sep = "\t")
gene_level_info_Delas_Luo <- full_join(tmap_delas_overlap_summary,
                                        Luo_intersect_summ)
table(gene_level_info_Delas_Luo$biotype,!is.na(gene_level_info_Delas_Luo$Luo_ids))

# write combined data ----
write.table(gene_level_info_Delas_Luo,"outputs/overlap_lncRNA_catalogs/overlap_StemLinc_genes_Delas_and_Luo_genes.txt",
            row.names = F,quote = F,sep = "\t")


delas_genes_not_in_StemLinc <- delas_genes %>% filter(!lnc%in%tmap_delas_overlap_summary$Delas_best_gene)
table(delas_genes_not_in_StemLinc$Delas_enrichment!="")
table(tmap_delas_overlap_summary$Delas_enrichment!="")
