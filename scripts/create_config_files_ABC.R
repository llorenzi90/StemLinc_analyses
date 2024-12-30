## Purpose of script: create mouse reference files to run ABC-Enhancer-Gene-Prediction
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-10-28
##
## Email: lucialorenzi90@gmail.com
##
#  Notes ---------------------------
##  https://abc-enhancer-gene-prediction.readthedocs.io/en/latest/usage/getting_started.html
##
##
#  Setup ---------------------------

options(scipen = 999)
# Load packages---------------------------
# load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
library(rtracklayer)
library(bedtoolsr)
library(httr)
library(biomaRt)
library(AnnotationHub)
library(trastools)
library(gtools)
#library(ggplot2)

source("scripts/source_all_functions.R")
# Create output dir ----
dir.create("outputs/ABC_files")

# Load all ABC hg38 genes ----
url_abc_all_genes <- "https://raw.githubusercontent.com/broadinstitute/ABC-Enhancer-Gene-Prediction/refs/heads/main/reference/hg38/CollapsedGeneBounds.hg38.bed"
abc_all_genes <- read.table(url_abc_all_genes)
table(abc_all_genes$V8)

# Load mm39 LSK annotation ----
#mm39_annotation <- read.table("data/annotated_tracking_file.updated_gene_names.txt",header=T,sep="\t")
#mm39_annotation <- read.table("~/Descargas/annotated_tracking_file.updated_gene_names.txt",header=T,sep="\t")

# Filtered genes:
mm39_LSK_genes <- read.table("outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.biotypes_of_interest.tsv", header = T)
table(mm39_LSK_genes$biotype)

# Format genes to sorted bed ----
mm39_LSK_genes_bed <- mm39_LSK_genes %>%
  mutate(V5 = 0, start=start - 1 , gene_id = gene_name,
         seqnames = factor(seqnames, levels=unique(mixedsort(seqnames)))) %>%
  dplyr::select(seqnames,start,end,
                gene_name,V5, strand, gene_id, biotype) %>% arrange(seqnames,start)

mm39_LSK_genes_bed$gene_id <- paste0("EN",mm39_LSK_genes_bed$gene_id)
# ABC code gives error otherwise:
# def read_gene_bed_file(bed_file):
# BED6 format with extra columns for ensemble info
# columns = BED6_COLS + ["Ensembl_ID", "gene_type"]
# skip = 1 if ("track" in open(bed_file, "r").readline()) else 0
# result = pd.read_table(
#   bed_file, names=columns, header=None, skiprows=skip, comment="#"
# )
# ensembl_id_col = result.loc[:, "Ensembl_ID"]
# if ensembl_id_col.isna().all() or not ensembl_id_col.str.contains("EN").all():
#   raise Exception("Gene file doesn't follow the correct format with Ensembl info")
# return result

# Define gene promoters ----
mm39_LSK_genes_bed_Prom <- get_promoters_from_bed(mm39_LSK_genes_bed,
                                                  before = 250,
                                                  after = 250)
# Define ubiquitous genes ----

# ABC ubiquitous human symbols
url_abc_ub_genes <- "https://raw.githubusercontent.com/broadinstitute/ABC-Enhancer-Gene-Prediction/refs/heads/main/reference/UbiquitouslyExpressedGenes.txt"
abc_ub_genes <- read.table(url_abc_ub_genes)

# HRT Atlas
url_hrt_ub_genes <- "https://housekeeping.unicamp.br/Housekeeping_GenesMouse.csv"
response <- GET(url_hrt_ub_genes, config = config(ssl_verifypeer = FALSE))
hrt_ub_genes <- read.csv(text = content(response, "text"), header = TRUE,sep = ";")

# Check overlap between HRT Atlas and ABC ubiquitous genes
url_hrt_human_ub_genes <- "https://housekeeping.unicamp.br/Housekeeping_GenesHuman.csv"
response <- GET(url_hrt_human_ub_genes, config = config(ssl_verifypeer = FALSE))
hrt_human_ub_genes <- read.csv(text = content(response, "text"), header = TRUE,sep = ";")

# first check how many of the HRT HK genes are in ABC all genes
table(unique(hrt_human_ub_genes$Gene.name)%in%abc_all_genes$V4)
# FALSE  TRUE
# 40  2136
# most of them are

# Now, check how many of the genes considered to be ubiquitous by ABC model
# are also marked as ubiquitous in HRT
table(abc_ub_genes$V1%in%hrt_human_ub_genes$Gene.name)
# only 497 out of 350... not great overlap

# I'll use a simple translation from the ABC model for now...
conversion_table <- fread("~/Descargas/Human_Ensembl_genes_104_GRCh38.p13_NCBI_and_mouse_orthologs.tsv",header = T,data.table=F)

table(abc_ub_genes$V1%in%c(conversion_table$`Gene name.x`,conversion_table$`Gene name.y`))
abc_ub_genes_tr_mouse <- unique(conversion_table %>%
                                   filter(`Gene name.x`%in%abc_ub_genes$V1|`Gene name.y`%in%abc_ub_genes$V1) %>%
                                   pull(`Mouse gene name`))
length(abc_ub_genes_tr_mouse)

table(abc_ub_genes_tr_mouse%in%mm39_LSK_genes_bed$gene_name)
table(abc_ub_genes$V1%in%abc_all_genes$V4)

not_found_abc_ub <- abc_ub_genes$V1[!abc_ub_genes$V1%in%c(conversion_table$`Gene name.x`,
                                      conversion_table$`Gene name.y`)]

# recover some genes with tolower()
recover_2lower <- mm39_LSK_genes_bed$gene_name[tolower(mm39_LSK_genes_bed$gene_name)%in%
                                                 tolower(not_found_abc_ub)]

table(mm39_LSK_genes_bed$gene_name%in%c(abc_ub_genes_tr_mouse,recover_2lower))

not_found_abc_ub <- not_found_abc_ub[!tolower(not_found_abc_ub)%in%tolower(mm39_LSK_genes_bed$gene_name)]

mm39_LSK_ub_genes <- mm39_LSK_genes_bed$gene_name[mm39_LSK_genes_bed$gene_name%in%
                                                    c(abc_ub_genes_tr_mouse,recover_2lower)]

length(mm39_LSK_ub_genes)
# I recovered 791 housekeeping genes in LSK transcriptome
table(mm39_LSK_genes_bed$biotype[mm39_LSK_genes_bed$gene_name%in%mm39_LSK_ub_genes])

mouse_ub_genes <-  c(abc_ub_genes_tr_mouse,recover_2lower)
length(mouse_ub_genes)

# Blacklist regions ----
download.file("https://drive.google.com/uc?export=download&id=18VWCR7QdpPCpSIwMoTTnMnFK4VuL1Uak",
              destfile = "~/mm39.excluderanges.bed",
              mode = "wb")
excluderanges_mm39 <- read.table("~/mm39.excluderanges.bed",sep = "\t")

# check overlap with gene regions
gene_BL_overlap <- bt.intersect(a=mm39_LSK_genes_bed,b=excluderanges_mm39,wo = T)

prom_BL_overlap <- bt.intersect(a=mm39_LSK_genes_bed_Prom,b=excluderanges_mm39,wo = T)

length(unique(gene_BL_overlap$V4))
# check overlap of abc example blacklist with their genes
url_abc_BL <- "https://raw.githubusercontent.com/broadinstitute/ABC-Enhancer-Gene-Prediction/refs/heads/main/reference/hg38/GRCh38_unified_blacklist.bed"
abc_BL <- read.table(url_abc_BL)
abc_gene_BL_overlap <- bt.intersect(abc_all_genes,abc_BL, wo = T)
length(unique(abc_gene_BL_overlap$V4))

genes_high_corr <- read.csv("gene_level_info_114.csv")

table(genes_high_corr$gene_name%in%gene_BL_overlap$V4)
genes_high_corr$gene_name[genes_high_corr$gene_name%in%gene_BL_overlap$V4]
table(genes_high_corr$gene_name%in%prom_BL_overlap$V4)
genes_high_corr$gene_name[genes_high_corr$gene_name%in%prom_BL_overlap$V4]

# only one of the potential novel high correlated
# genes has overlap with blacklisted regions: XLOC_013347

# Create chr sizes file ----
chrs_mm39_bed <- read.table("data/references/genomes/GRCm39_primary_chrs.bed",sep = "\t")
chrs_mm39_bed <- chrs_mm39_bed %>% arrange(match(V1,unique(mm39_LSK_genes_bed$seqnames))) %>%
  dplyr::select(V1,V3)

# write data ----
genome="mm39"
abc_path="~/ABC-Enhancer-Gene-Prediction/"
dir.create(paste0(abc_path,"reference/",genome))
# In the order they appear in the config file
## chrom_sizes ----
write.table(chrs_mm39_bed,"outputs/ABC_files/mm39_noalt.chrom.sizes.tsv",
            quote = F,col.names = F,row.names = F,sep = "\t")
# copy to ABC folder
system(paste0("cp outputs/ABC_files/mm39_noalt.chrom.sizes.tsv ",
              paste0(abc_path,"reference/",genome)))

## blacklist ----
write.table(excluderanges_mm39[,1:3],
            "outputs/ABC_files/mm39.excluderanges.short.bed",
            quote = F,col.names = F,row.names = F,sep = "\t")
# copy to ABC folder
system(paste0("cp outputs/ABC_files/mm39.excluderanges.short.bed ",
              paste0(abc_path,"reference/",genome)))

## ubiquitous genes ----
write.table(mouse_ub_genes,"outputs/ABC_files/UbiquitouslyExpressedGenes.Mouse.txt",
            quote = F,col.names = F,row.names = F,sep = "\t")

# copy to ABC folder
system(paste0("cp outputs/ABC_files/UbiquitouslyExpressedGenes.Mouse.txt ",
              paste0(abc_path,"reference/")))

## list of genes corresponding with the genome ----
write.table(mm39_LSK_genes_bed,"outputs/ABC_files/LSK_CollapsedGeneBounds.mm39.bed",
            quote = F,col.names = F,row.names = F,sep = "\t")

# copy to ABC folder
system(paste0("cp outputs/ABC_files/LSK_CollapsedGeneBounds.mm39.bed ",
              paste0(abc_path,"reference/",genome)))

## list of promoters corresponding with the genome ----
write.table(mm39_LSK_genes_bed_Prom,"outputs/ABC_files/LSK_CollapsedGeneBounds.mm39.TSS500bp.bed",
            quote = F,col.names = F,row.names = F,sep = "\t")

# copy to ABC folder
system(paste0("cp outputs/ABC_files/LSK_CollapsedGeneBounds.mm39.TSS500bp.bed ",
              paste0(abc_path,"reference/",genome)))
