options(scipen = 999)
## Load packages---------------------------
require(tidyverse)
require(data.table)
require(rtracklayer)
library(DESeq2)
library(pheatmap)
args <- commandArgs(trailingOnly = TRUE)

# Check if an input dir is provided
if (length(args) == 0) {
  stop("Please provide the working directory as a command-line argument.")
}

# Setup ----
dir_path <- args[1] # path to wdir
setwd(dir_path) # dir_path="outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation_1"
source("functions/nejm_palette.R")
dds_path <- "DGEA/dds_featureCounts.rds"
DGEA_results_path <- "DGEA/all_merged_genes_annotation.DGEA_Diff.tsv"
gene_annot_path <- "all_merged_genes_annotation.sel_biotypes.pass_filter_col.tsv"


# Load data ----
dds <- readRDS(dds_path)


# Normalize counts
normalized_counts <- counts(dds, normalized = TRUE)

# Fingerprint genes from Chambers et al 2007
HSCs=readxl::read_xls("/home/llorenzi/Descargas/mmc3.xls",sheet = 2)
table(HSCs$`Gene Symbol`%in%rownames(normalized_counts))
Monocytes=readxl::read_xls("/home/llorenzi/Descargas/mmc3.xls",sheet = 6)
table(Monocytes$`Gene Symbol`%in%rownames(normalized_counts))
Tcells=readxl::read_xls("/home/llorenzi/Descargas/mmc3.xls",sheet = 5)
table(Tcells$`Gene Symbol`%in%rownames(normalized_counts))
table(Tcells$`Gene Symbol`%in%HSCs$`Gene Symbol`)
HSCs_filtered=HSCs[!HSCs$`Gene Symbol`%in%c(Monocytes$`Gene Symbol`,Tcells$`Gene Symbol`),]
Monocytes_filtered=Monocytes[!Monocytes$`Gene Symbol`%in%c(HSCs$`Gene Symbol`,Tcells$`Gene Symbol`),]
Tcells_filtered=Tcells[!Tcells$`Gene Symbol`%in%c(HSCs$`Gene Symbol`,Monocytes$`Gene Symbol`),]

table(HSCs_filtered$`Gene Symbol`%in%rownames(normalized_counts))

genes_normalized_counts <- rownames(normalized_counts)

fingerprints_counts <- normalized_counts[rownames(normalized_counts)%in%c(HSCs_filtered$`Gene Symbol`,Monocytes_filtered$`Gene Symbol`,
                                                                          Tcells_filtered$`Gene Symbol`),]

annotation_row <- data.frame(Fingerprint=ifelse(rownames(fingerprints_counts)%in%HSCs_filtered$`Gene Symbol`,
                                                "HSC",
                                                ifelse(rownames(fingerprints_counts)%in%Monocytes_filtered$`Gene Symbol`,
                                                       "Monocyte",
                                                       "Tcell")))
rownames(annotation_row)=rownames(fingerprints_counts)
table(annotation_row$Fingerprint)

library(pheatmap)


# Create annotation for gene clusters


# Create an annotation color list
my_cluster_colors <- c("HSC" = "#E18727FF",
                       "Monocyte" =  "#7876B1FF",
                       "Tcell" = "#20854EFF")

# Create an annotation color list
annotation_colors <- list(Fingerprint = my_cluster_colors)
scaled_counts <- t(scale(t(fingerprints_counts)))  # Center genes across samples
scaled_counts <- scaled_counts[order(match(annotation_row$Fingerprint,c("HSC","Monocyte","Tcell"))),]
pheatmap(scaled_counts,
         cluster_rows = F, cluster_cols = FALSE,
         show_rownames = FALSE,  # Hide gene names for clarity
         annotation_row = annotation_row,  # Color genes by cluster,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(100),
         main = "")
