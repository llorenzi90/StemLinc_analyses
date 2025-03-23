# Load necessary libraries
library(biomaRt)
library(tidyverse)
library(GO.db)


# Get genes related with hematopoiesis GO terms ----
# Get all GO terms
all_go_terms <- as.list(GOTERM)

# Extract GO ID, term name, and definition
go_info <- data.frame(
  GO_ID = names(all_go_terms),
  Term = sapply(all_go_terms, Term),
  Definition = sapply(all_go_terms, Definition)
)

print( go_info$Term[grep("hemat|hemopoiesis",go_info$Term)])
# [1] "hematopoietic progenitor cell differentiation"
# [2] "hemopoiesis"
# [3] "embryonic hemopoiesis"
# [4] "post-embryonic hemopoiesis"
# [5] "larval lymph gland hemopoiesis"
# [6] "fibroblast growth factor receptor signaling pathway involved in hemopoiesis"
# [7] "hematopoietic stem cell migration"
# [8] "hematopoietic or lymphoid organ development"
# [9] "primitive hemopoiesis"
# [10] "definitive hemopoiesis"
# [11] "hematopoietic stem cell differentiation"
# [12] "hematopoietic stem cell homeostasis"
# [13] "hematopoietic stem cell proliferation"
# [14] "hematopoietic stem cell migration to bone marrow"
# [15] "endothelial to hematopoietic transition"
# [16] "regulation of hematopoietic progenitor cell differentiation"
# [17] "negative regulation of hematopoietic progenitor cell differentiation"
# [18] "positive regulation of hematopoietic progenitor cell differentiation"
# [19] "regulation of hematopoietic stem cell proliferation"
# [20] "negative regulation of hematopoietic stem cell proliferation"
# [21] "positive regulation of hematopoietic stem cell proliferation"
# [22] "regulation of hematopoietic stem cell differentiation"
# [23] "negative regulation of hematopoietic stem cell differentiation"
# [24] "positive regulation of hematopoietic stem cell differentiation"
# [25] "regulation of hemopoiesis"
# [26] "negative regulation of hemopoiesis"
# [27] "positive regulation of hemopoiesis"
# [28] "regulation of hematopoietic stem cell migration"
# [29] "negative regulation of hematopoietic stem cell migration"
# [30] "positive regulation of hematopoietic stem cell migration"

hemato_go_terms <- go_info$GO_ID[grep("hemat|hemopoiesis",go_info$Term)]

# Connect to the Ensembl BioMart database
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Retrieve gene information
hematopoiesis_genes <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "go",
  values = hemato_go_terms,
  mart = mart
)


# View results
head(hematopoiesis_genes)

# Load gene signature used for correlation analyses ----
SSignature <- read.table("data/gene_lists/HSC_signature3.txt")
table(SSignature$V1%in%hematopoiesis_genes$external_gene_name)


# Load fingerprint genes from Stuart 2007 ----
# 10.1016/j.stem.2007.10.003
path_Stuart2007 <- "data/public_data/fingerprint_genes_Stuart2007/mmc3.xls"
fingerprint_HSC <- readxl::read_xls(path_Stuart2007,sheet = 2)

table(hematopoiesis_genes$external_gene_name%in%fingerprint_HSC$`Gene Symbol`)
hematopoiesis_genes$external_gene_name[hematopoiesis_genes$external_gene_name%in%fingerprint_HSC$`Gene Symbol`]

# NOTE: there is not great overlap between genes defined by GO terms and
# genes in the Stuart 2007 fingerprints

# Define a big list with all genes ----
big_list <- unique(c(hematopoiesis_genes$external_gene_name,
                     fingerprint_HSC$`Gene Symbol`,SSignature$V1))

# Check how many are in our filtered gene annotation
gene_level_info <- read.csv("outputs/gene_level_info_all_evidence_oct24.csv")
table(big_list%in%gene_level_info$gene_name)
# FALSE  TRUE
# 185   437

# Check individually fingerprints and GO genes

table(hematopoiesis_genes$external_gene_name%in%gene_level_info$gene_name)
# FALSE  TRUE
# 68   247
table(fingerprint_HSC$`Gene Symbol`%in%gene_level_info$gene_name)
# FALSE  TRUE
# 143   230

# Try to recover genes by aliases ----
# Define your list of genes
genes <- big_list[!big_list%in%gene_level_info$gene_name]

# Query Ensembl for gene synonyms (aliases)
aliases <- getBM(
  attributes = c("external_gene_name", "external_synonym"),
  filters = "external_gene_name",
  values = genes,
  mart = mart
)

# View the results
print(aliases)

# Check if aliases are present in our gene_names

table(aliases$external_synonym%in%gene_level_info$gene_name)
# FALSE  TRUE
# 206     3

# Only 3 genes can be recovered in this way...

big_list <- c(big_list, aliases$external_synonym[aliases$external_synonym%in%gene_level_info$gene_name])

table(big_list%in%gene_level_info$gene_name)

# Write results ----
write.table(big_list,"data/gene_lists/HSC_genes.long.txt",quote = F,col.names = F,row.names = F)
write.table(big_list[big_list%in%gene_level_info$gene_name],
            "data/gene_lists/HSC_genes.match_filtered_genes_440.txt",
            quote = F,col.names = F,row.names = F)
