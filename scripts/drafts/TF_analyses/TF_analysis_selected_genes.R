gene_level_info <- read.csv("outputs/gene_level_info_all_evidence_oct24.csv")
mean_corr_co=0.65

gene_level_info_sel <- gene_level_info %>% filter(biotype!="protein_coding"&mean_corr_signature>mean_corr_co)

gene_level_info_sel$gene_name[!is.na(gene_level_info_sel$lnc_name)] <- "Spehd"
gene_level_info_sel$gene_name[gene_level_info_sel$best_Luo_name=="LncHSC-2"] <- "LncHSC-2"

TF_long=read.delim("outputs/TF_enrichment/TF_data_long.txt")

TF_long <- TF_long %>%filter(gene_name%in%gene_level_info_sel$gene_name)

# Sample data structure
# genes_df <- data.frame(
#   gene = c("Gene1", "Gene2", "Gene3", ...),
#   TFs = c("TF1,TF2", "TF2,TF3", "TF1,TF4", ...)
# )

# Example: Extract unique TFs from the data frame
library(dplyr)
library(tidyr)

# Assuming your data frame is named genes_df with columns 'gene' and 'TFs'
unique_TFs <- unique(TF_long$TF)
# View the unique TFs
print(unique_TFs)

N_TFs <- TF_long %>% group_by(gene_name) %>% summarise(N_TFs=length(unique(TF)))

# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")  # Mouse annotation package

library(org.Mm.eg.db)
library(clusterProfiler)

# Assuming your TFs are gene symbols
# Convert gene symbols to Entrez IDs
TF_symbols <- unique_TFs
write.table(TF_symbols,"outputs/TF_enrichment/TF_symbols_selected_candidates.txt",quote = F,row.names = F,col.names = F)
conversion_mouse_human=read.delim('/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/references/conversion_tables_human_mouse/Human_Ensembl_genes_104_GRCh38.p13_NCBI_and_mouse_orthologs.tsv')
table(TF_symbols%in%conversion_mouse_human$Gene.name.y)
TF_entrez <- bitr(TF_symbols,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)

# Remove any NA mappings
TF_entrez <- TF_entrez[!is.na(TF_entrez$ENTREZID), ]

# View the mapped Entrez IDs
head(TF_entrez)
