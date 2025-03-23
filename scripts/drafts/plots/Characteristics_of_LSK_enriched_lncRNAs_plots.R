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
library(plyr)

lncRNAs_enriched_in_LSK <- read.delim("lncRNAs_enriched_in_LSK_log2FC_2.tsv")

# exonic type per class with labels ----
df_csum <- lncRNAs_enriched_in_LSK %>% group_by(class,exonic_type) %>% summarise(count=n())
df_sorted <- arrange(df_csum, class, desc(exonic_type))
head(df_sorted)
df_cumsum <- ddply(df_sorted, "class",
                   transform, label_ypos=cumsum(count))
head(df_cumsum)
ggplot(df_cumsum,aes(x=class,y=count,fill=exonic_type)) + geom_bar(stat="identity") +
  scale_fill_manual(values = nejm_pal) +  geom_text(aes(y=label_ypos, label=count), vjust=1,
                                                    color="white", size=10)+  theme_minimal()+
  theme(text = element_text(size=18))

# genomic position per class ----
gene_classif <- read.delim("gene_classif/all_merged_genes.sel_biotypes.gene_classif.tsv")
lncRNAs_enriched_in_LSK <- left_join(lncRNAs_enriched_in_LSK,gene_classif,by="gene_name")

ggplot(lncRNAs_enriched_in_LSK,aes(class,fill = simpler_classif_to_PCG)) + geom_bar(position = "dodge") +
    theme_minimal() +
  geom_text(
  stat = "count",
  aes(label = after_stat(count)),
  position = position_dodge(width = 0.9),
  vjust = -0.3,
  size = 5
) +   scale_fill_manual(values = nejm_pal) +
  theme(text = element_text(size = 18))

# evidence categories ----
evidence_data <- read.delim("overlap_marks/expressed_merged_genes_annotation.logic_marks.tsv")
lncRNAs_enriched_in_LSK <- left_join(lncRNAs_enriched_in_LSK,evidence_data%>%select(gene_name,evidence_level),by="gene_name")

lncRNAs_enriched_in_LSK$evidence_level <- factor(lncRNAs_enriched_in_LSK$evidence_level,levels = c("none","both","DNA","RNA"))
ggplot(lncRNAs_enriched_in_LSK,aes(x=class,fill=evidence_level)) +
  geom_bar(col="black", position = "fill")+
  theme_minimal() + ylab("Fraction")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("white","#E18727FF","#0072B5FF", "#20854EFF"))+
  theme(text = element_text(size = 18))

ggplot(lncRNAs_enriched_in_LSK,aes(x=class,fill=evidence_level)) +
  geom_bar(col="black")+
  theme_minimal() + ylab("Count")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("white","#E18727FF","#0072B5FF", "#20854EFF"))


# Distance to closest PCG ----
ggplot(lncRNAs_enriched_in_LSK,aes(x=dist_closest_PCG+1,col=class)) +geom_density() + scale_x_log10() +
  theme_minimal()
ggplot(lncRNAs_enriched_in_LSK,aes(x=class,y=dist_closest_PCG+1)) +geom_boxplot() + scale_y_log10() +
  theme_minimal() +theme(text = element_text(size = 18))

lncRNAs_enriched_in_LSK %>% dplyr::group_by(class) %>% dplyr::summarise(median_dist_PCG=median(dist_closest_PCG),
                                                                        mean_dist_PCG=mean(dist_closest_PCG))
# GO enrichment of closest PCGs ----
library(clusterProfiler)
library(org.Mm.eg.db)
#library(enrichplot)

go_enrichment <- enrichGO(gene = lncRNAs_enriched_in_LSK$closest_PCG,
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)

# View the top GO terms
head(go_enrichment)

# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) +
  ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")


# Check TFs bound to these lncRNAs ----
TF_binding <- read.delim("overlap_marks/all_merged_transcripts_annotation.sel_biotypes.TFs_overlap_500_250.gene_level.tsv")

lncRNAs_enriched_in_LSK <- left_join(lncRNAs_enriched_in_LSK,TF_binding)

all_TFs=unlist(strsplit(TF_binding$TF_rel.lev_A,split = ","))
# Assuming your TFs are gene symbols
# Convert gene symbols to Entrez IDs
TF_symbols <- unlist(strsplit(lncRNAs_enriched_in_LSK$TF_rel.lev_A,split = ","))
#write.table(TF_symbols,"outputs/TF_enrichment/TF_symbols_selected_candidates.txt",quote = F,row.names = F,col.names = F)
conversion_mouse_human=read.delim('/run/user/1608803857/gvfs/smb-share:server=10.110.20.7,share=bdcuartero/references/conversion_tables_human_mouse/Human_Ensembl_genes_104_GRCh38.p13_NCBI_and_mouse_orthologs.tsv')
table(unique(TF_symbols)%in%conversion_mouse_human$Gene.name.y)
table(unique(all_TFs)%in%conversion_mouse_human$Gene.name.y)

TF_entrez <- bitr(conversion_mouse_human$Mouse.gene.name[match(TF_symbols,
                                                               conversion_mouse_human$Gene.name.y)],
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)

# Remove any NA mappings
TF_entrez <- TF_entrez[!is.na(TF_entrez$ENTREZID), ]


# Fisher + go enrichment ----
TFs_other_genes <- unlist(strsplit(TF_binding$TF_rel.lev_A[!TF_binding$gene_name%in%lncRNAs_enriched_in_LSK$gene_name],
                            split = ","))

TF_names= sort(unique(c(TF_symbols,TFs_other_genes)))

lncRNA_binding_counts <- table(factor(TF_symbols,levels = TF_names))
background_binding_counts <- table(factor(TFs_other_genes,levels = TF_names))
tf_table <- data.frame(
  TF = TF_names,
  lncRNA_count = as.vector(lncRNA_binding_counts),
  other_count = as.vector(background_binding_counts)
)

# Total lncRNAs and total background genes
total_lncRNAs <- sum(tf_table$lncRNA_count)
total_other_genes <- sum(tf_table$other_count)

library(purrr)
library(dplyr)

tf_enrichment <- tf_table %>%
  rowwise() %>%
  mutate(pval = fisher.test(matrix(c(
    lncRNA_count,
    other_count,
    total_lncRNAs - lncRNA_count,
    total_other_genes - other_count
  ), nrow = 2))$p.value) %>%
  ungroup() %>%
  mutate(padj = p.adjust(pval, method = "BH"))

library(clusterProfiler)
library(org.Mm.eg.db)

enriched_TFs <- tf_enrichment %>%
  filter(padj < 0.05) %>%
  pull(TF)

table(enriched_TFs%in%conversion_mouse_human$Gene.name.x)
enriched_TFs[!enriched_TFs%in%conversion_mouse_human$Gene.name.x]
human_to_mouse_TFs <- c(
  "AP2C"  = "Tfap2c",
  "COE1"  = "Ebf1",
  "GCR"   = "Nr3c1",
  "HXA9"  = "Hoxa9",
  "HXC9"  = "Hoxc9",
  "NF2L2" = "Nfe2l2",
  "NKX21" = "Nkx2-1",
  "NKX25" = "Nkx2-5",
  "P53"   = "Trp53",
  "PO5F1" = "Pou5f1",
  "TF65"  = "Rela",
  "TFE2"  = "Tcf3",
  "TYY1"  = "Zfp42"
)

table(human_to_mouse_TFs%in%TF_binding$gene_name)

mouse_enriched_TFs1=conversion_mouse_human$Mouse.gene.name[match(enriched_TFs,conversion_mouse_human$Gene.name.x)]
mouse_enriched_TFs2=human_to_mouse_TFs

mouse_enriched_TFs=c(mouse_enriched_TFs1[!is.na(mouse_enriched_TFs1)],mouse_enriched_TFs2)

ego_TFs <- enrichGO(
  gene         = mouse_enriched_TFs,
  OrgDb        = org.Mm.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP", # Or "MF", "CC"
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Plot results with a bar plot
barplot(ego_TFs, showCategory = 10) +
  ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(ego_TFs, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")

tf_enrichment %>% arrange(padj)
