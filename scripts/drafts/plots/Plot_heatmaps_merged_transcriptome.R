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
gene_annot_path <- "expressed_merged_genes_annotation.sel_biotypes.assembled_col.CodingProb.olClass.tsv"


# Load data ----
dds <- readRDS(dds_path)


# Normalize counts
normalized_counts <- counts(dds, normalized = TRUE)

# Select genes to plot
gene_annot <- read.delim(gene_annot_path)
DGEA_results <- read.delim(DGEA_results_path)
gene_annot <- left_join(gene_annot,DGEA_results)
# Genes to plot: pass filter and are differential in at least 1 comparison
gene_annot <- gene_annot %>%filter(!(biotype=="potNovel"&Coding_prob>0.44))
genes2plot <- gene_annot%>%filter(LSK_vs_T.cell!="NS"|LSK_vs_macrophage!="NS"|T.cell_vs_macrophage!="NS")
table(genes2plot$biotype)

# Clustered heatmap by biotype
biot="potNovel"

sig_genes <- genes2plot$gene_name[genes2plot$biotype==biot]

sig_gene_counts <- normalized_counts[sig_genes, ]

# Scale expression values (Z-score transformation)
scaled_counts <- t(scale(t(sig_gene_counts)))  # Center genes across samples

# Hierarchical clustering of genes
gene_dist <- dist(scaled_counts, method = "euclidean")  # Distance metric
gene_clust <- hclust(gene_dist, method = "ward.D2")  # Clustering method

# Cut tree into K clusters (e.g., 4 clusters)
k_clusters <- 4
gene_clusters <- cutree(gene_clust, k = k_clusters)

# Convert to dataframe
gene_cluster_df <- data.frame(Gene = rownames(scaled_counts), Cluster = gene_clusters)
head(gene_cluster_df)
gc_annot <- left_join(gene_cluster_df%>%select(gene_name=Gene,Cluster),gene_annot)
library(pheatmap)
scaled_counts <- scaled_counts[order(gene_cluster_df$Cluster),]
# Create annotation for gene clusters
annotation_row <- data.frame(Cluster = as.factor(gene_clusters))
#rownames(annotation_row) <- rownames(scaled_counts)
# Define colors manually
my_cluster_colors <- c("1" = "#E18727FF",
                       "2" = "#20854EFF",
                       "3" =  "#7876B1FF",
                       "4" =  "#FFDC91FF")

# Create an annotation color list
annotation_colors <- list(Cluster = my_cluster_colors)

colnames(scaled_counts) <- c(paste0("LSK_",1:3),
                             paste0("macrophage_",1:3),
                             paste0("T-cell_",1:3))

pheatmap(scaled_counts,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE,  # Hide gene names for clarity
         annotation_row = annotation_row,  # Color genes by cluster,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(100),
         main = "")

svg(paste0("plots/Heatmap_differential_genes_",biot,".svg"),width = 8,height = 10)
pheatmap(scaled_counts,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE,  # Hide gene names for clarity
         annotation_row = annotation_row,  # Color genes by cluster,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(100),
         main = "",fontsize = 30)
dev.off()

biots=c("potNovel","lncRNA","pseudogene","TEC","protein_coding")
for (biot in biots) {
  sig_genes <- genes2plot$gene_name[genes2plot$biotype==biot]

  sig_gene_counts <- normalized_counts[sig_genes, ]

  # Scale expression values (Z-score transformation)
  scaled_counts <- t(scale(t(sig_gene_counts)))  # Center genes across samples

  # Hierarchical clustering of genes
  gene_dist <- dist(scaled_counts, method = "euclidean")  # Distance metric
  gene_clust <- hclust(gene_dist, method = "ward.D2")  # Clustering method

  # Cut tree into K clusters (e.g., 4 clusters)
  k_clusters <- 4
  gene_clusters <- cutree(gene_clust, k = k_clusters)

  # Convert to dataframe
  gene_cluster_df <- data.frame(Gene = rownames(scaled_counts), Cluster = gene_clusters)
  head(gene_cluster_df)
  print(biot)
  print(table(gene_cluster_df$Cluster))
  scaled_counts <- scaled_counts[order(gene_cluster_df$Cluster),]

  # Create annotation for gene clusters
  annotation_row <- data.frame(Cluster = as.factor(gene_clusters))
  #annotation_row <- annotation_row[order(annotation_row$Cluster),]
  #rownames(annotation_row) <- rownames(scaled_counts)
  # Define colors manually
  my_cluster_colors <- c("1" = "#E18727FF",
                         "2" = "#20854EFF",
                         "3" =  "#7876B1FF",
                         "4" =  "#FFDC91FF")

  # Create an annotation color list
  annotation_colors <- list(Cluster = my_cluster_colors)

  colnames(scaled_counts) <- c(paste0("LSK_",1:3),
                               paste0("macrophage_",1:3),
                               paste0("T-cell_",1:3))

  # pheatmap(scaled_counts,
  #          cluster_rows = FALSE, cluster_cols = FALSE,
  #          show_rownames = FALSE,  # Hide gene names for clarity
  #          annotation_row = annotation_row,  # Color genes by cluster,
  #          annotation_colors = annotation_colors,
  #          color = colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(100),
  #          main = "")

  svg(paste0("plots/Heatmap_differential_genes_",biot,".svg"),width = 8,height = 10)
  pheatmap(scaled_counts,
           cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = FALSE,  # Hide gene names for clarity
           annotation_row = annotation_row,  # Color genes by cluster,
           annotation_colors = annotation_colors,
           color = colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(100),
           main = "",fontsize = 30)
  dev.off()

}



# add log2FC ----
res_LSKvsT.cell <- as.data.frame(results(dds, contrast = c("sample_type", "LSK","T.cell")))
res_LSKvsT.cell <- res_LSKvsT.cell%>%filter(!is.na(padj))
res_LSKvsT.cell$gene_name=rownames(res_LSKvsT.cell)
gene_annot <- left_join(gene_annot,res_LSKvsT.cell%>%select(gene_name,l2FC_LSKvsT.cell=log2FoldChange,padj_LSKvsT.cell=padj))


res_LSKvsmacrophage <- as.data.frame(results(dds, contrast = c("sample_type", "LSK","macrophage")))
res_LSKvsmacrophage <- res_LSKvsmacrophage%>%filter(!is.na(padj))
res_LSKvsmacrophage$gene_name=rownames(res_LSKvsmacrophage)
gene_annot <- left_join(gene_annot,res_LSKvsmacrophage%>%select(gene_name,l2FC_LSKvsmacrophage=log2FoldChange,padj_LSKvsmacrophage=padj))


res_T.cellvsmacrophage <- as.data.frame(results(dds, contrast = c("sample_type", "T.cell","macrophage")))
res_T.cellvsmacrophage <- res_T.cellvsmacrophage%>%filter(!is.na(padj))
res_T.cellvsmacrophage$gene_name=rownames(res_T.cellvsmacrophage)
gene_annot <- left_join(gene_annot,res_T.cellvsmacrophage%>%select(gene_name,l2FC_T.cellvsmacrophage=log2FoldChange,padj_T.cellvsmacrophage=padj))

lFCco=2
genes2plot_high_FC <- gene_annot%>%filter((LSK_vs_T.cell!="NS"|LSK_vs_macrophage!="NS"|T.cell_vs_macrophage!="NS"),
                                          abs(l2FC_LSKvsT.cell)>lFCco|abs(l2FC_LSKvsmacrophage)>lFCco|abs(l2FC_T.cellvsmacrophage)>lFCco )
table(genes2plot_high_FC$biotype,genes2plot_high_FC$exonic_type)
table(genes2plot_high_FC$class)
table(genes2plot$biotype)


# High FC genes ----
# Clustered heatmap by biotype
biot="potNovel"

sig_genes <- genes2plot_high_FC$gene_name[genes2plot_high_FC$biotype==biot]

sig_gene_counts <- normalized_counts[sig_genes, ]

# Scale expression values (Z-score transformation)
scaled_counts <- t(scale(t(sig_gene_counts)))  # Center genes across samples

# Hierarchical clustering of genes
gene_dist <- dist(scaled_counts, method = "euclidean")  # Distance metric
gene_clust <- hclust(gene_dist, method = "ward.D2")  # Clustering method

# Cut tree into K clusters (e.g., 4 clusters)
k_clusters <- 4
gene_clusters <- cutree(gene_clust, k = k_clusters)

# Convert to dataframe
gene_cluster_df <- data.frame(Gene = rownames(scaled_counts), Cluster = gene_clusters)
head(gene_cluster_df)
gc_annot <- left_join(gene_cluster_df%>%select(gene_name=Gene,Cluster),gene_annot)
library(pheatmap)
scaled_counts <- scaled_counts[order(gene_cluster_df$Cluster),]
# Create annotation for gene clusters
annotation_row <- data.frame(Cluster = as.factor(gene_clusters))
#rownames(annotation_row) <- rownames(scaled_counts)
# Define colors manually
my_cluster_colors <- c("1" = "#E18727FF",
                       "2" = "#20854EFF",
                       "3" =  "#7876B1FF",
                       "4" =  "#FFDC91FF")

# Create an annotation color list
annotation_colors <- list(Cluster = my_cluster_colors)

colnames(scaled_counts) <- c(paste0("LSK_",1:3),
                             paste0("macrophage_",1:3),
                             paste0("T-cell_",1:3))

pheatmap(scaled_counts,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE,  # Hide gene names for clarity
         annotation_row = annotation_row,  # Color genes by cluster,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(100),
         main = "")

svg(paste0("plots/Heatmap_differential_genes_logFC2",biot,".svg"),width = 8,height = 10)
pheatmap(scaled_counts,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE,  # Hide gene names for clarity
         annotation_row = annotation_row,  # Color genes by cluster,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(100),
         main = "",fontsize = 30)
dev.off()

biots=c("potNovel","lncRNA","pseudogene","TEC","protein_coding")
gcs=list()
for (biot in biots) {
  sig_genes <- genes2plot_high_FC$gene_name[genes2plot_high_FC$biotype==biot]

  sig_gene_counts <- normalized_counts[sig_genes, ]

  # Scale expression values (Z-score transformation)
  scaled_counts <- t(scale(t(sig_gene_counts)))  # Center genes across samples

  # Hierarchical clustering of genes
  gene_dist <- dist(scaled_counts, method = "euclidean")  # Distance metric
  gene_clust <- hclust(gene_dist, method = "ward.D2")  # Clustering method

  # Cut tree into K clusters (e.g., 4 clusters)
  k_clusters <- 4
  gene_clusters <- cutree(gene_clust, k = k_clusters)
  # Convert to dataframe
  gene_cluster_df <- data.frame(Gene = rownames(scaled_counts), Cluster = gene_clusters)
  gcs[[biot]] <- left_join(gene_cluster_df%>%select(gene_name=Gene,Cluster),gene_annot)

  head(gene_cluster_df)
  print(biot)
  print(table(gene_cluster_df$Cluster))
  scaled_counts <- scaled_counts[order(gene_cluster_df$Cluster),]

  # Create annotation for gene clusters
  annotation_row <- data.frame(Cluster = as.factor(gene_clusters))
  #annotation_row <- annotation_row[order(annotation_row$Cluster),]
  #rownames(annotation_row) <- rownames(scaled_counts)
  # Define colors manually
  my_cluster_colors <- c("1" = "#E18727FF",
                         "2" = "#20854EFF",
                         "3" =  "#7876B1FF",
                         "4" =  "#FFDC91FF")

  # Create an annotation color list
  annotation_colors <- list(Cluster = my_cluster_colors)

  colnames(scaled_counts) <- c(paste0("LSK_",1:3),
                               paste0("macrophage_",1:3),
                               paste0("T-cell_",1:3))

  # pheatmap(scaled_counts,
  #          cluster_rows = FALSE, cluster_cols = FALSE,
  #          show_rownames = FALSE,  # Hide gene names for clarity
  #          annotation_row = annotation_row,  # Color genes by cluster,
  #          annotation_colors = annotation_colors,
  #          color = colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(100),
  #          main = "")

  svg(paste0("plots/Heatmap_differential_genes_logFC2",biot,".svg"),width = 8,height = 10)
  pheatmap(scaled_counts,
           cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = FALSE,  # Hide gene names for clarity
           annotation_row = annotation_row,  # Color genes by cluster,
           annotation_colors = annotation_colors,
           color = colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(100),
           main = "",fontsize = 30)
  dev.off()

}
table(gcs$potNovel$Cluster,gcs$potNovel$exonic_type)
table(gcs$potNovel$Cluster,gcs$potNovel$class)

table(gcs$lncRNA$Cluster,gcs$lncRNA$exonic_type)

lncRNAs_enriched_in_LSK <- rbind(gcs$potNovel%>%filter(Cluster==1),
                                 gcs$lncRNA%>%filter(Cluster==2))
lncRNAs_enriched_in_LSK <- rbind(lncRNAs_enriched_in_LSK,gcs$TEC%>%filter(Cluster==3))
lncRNAs_enriched_in_LSK <- lncRNAs_enriched_in_LSK%>%filter(LSK_vs_T.cell!="DOWN"&LSK_vs_macrophage!="DOWN")
write.table(lncRNAs_enriched_in_LSK,"lncRNAs_enriched_in_LSK_log2FC_2.tsv",
            quote = F,row.names = F,sep = "\t")


