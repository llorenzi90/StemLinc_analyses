# analyse combined transcriptome:
# annotate genes in terms of overlap Ref, biotype and classification relative to closest PCG
# make script for this?? for now is the combination of rmd script,
# For the ones in LSK, track the name to this file, for the new names add a prefix to avoid confusions
# check how many genes in each sample
# perform PCA and DGEA
library(trastools)
library(tidyverse)
source("scripts/functions/nejm_palette.R")
get_n_samples_per_gene_from_tracking_gene_name <- function(tracking,gene_id="gene_name",cols=6:8){
  tracking <- as.data.frame(tracking)
  tracking$V2=tracking[,gene_id]
  get_n_samples_per_gene_from_tracking(tracking,cols = cols)
}
ordered_classcodes <- c("=", "j", "k", "c", "m", "n", "e", "o", "p", "s", "x", "i", "y", "u", "r", ".")

LSK_annotated_tracking_path="outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.tsv"
tracking_path="data/hpc_data/filtered_LSK_T-cell_macrophage.tracking"
annotated_tracking_path <- "outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/filtered_LSK_T-cell_macrophage.combined_annotated_tracking.tsv"
ref_annot_path="data/references/merged_refs_annotation/annotated_tracking_file.updated_gene_names.txt"
tpm_path <- "data/hpc_data/filtered_assembly_StemLinc_samples_kallisto.TPM.tsv"
counts_path <- "data/hpc_data/filtered_assembly_StemLinc_samples_kallisto.counts.tsv"

tracking=read.table(tracking_path)
ref_annot=read.table(ref_annot_path,header = T)
LSK_tracking=read.table(LSK_annotated_tracking_path,header = T)

genes_in_samples <- trastools::get_sample_occurrence_per_gene_from_tracking(tracking)
colSums(genes_in_samples[,-4])
table(rowSums(genes_in_samples[,-4]))

colSums(tracking[,5:ncol(tracking)]!="-")

# Read annotated tracking ----
annotated_tracking <- read.table(annotated_tracking_path,header = T)
# change gene_id and transcript_id to avoid mixing with old ones ----
annotated_tracking$V1 <- gsub("TCONS","STTRA",annotated_tracking$V1)

annotated_tracking$V2 <- gsub("XLOC","STLOC",annotated_tracking$V2)

# For the ones in LSK, track the name to this file ----
splitted_names_LSK <- separate(annotated_tracking,col = "V5",into = paste0("f",1:7),sep = "\\|")
splitted_names_LSK$f1 <- gsub("q1:","",splitted_names_LSK$f1)
LSK_conversion_table <- splitted_names_LSK %>% select(transcript_id=V1, overlap_Ref,LSK_transcript=f2,
                                                      ref_gene_name, LSK_gene_name=f1)
# LSK_conversion_table_simpl <- LSK_conversion_table%>%filter(gene_name!=LSK_gene_name)
# View(LSK_conversion_table_simpl[!grepl("STLOC",LSK_conversion_table_simpl$gene_name)&LSK_conversion_table_simpl$LSK_gene_name!="-",])
# View(LSK_conversion_table_simpl[grepl("STLOC",LSK_conversion_table_simpl$gene_name)&LSK_conversion_table_simpl$LSK_gene_name!="-",])
# new_genes_conversion <- LSK_conversion_table_simpl[grepl("STLOC",LSK_conversion_table_simpl$gene_name)&LSK_conversion_table_simpl$LSK_gene_name!="-",]
# table(!duplicated(new_genes_conversion[,2:3]))
# new_genes_conversion_undup <- new_genes_conversion[!duplicated(new_genes_conversion[,2:3]),]

# Assign meaningful gene names
# NOTE: I'll take the approach where I select first the assigned
# reference gene name for the consensus transcript if classcode is overlapping
# if not overlapping, I take the LSK name first,
annotated_tracking <- left_join(annotated_tracking, LSK_conversion_table%>%select(transcript_id,LSK_gene_name),
                                by=c(V1="transcript_id"))
annotated_tracking <- annotated_tracking %>% mutate(gene_name = ifelse(V4%in%trastools::overlapping_class_codes,
                                                                       ref_gene_name,ifelse(LSK_gene_name!="-",LSK_gene_name,V2)))


# annotate genes ----
N_samps_per_gene=get_n_samples_per_gene_from_tracking_gene_name(annotated_tracking)
N_samps_per_gene <- data.frame(gene_name=names(N_samps_per_gene),Nsamps=N_samps_per_gene)

best_cc_per_gene=annotated_tracking %>% group_by(gene_name) %>% arrange(match(V4,ordered_classcodes)) %>%
  summarise(best_cc=V4[1])

gene_level_info <- left_join(N_samps_per_gene,best_cc_per_gene)
gene_level_info$biotype=ref_annot$simplified_gene_biotype[match(gene_level_info$gene_name,
                                                                ref_annot$gene_name)]
gene_level_info <- gene_level_info %>% mutate(gene_class=ifelse(is.na(biotype),best_cc,biotype))

gene_level_genomic_range <- trastools::get_genomic_range_by_gene(annotated_tracking,by = "gene_name")
gene_level_info <- left_join(gene_level_info, gene_level_genomic_range)
gene_level_info <- gene_level_info %>% mutate(biotype = ifelse(is.na(biotype),"potNovel",biotype))
table(gene_level_info$biotype)

# load TPM expression ----
TPMs <- read.table(tpm_path,header = T)
TPMs$transcript_id <- gsub("TCONS","STTRA",TPMs$transcript_id)
TPMs$gene_name <- annotated_tracking$gene_name[match(TPMs$transcript_id,
                                                     annotated_tracking$V1)]

# Aggregate at gene level ----
TPM_gene_level <- aggregate(TPMs[,2:(ncol(TPMs)-1)],by = list(gene_name=TPMs$gene_name),sum)



# load count expression ----
counts <- read.table(counts_path,header = T)
counts$transcript_id <- gsub("TCONS","STTRA",counts$transcript_id)
counts$gene_name <- annotated_tracking$gene_name[match(counts$transcript_id,
                                                     annotated_tracking$V1)]

# Aggregate at gene level ----
counts_gene_level <- aggregate(counts[,2:(ncol(counts)-1)],by = list(gene_name=counts$gene_name),sum)

# Select for biotypes of interest only
biots=c("protein_coding","lncRNA","potNovel","TEC","pseudogene")

matched_biotypes <- gene_level_info$biotype[match(counts_gene_level$gene_name,
                                                  gene_level_info$gene_name)]

counts_gene_level <- counts_gene_level %>%filter(matched_biotypes%in%biots)
# DEGEA ----
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Load the count matrix and metadata
gene_names <- counts_gene_level$gene_name
counts_gene_level <- counts_gene_level[,-1]
coldata <- data.frame(sample_id=colnames(counts_gene_level),
                      sample=rep(c("LSK","macrophage","T-cell"),each=3),
                      replicate=rep(1:3,3))
rownames(coldata) <- coldata$sample_id
# Ensure sample names in metadata match column names in count matrix
all(colnames(counts_gene_level) == rownames(coldata))

# convert counts to integer
counts_gene_level <- apply(counts_gene_level,2,as.integer)
rownames(counts_gene_level) <- gene_names
# generate DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts_gene_level, colData = coldata, design = ~ sample)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
resultsNames(dds)
rowData(dds)$biotype <- gene_level_info$biotype[match(rownames(dds),
                                                      gene_level_info$gene_name)]

# Transform data with variance stabilizing transformation (VST):

vsd <- vst(dds, blind = FALSE)

# Perform PCA and visualize:

pca_data <- plotPCA(vsd, intgroup = "sample", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = sample)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_minimal()


# PCA per biotype
for (biot in biots) {
  pca_tmp <- plotPCA(vsd[rowData(dds)$biotype == biot,],
                     intgroup = "sample", returnData = TRUE)
  percent_var <- round(100 * attr(pca_tmp, "percentVar"))

 g = ggplot(pca_tmp, aes(PC1, PC2, color = sample)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    theme_minimal() + ggtitle(paste("Biotype :",biot))

 print(g)

}


# Sort by significance and plot heatmap for each biotype:
results_list=list()
for(biot in biots){

  # T-cell
  res_T.cellvsLSK <- results(dds[rowData(dds)$biotype == biot,], contrast = c("sample",
                                                                              "T.cell", "LSK"))
  res_name=paste0(biot,"T.cell")
  res_T.cellvsLSK <- res_T.cellvsLSK[order(res_T.cellvsLSK$padj), ]
  summary(res_T.cellvsLSK)
  res_T.cellvsLSK <- as.data.frame(res_T.cellvsLSK)
  res_T.cellvsLSK <- res_T.cellvsLSK %>% mutate(Diff=ifelse(padj<0.05,ifelse(log2FoldChange<0,"DOWN","UP"),"NS"))
  res_T.cellvsLSK <- res_T.cellvsLSK %>% filter(!is.na(padj))
  results_list[[res_name]] <- as.data.frame(res_T.cellvsLSK)
  #write.csv(as.data.frame(res), "differential_genes.csv")
  g <- ggplot(res_T.cellvsLSK, aes(x = log2FoldChange, y = -log10(padj), col=Diff)) +
    geom_point() +
    #geom_hline(y = -log10(sig_threshold), linetype = "dashed") +
    #geom_vline(x = c(-fc_threshold, fc_threshold), linetype = "dashed") +
    xlab("Log2 Fold Change") +
    ylab("-log10 Adjusted p-value") +
    ggtitle(paste("Volcano Plot T-cell vs LSK",biot)) +
    scale_color_manual(values = nejm_pal[c(2,6,1)]) +
    theme_bw()
  print(g)
  # macrophage
  res_macrophagevsLSK <- results(dds[rowData(dds)$biotype == biot,], contrast = c("sample",
                                                                                  "macrophage", "LSK"))
  res_name=paste0(biot,"macrophage")
  res_macrophagevsLSK <- res_macrophagevsLSK[order(res_macrophagevsLSK$padj), ]
  summary(res_macrophagevsLSK)
  res_macrophagevsLSK <- as.data.frame(res_macrophagevsLSK)
  res_macrophagevsLSK <- res_macrophagevsLSK %>% mutate(Diff=ifelse(padj<0.05,ifelse(log2FoldChange<0,"DOWN","UP"),"NS"))
  res_macrophagevsLSK <- res_macrophagevsLSK %>% filter(!is.na(padj))
  results_list[[res_name]] <- as.data.frame(res_macrophagevsLSK)
  #write.csv(as.data.frame(res), "differential_genes.csv")
  g <- ggplot(res_macrophagevsLSK, aes(x = log2FoldChange, y = -log10(padj), col=Diff)) +
    geom_point() +
    #geom_hline(y = -log10(sig_threshold), linetype = "dashed") +
    #geom_vline(x = c(-fc_threshold, fc_threshold), linetype = "dashed") +
    xlab("Log2 Fold Change") +
    ylab("-log10 Adjusted p-value") +
    ggtitle(paste("Volcano Plot macrophage vs LSK",biot)) +
    scale_color_manual(values = nejm_pal[c(2,6,1)]) +
    theme_bw()
  print(g)
  #Filter significant genes (adjusted p-value < 0.05):

  sig_genes1 <- res_T.cellvsLSK[which(res_T.cellvsLSK$padj < 0.05), ]
  sig_genes2 <- res_macrophagevsLSK[which(res_macrophagevsLSK$padj<0.05),]

  table(rownames(sig_genes1)%in%rownames(sig_genes2))
  sig_genes <- unique(rownames(sig_genes1),rownames(sig_genes2))

  #sig_genes1 <- as.data.frame(sig_genes1)
  # Step 5: Generate a Heatmap
  # Extract normalized counts for significant genes:


  normalized_counts <- counts(dds, normalized = TRUE)
  sig_gene_counts <- normalized_counts[sig_genes, ]

  #Scale data and plot heatmap:

  pheatmap(sig_gene_counts,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = FALSE,
           scale = "row",
           #annotation_col = coldata,
           main = paste("Heatmap of",biot, "Differentially Expressed Genes"))

}

# mañana:
# hacer heatmap para cada biotype y tambien ver cuantos son
# differentially expressed
# hacer volcano plots
# empezar a hacer figura 1
# cuantos son exclusivos de LSK??
# agregar MPPs?

table(rownames(results_list$potNovelT.cell[results_list$potNovelT.cell$Diff=="DOWN",])%in%
        rownames(results_list$potNovelmacrophage[results_list$potNovelmacrophage$Diff=="DOWN",]))
table(rownames(results_list$lncRNAT.cell[results_list$lncRNAT.cell$Diff=="DOWN",])%in%
        rownames(results_list$lncRNAmacrophage[results_list$lncRNAmacrophage$Diff=="DOWN",]))

# 529 lncRNAs and 760 potential novel are enriched in LSK compared to macrophages and T-cells
#
res_LSKvsT.cell <- as.data.frame(results(dds, contrast = c("sample", "LSK","T.cell")))
res_LSKvsT.cell <- res_LSKvsT.cell %>% mutate(Diff=ifelse(padj<0.05,
                                                          ifelse(log2FoldChange<0,
                                                                 "DOWN","UP"),"NS"))

res_LSKvsmacrophage <- as.data.frame(results(dds, contrast = c("sample", "LSK","macrophage")))
res_LSKvsmacrophage <- res_LSKvsmacrophage %>% mutate(Diff=ifelse(padj<0.05,
                                                          ifelse(log2FoldChange<0,
                                                                 "DOWN","UP"),"NS"))

gene_level_info$LSK_vs_T.cell <- res_LSKvsT.cell$Diff[match(gene_level_info$gene_name,
                                                            rownames(res_LSKvsT.cell))]

gene_level_info$LSK_vs_macrophage <- res_LSKvsmacrophage$Diff[match(gene_level_info$gene_name,
                                                            rownames(res_LSKvsmacrophage))]

table(gene_level_info$LSK_vs_T.cell=="UP"&gene_level_info$LSK_vs_macrophage=="UP",
      gene_level_info$biotype)

View(gene_level_info%>%filter(LSK_vs_T.cell=="UP",LSK_vs_macrophage=="UP",biotype=="potNovel"))
# plot segments
ggplot(annotated_tracking%>%filter(V2=="STLOC_000372"), aes(x = start, xend = end, y = V1)) +
  geom_segment(color = "blue") +
  labs(x = "Position", y = "Chromosome") +
  ggtitle("Genomic Ranges") +
  theme_bw()

# I see many examples of pot novel genes that are not
# assembled in LSK data but however they appear as enriched in LSK
# compared to both samples when I perform kallisto and DESeq2
# What to do?
# I think I might clean up a bit the transcriptome before doing DGEA
# Advantage of having performed kallisto (transcript-level quantification)
# is that I can summarize at gene level in different ways:
#  use the same gene name for loci that have the same XLOC, especially for new genes
# 2) remove spurious classcodes: p, s, e and o
# 3) for known genes, keep only classcodes "=", "j", "k","m","n" (check my paper)
# 4) collapse the non-overlapping isoforms keeping the same gene id for a given XLOC
# 5) consider removing intronic things to PCGs in sense and also "p" classcodes (check my paper)

