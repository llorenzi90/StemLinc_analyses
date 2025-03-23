library(tidyverse)
setwd("/home/llorenzi/Rprojects/StemLinc_analyses/outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation")

DGEA=read.delim("DGEA/all_merged_genes_annotation.DGEA_Diff.tsv")

table(DGEA$LSK_vs_T.cell=="UP"&DGEA$LSK_vs_macrophage=="UP",
      DGEA$biotype,DGEA$exonic_type)

filtered_assembled_loci <- read.delim("filtered_assembled_loci_annotated_and_novel.tsv")
counts <- read.table("expression/featureCounts/all_genes.counts.tsv")
TPM_data <- read.table("expression/featureCounts/all_genes.TPM.tsv")
gene_annot <- read.delim("all_merged_genes_annotation.sel_biotypes.pass_filter_col.tsv")

lncRNA_counts <- counts[rownames(counts)%in%gene_annot$gene_id[gene_annot$biotype=="lncRNA"],]

table(rowSums(lncRNA_counts>10)>2)
table(rowSums(lncRNA_counts>30)>2)
low_lncRNAs <- rownames(lncRNA_counts)[rowSums(lncRNA_counts>10)<2]

TPM_data <- TPM_data[match(gene_annot$gene_id,rownames(TPM_data)),]
table(rownames(TPM_data)==gene_annot$gene_id)

table(rowSums(TPM_data>=0.2)>2)

biots <- c("protein_coding", "lncRNA",  "potNovel",       "pseudogene",     "TEC")

table(gene_annot$biotype%in%biots)
table(gene_annot$biotype[gene_annot$biotype%in%biots])
table(gene_annot$biotype[gene_annot$biotype%in%biots&rowSums(TPM_data>=0.2)>2])

# check differences between assembled and "expressed" lncRNAs ----
assembled_lncRNAs <- filtered_assembled_loci$gene_id[filtered_assembled_loci$gene_biotype=="lncRNA"]
expressed_lncRNAs <- gene_annot$gene_id[gene_annot$biotype=="lncRNA"&rowSums(TPM_data>=0.2)>2]
table(assembled_lncRNAs%in%expressed_lncRNAs)
table(expressed_lncRNAs%in%assembled_lncRNAs)

expresesd_but_not_assembled <- expressed_lncRNAs[!expressed_lncRNAs%in%assembled_lncRNAs]
assembled_but_not_expressed <- assembled_lncRNAs[!assembled_lncRNAs%in%expressed_lncRNAs]

# check differences between assembled and expressed genes ----
assembled_genes <- filtered_assembled_loci$gene_id
expressed_genes <- gene_annot$gene_id[rowSums(TPM_data>=0.2)>2]

table(assembled_genes%in%expressed_genes)
table(expressed_genes%in%assembled_genes)

gene_annot_expressed <- gene_annot[rowSums(TPM_data>=0.2)>2,]
gene_annot_expressed$assembled <- gene_annot_expressed$gene_id%in%assembled_genes
gene_annot_expressed <- gene_annot_expressed %>% select(-pass_filter)
gene_annot_expressed$assembled[gene_annot_expressed$biotype=="potNovel"] <- TRUE
write.table(gene_annot_expressed,"expressed_merged_genes_annotation.sel_biotypes.assembled_col.tsv",
            quote = F,row.names = F,sep = "\t")
table(gene_annot_expressed$biotype[gene_annot_expressed$biotype!="potNovel"],gene_annot_expressed$assembled[gene_annot_expressed$biotype!="potNovel"])

# keep only annotated genes if they are assembled in our transcriptome ----
gene_annot <- gene_annot_expressed

gene_annot <- gene_annot[gene_annot$biotype=="potNovel"|gene_annot$assembled,]

table(gene_annot$exonic_type,gene_annot$biotype)

# Filter out novel genes with high coding potential ----
nejm_pal=ggsci::pal_nejm()(8)

ggplot(gene_annot,aes(x=Coding_prob,col=biotype)) + geom_density(linewidth = 2) +
  theme_minimal() +scale_color_manual(values = nejm_pal) + theme(text = element_text(size = 18))

table(gene_annot$Coding_prob>0.44,gene_annot$biotype)

ggplot(gene_annot,aes(x=biotype,fill=Coding_prob>0.44)) + geom_bar() +
  theme_minimal() +scale_fill_manual(values = nejm_pal) +
  theme(text = element_text(size = 18),axis.text.x = element_text( angle = 45,hjust = 1))

gene_annot <- gene_annot[!(gene_annot$Coding_prob>0.44&gene_annot$biotype=="potNovel"),]

table(gene_annot$Coding_prob>0.44, gene_annot$biotype)
View(gene_annot%>%filter(Coding_prob>0.44, biotype=="lncRNA"))


CPAT_data <- read.delim("CPAT_results/genes_best_ORF_matched_strand.txt")
CPAT_data <- CPAT_data[!duplicated(CPAT_data$gene_name),]

table(gene_annot$gene_id%in%CPAT_data$gene_name)
colnames(CPAT_data)[colnames(CPAT_data)=="gene_name"] <- "gene_id"
table(gene_annot$gene_id%in%CPAT_data$gene_id)
gene_annot <- left_join(gene_annot, CPAT_data%>%dplyr::select(gene_id,Coding_prob))
gene_annot <- gene_annot%>%filter(biotype%in%biots)
gene_annot$Coding_prob[is.na(gene_annot$Coding_prob)] <- 0
table(gene_annot$biotype[gene_annot$biotype%in%biots&rowSums(TPM_data>=0.2)>2])

genes2keep <- gene_annot %>% filter(biotype=="protein_coding"|Coding_prob<0.44)
table(gene_annot$biotype[gene_annot$gene_id%in%genes2keep$gene_id&rowSums(TPM_data>=0.2)>2])

biots_TPM <- TPM_data[rownames(TPM_data)%in%gene_annot$gene_id,]
table(rownames(biots_TPM)==gene_annot$gene_id)

table(genes2keep$biotype)
table(gene_annot$biotype[gene_annot$gene_id%in%genes2keep$gene_id&rowSums(biots_TPM>=0.2)>2])
table(gene_annot$biotype[gene_annot$gene_id%in%genes2keep$gene_id&rowSums(biots_TPM>=0.2)>2],
      gene_annot$exonic_type[gene_annot$gene_id%in%genes2keep$gene_id&rowSums(biots_TPM>=0.2)>2])


gene_ids_2_keep <- gene_annot$gene_id[gene_annot$gene_id%in%genes2keep$gene_id&rowSums(biots_TPM>=0.2)>2]

DGEA <- DGEA %>% filter(gene_id%in%gene_ids_2_keep)
DGEA %>% summarise(LSK_vs_T.cell=sum(LSK_vs_T.cell!="NS",na.rm = T),
                                                      LSK_vs_macrophage=sum(LSK_vs_macrophage!="NS",na.rm = T),
                                                      T.cell_vs_macrophage=sum(T.cell_vs_macrophage!="NS",na.rm = T))

DGEA <- DGEA %>% mutate(HSPC_enriched=ifelse(!is.na(LSK_vs_T.cell)&LSK_vs_T.cell=="UP",
                                             ifelse(!is.na(LSK_vs_macrophage)&LSK_vs_macrophage=="UP",
                                                    "both",
                                                    "T.cell"),
                                             ifelse(!is.na(LSK_vs_macrophage)&LSK_vs_macrophage=="UP",
                                                    "macrophage",
                                                    "none")))
# DGE distribution per biotype
table(DGEA$biotype,DGEA$HSPC_enriched)
table(DGEA$biotype,DGEA$exonic_type,DGEA$HSPC_enriched)
