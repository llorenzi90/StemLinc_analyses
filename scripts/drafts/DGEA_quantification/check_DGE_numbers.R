setwd("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation")

DGEA=read.delim("DGEA/filtered_genes.CPAT.CAGE.DGEA.tsv")

table(DGEA$LSK_vs_T.cell=="UP"&DGEA$LSK_vs_macrophage=="UP",
      DGEA$biotype,DGEA$exonic_type)

counts <- read.table("expression/featureCounts/all_genes.counts.tsv")
gene_annot <- read.delim("DGEA/filtered_genes.CPAT.CAGE.DGEA.tsv")

lncRNA_counts <- counts[rownames(counts)%in%gene_annot$gene_id[gene_annot$biotype=="lncRNA"],]

table(rowSums(lncRNA_counts>10)>2)
table(rowSums(lncRNA_counts>30)>2)
low_lncRNAs <- rownames(lncRNA_counts)[rowSums(lncRNA_counts>10)<2]

DGEA %>% filter(gene_id%in%low_lncRNAs) %>% summarise(LSK_vs_T.cell=sum(LSK_vs_T.cell!="NS",na.rm = T),
                                                      LSK_vs_macrophage=sum(LSK_vs_macrophage!="NS",na.rm = T),
                                                      T.cell_vs_macrophage=sum(T.cell_vs_macrophage!="NS",na.rm = T))

# DGE distribution per biotype
table(DGEA$biotype,DGEA$exonic_type,DGEA$LSK_vs_macrophage!="NS"|DGEA$LSK_vs_macrophage!="NS"|DGEA$T.cell_vs_macrophage!="NS")

# What if I apply a minimum expression cutoff?
counts_pass <- counts[rowSums(counts>30)>2,]

DGEA_pass <- DGEA %>%filter(gene_id%in%rownames(counts_pass))

table(DGEA_pass$biotype,
      DGEA_pass$exonic_type,
      DGEA_pass$LSK_vs_macrophage!="NS"|DGEA_pass$LSK_vs_macrophage!="NS"|DGEA_pass$T.cell_vs_macrophage!="NS")
