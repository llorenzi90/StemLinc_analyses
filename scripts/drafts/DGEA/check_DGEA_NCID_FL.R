fc_FL= read.table('/run/user/1608803857/gvfs/smb-share:server=10.110.20.7,share=bdcuartero/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/featureCounts/PCGslncRNAsANDpotNovelLSK.counts.HPC7_FL.HSC.txt' ,header = T)
biots=read.table('/run/user/1608803857/gvfs/smb-share:server=10.110.20.7,share=bdcuartero/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/featureCounts/PCGslncRNAsANDpotNovelLSK.biotypes.txt' ,header = T)

counts <- fc_FL[,1:6]
colnames(counts) <- c("FL_1","FL_2","FL_3",
                           "FL.NICD_1","FL.NICD_2","FL.NICD_3")

library(DESeq2)
library(tidyverse)
coldata=data.frame(sample_id=colnames(counts),
                   sample=c(rep("FL",3),rep("FL.NICD",3)),
                   rep=rep(1:3,2))
rownames(coldata)=colnames(counts)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ sample)

table(rowSums(counts(dds)) > 10)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, contrast = c("sample", "FL.NICD","FL"))
summary(res)
res=as.data.frame(res)
res=res[!is.na(res$padj),]
res <- res %>% mutate(Diff=ifelse(padj<0.05,
                                                          ifelse(log2FoldChange<0,
                                                                 "DOWN","UP"),"NS"))
res$gene_id=rownames(res)
res <- left_join(res,biots%>%select(gene_id,biotype=simpl_biotype))
