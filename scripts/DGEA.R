source("~/Rprojects/StemLinc_analyses/scripts/correlation_with_signature.R")

coldata <- read.table(coldata_path,header = T)

## DGEA ----
library(DESeq2)
counts=read.table(counts_path,header = T)
all(counts$Geneid==gene_level_info$gene_name) # this has to be TRUE!

# remove gene id column
rownames(counts) <- counts$Geneid
counts <- counts[,-1]

rownames(coldata)=coldata$sample


dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design= ~ study + cell_class)

#dds <- DESeq(dds)
tstamp=get_timestamp(counts_path)
dir.create("outputs/DESeq2")

#saveRDS(dds,paste0("outputs/DESeq2/dds_DESeq2_all_blood_cells.",
#                  tstamp,"RDS"))
