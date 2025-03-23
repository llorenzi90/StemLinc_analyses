library(rtracklayer)
library(trastools)
# gtf_path="data/raw/T_cell.combined.gtf"
# gtf=readGFF(gtf_path)

tracking_path="data/raw/T_cell.tracking"
tracking <- read.table(tracking_path)
genes_in_samples <- trastools::get_sample_occurrence_per_gene_from_tracking(tracking)
colSums(genes_in_samples[,-4])
table(rowSums(genes_in_samples[,-4]))

colSums(tracking[,5:ncol(tracking)]!="-")


tracking_path="data/raw/Macro.tracking"
tracking <- read.table(tracking_path)
genes_in_samples <- trastools::get_sample_occurrence_per_gene_from_tracking(tracking)
colSums(genes_in_samples[,-4])
table(rowSums(genes_in_samples[,-4]))

colSums(tracking[,5:ncol(tracking)]!="-")


tracking_path="data/raw/LSK_StemLinc.tracking"
tracking <- read.table(tracking_path)
genes_in_samples <- trastools::get_sample_occurrence_per_gene_from_tracking(tracking)
colSums(genes_in_samples[,-4])
table(rowSums(genes_in_samples[,-4]))

colSums(tracking[,5:ncol(tracking)]!="-")
