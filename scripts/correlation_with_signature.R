options(scipen = 999)
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
## Load data---------------------------
source("scripts/source_all_functions.R")
mycolors=nejm_pal
source("scripts/load_gene_level_data_biotypes_of_interest.R")

# define paths ----
HSC_signature_path="data/gene_lists/HSC_signature3.txt"
vst_path="outputs/expression_data/featureCounts_20240930_173216.ordered_samples.selected_biotypes.vst.txt"
counts_path="outputs/expression_data/featureCounts_20240930_173216.ordered_samples.selected_biotypes.counts.txt"
coldata_path="outputs/expression_data/coldata_20240930_173216.ordered_samples.txt"

# define outdir
outdir="outputs/plots"

HSC_signature=read.table(HSC_signature_path)
table(HSC_signature$V1%in%gene_level_info$gene_name)

# Correlation with HSC signature ----
vst=read.table(vst_path,header = T)
all(vst$Geneid==gene_level_info$gene_name) # this has to be TRUE!

# remove gene id column
rownames(vst) <- vst$Geneid
vst <- vst[,-1]

# select non-protein_coding biotypes
vst_lncRNAs=vst[rownames(vst)%in%gene_level_info$gene_name[gene_level_info$biotype!="protein_coding"],]

# select HSC signature
vst_signature=vst[rownames(vst)%in%HSC_signature$V1,]

# calculate correlation between non-coding genes and HSC signature
corr=cor(t(vst_lncRNAs),t(vst_signature))

corr=as.data.frame(corr)
corr$ncGene=rownames(corr)
corr=pivot_longer(corr,cols=1:(ncol(corr)-1),
                  values_to = "corr",names_to = "PCG_name")

mean_corr=corr%>%group_by(ncGene)%>%summarise(mean_corr=mean(corr))

# correlation signature with PCGs
vst_PCG <- vst[rownames(vst)%in%gene_level_info$gene_name[gene_level_info$biotype=="protein_coding"],]

corr_PCG <- cor(t(vst_PCG),t(vst_signature))

corr_PCG=as.data.frame(corr_PCG)
corr_PCG$PCG=rownames(corr_PCG)
corr_PCG=pivot_longer(corr_PCG,cols=1:(ncol(corr_PCG)-1),
                      values_to = "corr",names_to = "PCG_name")


mean_corr_PCG=corr_PCG%>%group_by(PCG)%>%summarise(mean_corr=mean(corr))

colnames(mean_corr)[1]="gene_name"
colnames(mean_corr_PCG)[1]="gene_name"

mean_corr=rbind(mean_corr,mean_corr_PCG)
colnames(mean_corr)[2]="mean_corr_signature"
#gene_level_info <- left_join(gene_level_info,mean_corr)

#write_info_table(gene_level_info,gene_level_data_path,modify_path = F)

# collapse both corr table:
colnames(corr) <- c("gene_name","Signature_gene","corr")
colnames(corr_PCG) <- c("gene_name","Signature_gene","corr")

corr=rbind(corr,corr_PCG)
max_corr_signature <- corr %>% arrange(desc(corr))
max_corr_signature <- max_corr_signature %>% filter(Signature_gene!=gene_name)
max_corr_signature <- max_corr_signature[!duplicated(max_corr_signature$gene_name),]
colnames(max_corr_signature)[2:3] <- c("max_Signature_gene","max_corr_signature")
# compute correlation with any PCG

# corr_nc_PCG=cor(t(vst_lncRNAs),t(vst_PCG))
# corr_nc_PCG=as.data.frame(corr_nc_PCG)
# corr_nc_PCG$gene_name=rownames(corr_nc_PCG)
# corr_nc_PCG=pivot_longer(corr_nc_PCG,cols=1:(ncol(corr_nc_PCG)-1),
#                       values_to = "corr",names_to = "PCG_name")
# corr_nc_PCG <- corr_nc_PCG %>% arrange(desc(corr))
#
# corr_nc_PCG <- corr_nc_PCG[!duplicated(corr_nc_PCG$gene_name),]
#dir.create("outputs/correlation_analysis")
#write.table(corr_nc_PCG,"outputs/correlation_analysis/corr_ncRNAs_vs_allPCG.txt",quote = F,row.names = F,sep = "\t")
