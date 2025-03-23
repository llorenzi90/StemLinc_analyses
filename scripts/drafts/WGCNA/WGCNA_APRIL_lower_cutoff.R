#### ### ### ### ### ### ### ### ###
## 
## Purpose of script: WGCNA
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-02-05 - Updated 2024-04-13
##
## Email: lucialorenzi90@gmail.com
##
# Notes ---------------------------
##
## Script adapted from K Patel's tutorial
##  https://github.com/kpatel427/YouTubeTutorials/blob/main/WGCNA.R
##   
##
# Setup ---------------------------

options(scipen = 999) 
# Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
#devtools::install_github("kevinblighe/CorLevelPlot")
library(gridExtra)
library(sva)
library(ggpubr)
library(GOfuncR)
library(data.table)

# Load data---------------------------
setwd("/home/llorenzi/work_local/APRIL/")


setwd("WGCNA")
dir.create("WGCNA_lower_filter")
# setwd("~/Desktop/demo/WGCNA")


gene_translation_file <- "/home/llorenzi/references/conversion_table_ENSEMBL_NCBI_ids.csv"
gene_translation <- fread(gene_translation_file)


allowWGCNAThreads()          # allow multi-threading (optional)

# 1. Fetch Data ------------------------------------------------
data <- read.delim("/home/llorenzi/work_local/expression_files/featureCounts/PCGslncRNAsPseudogenes.110424.counts.txt")

gene_id_biotype=read.delim("/home/llorenzi/work_local/expression_files/featureCounts/PCGslncRNAsPseudogenes.110424.annot.txt")

table(gene_id_biotype$simpl_biotype)


# create metadata
sample_type=sapply(strsplit(colnames(data),split = "_"),function(x)x[[1]])
study=sapply(strsplit(colnames(data),split = "_"),function(x)x[[2]])
rep=sapply(strsplit(colnames(data),split = "_"),function(x)x[[3]])
coldata=data.frame(sample=colnames(data),
                   sample_type=sample_type,
                   study=study,
                   rep=rep)
rownames(coldata)=colnames(data)

#sort samples first by sample type

sort(sample_type)
sample_type_order=c("LT.HSC","LSK","MPP1","MPP2",
                    "MPP3","MPP3.4","MPP4","CMP",
                    "CLP","MEP","GMP","megakaryocyte",
              "erythrocyte","granulocyte","monocyte",
              "macrophage","dendritic","NK","PreB",
              "ProB","B.cell","T.cell")

table(sample_type%in%sample_type_order)
#sample_type[!sample_type%in%sample_type_order]#run this 
#and modify until sample_type is all in sample_type_order

coldata=coldata[order(match(sample_type,sample_type_order),
                      study,rep),]

data <- data[,match(coldata$sample,colnames(data))]

# prepare data
gene_ids=rownames(data)



# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

table(grepl("XLOC",gene_ids[!gsg$goodGenes]))

# remove genes that are detectd as outliers
total=table(gene_id_biotype$biotype)
pass=table(as.factor(gene_id_biotype$biotype)[gene_id_biotype$Geneid%in%
                                          rownames(data)[gsg$goodGenes]])
nopass=table(as.factor(gene_id_biotype$biotype)[gene_id_biotype$Geneid%in%
                                rownames(data)[!gsg$goodGenes]])



data_noMP <- data[,coldata$sample_type!="macrophage"]
coldata_noMP=coldata[coldata$sample_type!="macrophage",]

data_noMP=as.matrix(data_noMP)
coldata_noMP=coldata[coldata$sample_type!="macrophage",]

coldata <- coldata_noMP

# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

# exclude outlier samples

data_noMP_batchcorrected=read.table("knownANDpotNovellncRNAs.PCGs.pseudogenes.counts.batch_corrected.notextremenumbers.txt")

class(data_noMP_batchcorrected)
cn=colnames(data_noMP_batchcorrected)
rn=rownames(data_noMP_batchcorrected)
data_noMP_batchcorrected <- matrix(as.integer(as.matrix(data_noMP_batchcorrected)),
                                   nrow = dim(data_noMP_batchcorrected)[1])

colnames(data_noMP_batchcorrected)=cn
rownames(data_noMP_batchcorrected)=rn


# create DESeq2 data
dds <- DESeqDataSetFromMatrix(countData = data_noMP_batchcorrected,
                              colData = coldata,
                              design = ~ 1) # not spcifying model



## remove all genes with counts < 15 in more 
# than 75% of primary samples
primary_cells=grep("LT|LSK|MPP",coldata$sample_type,value = T)
length(primary_cells)
samps75=0.75*length(primary_cells)
## suggested by WGCNA on RNAseq FAQ

count_filter=rowSums(counts(dds) >= 15) >= samps75
count_filter=cbind(count_filter,gene_id_biotype[match(names(count_filter),
                                         gene_id_biotype$Geneid),c("biotype","simpl_biotype")])

setwd("WGCNA_lower_filter/")
write.table(count_filter,"Genes_pass_filter_15_counts_75percent_primarysamples.txt",col.names = F,quote = F,sep = "\t")

selected_genes_info=read.delim("/home/llorenzi/work_local/APRIL/selected_genes_based_on_corr.txt")
table(selected_genes_info$gene_id%in%
        rownames(count_filter)[count_filter$count_filter])
table(count_filter$count_filter,count_filter$biotype)

#815 potential novel pass this filter

# what about genes that pass the filter of 20 counts in StemLinc 
# # samples?
# StemLinc_samps=grep("LSK_StL",colnames(counts(dds)))
# StemLinc_filter=rowSums(counts(dds)[,StemLinc_samps] >= 20) ==3
# StemLinc_filter=cbind(StemLinc_filter,gene_id_biotype[match(names(StemLinc_filter),
#                                                       gene_id_biotype$Geneid),c("biotype","simpl_biotype")])
# 
# write.table(StemLinc_filter,"Genes_pass_filter_20_counts_in_all3_StemLinc_LSK.txt",col.names = F,quote = F,sep = "\t")
# table(StemLinc_filter$StemLinc_filter,StemLinc_filter$biotype)

dds75 <- dds[rowSums(counts(dds) >= 15) >= samps75,]
#dds75 <- dds[rowSums(counts(dds) >= 15) >= 61,]

nrow(dds75) # 19178 (13284 genes in tutorial)


# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

saveRDS(list(dds=dds,norm.counts=norm.counts),"dds_norm.counts.RDS")
mlist=readRDS("dds_norm.counts.RDS")

norm.counts=mlist$norm.counts
coldata_noMP=coldata[coldata$sample_type!="macrophage",]

coldata <- coldata_noMP
# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()
  

grid.arrange(a1, a2, nrow = 2)

tiff("scaleFreeTopology_Connectivity_vs_power.tiff",
     units = "in",res = 300,width = 10,height = 5)
grid.arrange(a1, a2, nrow = 2)
dev.off()


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                 maxBlockSize = 14000,
                 TOMType = "signed",
                 power = soft_power,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed = 1234,
                 verbose = 3)

cor <- temp_cor

saveRDS(bwnet,"bwnet.rds")
bwnet <- readRDS("bwnet.rds")


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], 
                    cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


tiff("Dendrogram_modules.tiff",width = 10,height = 5,res = 300,units = "in")
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)
dev.off()

# grey module = all genes that doesn't fall into other modules were assigned to the grey module


# 6 Relate modules to traits --------------------------------------------------
# module trait associations

# In my case I want to find genes or clusters of genes that
# have significant association with the more primary
# (less differentiated) blood cells

# I will start by defining primary cell as LT.HSC, LSK and MPPs
# create traits file - binarize categorical variables
primary_cells=grep("LT|LSK|MPP",coldata$sample_type,value = T)
coldata$sample_type[!coldata$sample_type%in%primary_cells]
traits <-  coldata %>% mutate(is_primary=ifelse(grepl("LT|LSK|MPP",
                                                           sample_type),1,0)) %>% 
  select(5)


# binarize categorical variables
# K Patel uses severity, I can use the type of cell instead


coldata$sample_type <- factor(coldata$sample_type, 
                           levels = unique(coldata$sample_type))

sample_type.out <- binarizeCategoricalColumns(coldata$sample_type,
                           includePairwise = FALSE,
                           includeLevelVsAll = TRUE,
                           dropFirstLevelVsAll = F, #I added this otherwise LT.HSC dissapears, not sure if this is correct, but for now I'll keep it
                           minCount = 1)


traits <- cbind(traits, sample_type.out)


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)



# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, 
                      traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')




CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[grep("^ME",names(heatmap.data),invert = T)],
             y = names(heatmap.data)[grep("^ME",names(heatmap.data))],
             rotLabX = 90,
             col = c("blue1", "skyblue", "white", "pink", "red"))


tiff("CorrLevelPlot.tiff",width = 15,height = 10,units = "in",res = 300)

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[grep("^ME",names(heatmap.data),invert = T)],
             y = names(heatmap.data)[grep("^ME",names(heatmap.data))],
             rotLabX = 90,
             col = c("blue1", "skyblue", "white", "pink", "red"))

dev.off()

module.gene.mapping <- as.data.frame(bwnet$colors)
rownames(module.gene.mapping) <- colnames(norm.counts)

# almost all potential novel genes are in grey module
module.gene.mapping$biotype=gene_id_biotype$simpl_biotype[match(rownames(module.gene.mapping),
                                                                    gene_id_biotype$Geneid)]

module.gene.mapping$gene_id=rownames(module.gene.mapping)

write.table(module.gene.mapping,"module.gene.mapping.txt",quote = F,row.names = F,sep = "\t")


                                                                             
table(module.gene.mapping$`bwnet$colors`,module.gene.mapping$biotype) 

## 6.1 NOTE ----

# This might be not the best way of setting the traits, 
# I would encourage to play around with it


# 7. GO enrichment with genes of each module ----

# GO enrichment with PCGs of brown module

gene_universe=colnames(norm.counts)
gene_universe <- gene_translation$gene_name[match(gene_universe,
                                                  gene_translation$gene_id)]
gene_universe <- na.omit(gene_universe)

go_enrichment_module <- function(col){
  mod_genes <- na.omit(gene_translation$gene_name[match(rownames(module.gene.mapping[module.gene.mapping$`bwnet$colors`==col,]),
                                                        gene_translation$gene_id)])
  
  input_hyper_mouse = data.frame(gene_id=gene_universe, 
                                 is_candidate=ifelse(gene_universe%in%mod_genes,1,0))
  
  
  res_hyper_mouse = go_enrich(input_hyper_mouse, orgDb = "org.Mm.eg.db")
  
  go_res=res_hyper_mouse$results
  
  
  anno_genes = get_anno_genes(go_ids=go_res$node_id, 
                              genes=mod_genes,database = "org.Mm.eg.db")
  
  anno_genes$GO_name=go_res$node_name[match(anno_genes$go_id,
                                            go_res$node_id)]
  
  anno_genes <- anno_genes %>% filter(!GO_name%in%c("molecular_function",
                                                    "cellular_component",
                                                    "biological_process"))
  
  
  #add info on genes associated with each GO term
  
  
  go_id_gene_ids=anno_genes %>% group_by(go_id) %>%
    summarise(genes=paste0(gene,collapse = ","),
              Ngenes=length(gene))
  
  
  go_res <- left_join(go_res,go_id_gene_ids,by=c("node_id" = "go_id"))
  
  return(go_res)
}

interesting_modules=c("yellow","lightcyan","grey60","lightyellow","cyan","turquoise","salmon","tan","magenta","pink","red",
                      "green","darkred")
go_res_interesting_modules=lapply(interesting_modules,go_enrichment_module)
names(go_res_interesting_modules) <- interesting_modules
View(go_res_interesting_modules$yellow)
View(go_res_interesting_modules$lightcyan)
View(go_res_interesting_modules$grey60)
View(go_res_interesting_modules$green)
View(go_res_interesting_modules$darkred)
saveRDS(go_res_interesting_modules,"go_res_interesting_modules_list.RDS")
colors2add=c("brown","blue")
go_res_colors2add=lapply(colors2add,go_enrichment_module)
names(go_res_colors2add)=c("brown","blue")
go_res_interesting_modules=readRDS("go_res_interesting_modules_list.RDS")
go_res_interesting_modules=c(go_res_interesting_modules,
                             go_res_colors2add)

saveRDS(go_res_interesting_modules,"go_res_interesting_modules_list.RDS")

# 8. Module-trait association ---------------

# plot module eigengenes vs is_primary

yellow.vs.is_primary <- data.frame(is_primary=as.factor(traits$is_primary),
                                 MEyellow=module_eigengenes$MEyellow)

g=ggplot(yellow.vs.is_primary,aes(x=is_primary, 
                              y=MEyellow,
                              col=is_primary)) + 
  geom_boxplot() +
  geom_point() +
  theme_classic()+ 
  stat_compare_means(method = "t.test")  

tiff("yellow_eigengenes_vs_is_primary_trait.boxplot.tiff",res=300,
     width = 5,height = 5,units = "in")
print(g)
dev.off()

yellow_genes=module.gene.mapping$gene_id[module.gene.mapping$`bwnet$colors`=="yellow"]

# same but with gene expression values

# NOTE: In this plot I am showing the mean expression of
# yellow genes for each sample
# I don't think this is very useful, it only shows that,
# in primary samples, the yellow genes are on average increased in expression
# relative to the non-primary samples,
# in other words, on average yellow genes are up-regulated in primary cells
# compared to non-primary cells. But it doesn't mean that ALL genes in the 
# yellow module are up-regulated, some of them are down-regulated, or 
# which is the same, negatively correlated with the trait
# Most genes possitively correlated though

# to see the trend per gene I have to calculate the mean expression 
# across primary and secondary samples

# 
yellow.vs.is_primary_exp <- data.frame(is_primary=
                                        as.factor(traits$is_primary),
                             MEyellow_mean_norm_counts=
                               apply(norm.counts[,colnames(norm.counts)%in%yellow_genes],1,mean))

g=ggplot(yellow.vs.is_primary_exp,aes(x=is_primary, 
                              y=MEyellow_mean_norm_counts,
                              col=is_primary)) + 
  geom_boxplot() +
  geom_point() +
  theme_classic()+ 
  stat_compare_means(method = "t.test")  

tiff("yellow_meanGeneexp_vs_is_primary_trait.boxplot.tiff",res=300,
     width = 5,height = 5,units = "in")
print(g)
dev.off()

# mean expression PER GENE
# how can I weight by the total number of samples in each
# sample-type?
# take random samples of equal number in each gene?
# I can pasar por alto this now because 
# in this case the samples are not that unbalanced
# 38 vs 46
meanyellowgeneExp_perisprimary <- data.frame(is_primary=c(
                                        rep(0,length(yellow_genes)),
                                        rep(1,length(yellow_genes))),
                                      meanMEyellow=c(
                                        apply(norm.counts[
                                          rownames(norm.counts)%in%
                                            rownames(traits[
                                          traits$is_primary==0,]),
                                          colnames(norm.counts)%in%
                                            yellow_genes],
                                          2,mean),
                                          apply(norm.counts[
                                            rownames(norm.counts)%in%
                                              rownames(traits[
                                                traits$is_primary==1,]),
                                            colnames(norm.counts)%in%
                                              yellow_genes],
                                            2,mean)),
                                      gene=rep(colnames(norm.counts[,colnames(norm.counts)%in%yellow_genes]),2)
                                      )


meanyellowgeneExp_perisprimary$biotype=gene_id_biotype$simpl_biotype[
  match(meanyellowgeneExp_perisprimary$gene,
        gene_id_biotype$gene_id)
]
g=ggplot(meanyellowgeneExp_perisprimary,aes(x=as.factor(is_primary), 
                                     y=meanMEyellow,
                                     col=as.factor(is_primary))) + 
  geom_boxplot() +
  geom_point() +
  theme_classic()+ 
  stat_compare_means(method = "t.test",paired = T)  

print(g)

tiff("meanBronwGeneExp_across_primary_vs_non_primary.boxplot.tiff",res=300,
     width = 10,height = 7,units = "in")
print(g)
dev.off()

# same but plot FC
FC_primary_non_primary=meanyellowgeneExp_perisprimary %>% group_by(gene) %>%
  summarise(log2FC=log2(meanMEyellow[is_primary==1]/meanMEyellow[is_primary==0]),
            FC=meanMEyellow[is_primary==1]/meanMEyellow[is_primary==0])

FC_primary_non_primary$biotype=gene_id_biotype$biotype[match(FC_primary_non_primary$gene,
                                                             gene_id_biotype$Geneid)]

g=ggplot(FC_primary_non_primary,aes(y=log2FC,fill=biotype,x=biotype)) + 
  geom_boxplot() + geom_jitter() + theme_classic() + ggtitle("log2 FC primary vs non-primary cells yellow module genes")

tiff("log2FC_yellowGenes_primary_vs_nonprimary.boxplot.perbiotype.tiff",
     width = 5,height = 5,res = 300,units = "in")
print(g)
dev.off()

g=ggplot(FC_primary_non_primary,aes(y=FC,fill=biotype,x=biotype)) + 
  geom_boxplot() + geom_jitter() + theme_classic() + ggtitle("FC primary vs non-primary cells yellow module genes")

tiff("FC_yellowGenes_primary_vs_nonprimary.boxplot.perbiotype.tiff",
     width = 5,height = 5,res = 300,units = "in")
print(g)
dev.off()

# 9. Identify potential driver genes --------------------------------------

# Identify those genes within a module that are
# 1) Highly connected within the module (hub genes)
# AND
# 2) Most strongly correlated with a clinical/phenotypical trait of interest

# Module membership: Correlation of a gene to a module eigengene
# • Genes with high module membership are good representatives of the
# overall expression profile in the module
# • Genes with high module membership tend to be “hub” genes in the
# module (high intramodule connectivity)
# • A gene can have high membership in several modules
# (not just the one to which it is assigned)

# Calculate the module membership and the associated p-values

all(rownames(module_eigengenes)==rownames(norm.counts))
module.membership.measure <- 
  cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- 
  corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure=as.data.frame(module.membership.measure)
module.membership.measure$module=rownames(module.membership.measure)
module.membership <- pivot_longer(module.membership.measure,
                                  cols = 1:(ncol(module.membership.measure)-1),
                                  names_to = "gene",
                                  values_to = "module_membership")

module.membership.measure.pvals=as.data.frame(module.membership.measure.pvals)
module.membership.measure.pvals$module=rownames(module.membership.measure.pvals)
module.membership.pvals <- pivot_longer(module.membership.measure.pvals,
                                  cols = 1:(ncol(module.membership.measure.pvals)-1),
                                  names_to = "gene",
                                  values_to = "module_membership_pval")

module.membership <- left_join(module.membership,module.membership.pvals)
module.membership <- left_join(module.membership,module.gene.mapping,by=c(gene="gene_id"))

colnames(module.membership)[5]="assigned_module"


module.membership$gene_name=gene_translation$gene_name[match(module.membership$gene,
                                                             gene_translation$gene_id_vM31)]
module.membership$gene_name=ifelse(is.na(module.membership$gene_name),module.membership$gene,module.membership$gene_name)

View(module.membership%>%filter(assigned_module=="yellow",module=="MEyellow"))

write.table(module.membership,"module.membership.all.txt",sep = "\t",quote = F,row.names = F)

mod="yellow"

g=ggplot(module.membership%>%filter(assigned_module==mod,module==paste0("ME",mod)),
       aes(x=biotype, 
           y=module_membership)) +
  geom_boxplot() + geom_jitter() + theme_classic() + 
  ggtitle(paste("Module",mod))

print(g)

tiff(paste0("module_",mod,".module_membership.biotype.boxplot.tiff"),
     units = "in",res = 300,width = 10,height = 5)
print(g)
dev.off()


# Calculate the gene-trait correlation and associated p-values

gene.signf.corr <- cor(norm.counts, 
                       traits$is_primary, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, 
                                          nSamples)

gene.signf.corr=as.data.frame(gene.signf.corr)
gene.signf.corr$pval=gene.signf.corr.pvals[,1]
colnames(gene.signf.corr)[1]="is_primary_corr"

gene.signf.corr$gene=rownames(gene.signf.corr)
gene.signf.corr= left_join(gene.signf.corr,
                           module.membership
                           %>% filter(module=="MEyellow") 
                           %>% dplyr::select(gene,module_membership,module_membership_pval,assigned_module,biotype))


colnames(gene.signf.corr)[4:5]=c("yellow_module_membership",
                                 "yellow_module_memb_pval")

assigned_module_membership=module.membership%>% filter(module==paste0("ME",assigned_module))
colnames(assigned_module_membership)[3:4]=paste0("assigned_",colnames(assigned_module_membership)[3:4])

gene.signf.corr=left_join(gene.signf.corr,assigned_module_membership %>% dplyr::select(gene,assigned_module_membership,assigned_module_membership_pval))
gene.signf.corr$gene_name=gene_translation$gene_name[match(gene.signf.corr$gene,
                                                           gene_translation$gene_id_vM31)]

gene.signf.corr$gene_name=ifelse(is.na(gene.signf.corr$gene_name),gene.signf.corr$gene,gene.signf.corr$gene_name)

write.table(gene.signf.corr,"isprimaryCorrANDyellowmembership.txt",sep = "\t",row.names = F,quote = F)


g=ggplot(gene.signf.corr,aes(x=biotype,y=is_primary_corr,fill=biotype)) +
  geom_boxplot() + geom_jitter()+theme_classic() + 
  geom_hline(yintercept = 0.6)

tiff("Corr_isprimary_vs_all_biotypes.boxplot.tiff",
     units = "in",res = 300,width = 10,height = 5)
print(g)
dev.off()


g=ggplot(gene.signf.corr,aes(x=yellow_module_membership ,
                           y=is_primary_corr,
                           col=assigned_module,
                           )
       ) + geom_point()+ 
theme_classic() +
  scale_color_manual(values = sort(unique(gene.signf.corr$assigned_module)))


tiff("Corr_isprimary_vs_yellow_membership.all_genes.tiff",
     units = "in",res = 300,width = 10,height = 5)
print(g)
dev.off()

g=ggplot(gene.signf.corr%>% filter(assigned_module=="yellow"),
         aes(x=yellow_module_membership ,
                             y=is_primary_corr,
                             col=biotype,
)
) + geom_point()+
  theme_classic() 

tiff("Corr_isprimary_vs_yellow_membership.yellowmodule.tiff",
     units = "in",res = 300,width = 10,height = 5)
print(g)
dev.off()





# 10. Evaluate module quality ----------------------------------------------

# Are modules better than random groupings of genes?
# a.Connectivity
# mean intra-module connectivity
# mean ratio of intra-module / total connectivity
# b.Trait correlations
# strong correlation between module eigengenes and traits of interest
# strong correlation between gene module membership and gene-trait corr.
# c.Functional enrichment
# many functionally related genes in the same module
# GO enrichment analysis

intramodular.connectivity=intramodularConnectivity.fromExpr(datExpr = norm.counts,colors = module.gene.mapping$`bwnet$colors`)
intramodular.connectivity$gene=colnames(norm.counts)
intramodular.connectivity=left_join(intramodular.connectivity,gene.signf.corr)
write.table(intramodular.connectivity,"intramodular_connectivity.txt",sep = "\t",row.names = F,quote = F)

gene_level_info_potnovel=read.delim("/home/llorenzi/work_local/APRIL/selected_genes_all.27.txt")
View(intramodular.connectivity%>%filter(gene%in%gene_level_info_potnovel$gene_id))