# Collect all data and make summary heatmap
library(DESeq2)

source("scripts/source_all_functions.R")
source("scripts/load_gene_level_data_biotypes_of_interest.R")

# data to add
all_marks_path <- "outputs/overlap_marks/gene_level_info.20240930_173216.logicmarks.tsv"
deseq2_res_path <- "outputs/DESeq2/dds_DESeq2_all_blood_cells.20240930_173216.RDS"
TF_data_path <- "outputs/overlap_marks/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.TFs_rel_level_500_250.tsv"
TF_data_long_path <- "outputs/TF_enrichment/TF_data_long.txt"
Delas_Luo_overlap_path <- "outputs/overlap_lncRNA_catalogs/overlap_StemLinc_genes_Delas_and_Luo_genes.txt"
expression_data_path <- "outputs/expression_data/featureCounts_20240930_173216.ordered_samples.selected_biotypes.TPM.txt"
closest_PCG_path <- "outputs/closest_PCGs.txt"
corr_allPCGs_path <- "outputs/correlation_analysis/corr_ncRNAs_vs_allPCG.txt"
# collect all data  ----

# max corr with signature ----
source("scripts/correlation_with_signature.R")
colnames(gene_level_info)
gene_level_info <- left_join(gene_level_info,max_corr_signature)

corr_allPCGs <- read.delim(corr_allPCGs_path)
colnames(corr_allPCGs)[2:3]=c("maxcorrPCG","max_corr")

gene_level_info <- left_join(gene_level_info,corr_allPCGs)

# chromatin marks ----
all_marks <- read.table(all_marks_path,header = T)
# NOTE these logic marks were originally generated with script "create_logical_marks_and_plots.R"
gene_level_info <- left_join(gene_level_info,all_marks%>% dplyr::select(gene_name,evidence_level,Enhancer,
                                                                        CAGE_within_100bp,
                                                                        polyAsite_within_100bp,
                                                                        H3K4me3_at_promoter,
                                                                        H3K4me1_at_promoter,
                                                                        H3K27ac_at_promoter,
                                                                        H3K36me3_at_geneBody,
                                                                        Enhancer_Atlas_at_promoter,
                                                                        FANTOM_enhancers_at_promoter,
                                                                        Enhancer_from_marks_at_promoter,
                                                                        H3K27me3_at_promoter))
# DGEA results ----
dds <- readRDS(deseq2_res_path)
contrast1=c("HSPC","differentiated")
contrast2=c("HSPC","progenitor")
contrast3=c("progenitor","differentiated")
clist <- list("HSPC_vs_diff"=contrast1,
              "HSPC_vs_prog"=contrast2,
              "prog_vs_diff"=contrast3)

FCco=2
for (cname in names(clist)) {
  cont=clist[[cname]]
  res <- results(dds,contrast = c("cell_class",
                                  cont))
  res <- as.data.frame(res) %>%
    mutate(diff=ifelse(padj<0.05&!is.na(padj),
                       ifelse(log2FoldChange>FCco,"UP","DOWN"),"NS"))

  res$gene_name <- rownames(res)
  res <- res%>%dplyr::select(gene_name,log2FoldChange,diff)
  colnames(res)[2:3] <- c(paste0(cname,".logFC"),cname)
  gene_level_info <- left_join(gene_level_info,res)
}

# TF binding ----
TF_data_long <- read.table(TF_data_long_path,header = T)
TF_data <- read.table(TF_data_path,header = T)

gene_level_info <- left_join(gene_level_info , TF_data %>% select(gene_name,TF_rel.lev_A))

# overlap Delas/Luo ----
Delas_Luo_overlap <- read.delim(Delas_Luo_overlap_path,header = T)
gene_level_info <- left_join(gene_level_info,Delas_Luo_overlap %>%select(gene_name,Delas_best_gene,Delas_best_cc,lnc_name,orientation,
                                                                         Delas_enrichment,best_Luo_name,best_Luo_enrichment))


# distance to closest PCG ----
closest_class_to_PCG <- read.delim(closest_PCG_path)
query_target <- paste0(gene_level_info$gene_name,"_",gene_level_info$best_closest_PCG)

gene_level_info$dist_PCG <- closest_class_to_PCG$dist[match(query_target,closest_class_to_PCG$query_target)]

# non-coding genes ----
gene_level_info_nc <- gene_level_info %>%filter(biotype!="protein_coding")
gene_level_info_nc %>% group_by(biotype) %>% summarise(Nhigh_corr=sum(mean_corr_signature>0.65))
gene_level_info_nc %>% group_by(classif) %>% summarise(Nhigh_corr=sum(mean_corr_signature>0.65))


# N TFs ----
TF_long=read.delim("outputs/TF_enrichment/TF_data_long.txt")

# N TFs at promoters all genes ----
N_TFs=TF_long %>% group_by(gene_name) %>% summarise(N_TFs=length(unique(TF)))
gene_level_info <- left_join(gene_level_info , N_TFs)
gene_level_info$N_TFs[is.na(gene_level_info$N_TFs)]=0


# write gene level data with all evidence ----
write.csv(gene_level_info,"outputs/gene_level_info_all_evidence_oct24.csv",
          row.names = F)
# select for highly corr nc genes ----
gene_level_info_sel <- gene_level_info %>%
  filter(biotype!="protein_coding"&mean_corr_signature>0.65) %>% arrange(desc(mean_corr_signature))
# assign names to previously validated lncRNAs ----
gene_level_info_sel$new_gene_name=gene_level_info_sel$gene_name
gene_level_info_sel$new_gene_name[!is.na(gene_level_info_sel$lnc_name)&gene_level_info_sel$lnc_name=="Sphed"]="Sphed"
gene_level_info_sel$new_gene_name[!is.na(gene_level_info_sel$best_Luo_name)&gene_level_info_sel$best_Luo_name=="LncHSC-2"]="LncHSC-2"
# write selected 114 ncRNAs based on high corr ----
write.csv(gene_level_info_sel,"outputs/gene_level_info_all_evidence_oct24.selected114.csv",row.names = F)
