source("scripts/source_all_functions.R")
library(ggplot2)

gene_level_info <- read.csv("outputs/gene_level_info_all_evidence_oct24.csv")
mean_corr_co=0.65

gene_level_info_sel <- gene_level_info %>% filter(biotype!="protein_coding"&mean_corr_signature>mean_corr_co)



TF_long=read.delim("outputs/TF_enrichment/TF_data_long.txt")

TF_long <- TF_long %>%filter(gene_name%in%gene_level_info_sel$gene_name)

TF_long=read.delim("outputs/TF_enrichment/TF_data_long.txt")

TF_long <- TF_long %>%filter(gene_name%in%gene_level_info_sel$gene_name)

N_TFs <- TF_long %>% group_by(gene_name) %>% summarise(N_TFs=length(unique(TF)))

gene_level_info_sel <- left_join(gene_level_info_sel, N_TFs)


gene_level_info_sel$gene_name[!is.na(gene_level_info_sel$lnc_name)] <- "Spehd"
gene_level_info_sel$gene_name[gene_level_info_sel$best_Luo_name=="LncHSC-2"] <- "LncHSC-2"


g=gene_level_info_sel %>% ggplot(aes(classif,fill = exonic_type)) + geom_bar() +
  facet_wrap(~ biotype, scales = "free_x") + scale_fill_manual(values = nejm_pal)+theme(axis.text.x = element_text(angle = 45, hjust = 1))
g

save_tiff_svg(g,outdir = "outputs/plots/",filename = "exonic_type_classif_facet_biotype_selected")

N_bcet=gene_level_info_sel %>% group_by(biotype,classif ,exonic_type) %>% summarise(n())


# parse expression data ----
vst <- read.delim("outputs/expression_data/featureCounts_20240930_173216.ordered_samples.selected_biotypes.vst.txt")
coldata <- read.delim("outputs/expression_data/coldata_20240930_173216.ordered_samples.txt")

# divide samples in LT-HSC, LSK, MPPs, CMP, CLP, lymphoid and myeloid

coldata$new_sample_type <- coldata$sample_type
coldata$new_sample_type[grepl("MPP",coldata$sample_type)] <- "MPPs"
coldata$new_sample_type[61:76]
# [1] "MEP"           "MEP"           "MEP"           "GMP"           "GMP"           "GMP"
# [7] "megakaryocyte" "erythrocyte"   "granulocyte"   "granulocyte"   "granulocyte"   "granulocyte"
# [13] "monocyte"      "macrophage"    "macrophage"    "macrophage"
coldata$new_sample_type[61:76] <- "myeloid"
coldata$new_sample_type[77:87]
# [1] "dendritic" "NK"        "PreB"      "PreB"      "ProB"      "ProB"      "B.cell"    "B.cell"    "T.cell"
# [10] "T.cell"    "T.cell"
coldata$new_sample_type[77:87] <- "lymphoid"

vst_LSK <- vst[,c(1,25:27)]
vst_LSK <- pivot_longer(vst_LSK,cols = 2:4, values_to = "vst", names_to = "rep")
colnames(vst_LSK)[1]="gene_name"
vst_LSK <- left_join(vst_LSK,gene_level_info %>% select(gene_name,biotype))
ggplot(vst_LSK,aes(biotype,vst,fill=rep)) + geom_boxplot()

selected_samples=coldata$sample[c(1:)]
