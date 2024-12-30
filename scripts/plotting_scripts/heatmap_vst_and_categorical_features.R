gene_level_info <- read.csv("outputs/gene_level_info_all_evidence_oct24.csv")
mean_corr_co=0.65

gene_level_info_sel <- gene_level_info %>% filter(biotype!="protein_coding"&mean_corr_signature>mean_corr_co)

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

vst_selected <- vst %>% filter(Geneid %in% gene_level_info_sel$gene_name)

vst_selected_long <- pivot_longer(vst_selected,cols = 2:ncol(vst_selected),names_to = "sample",values_to = "vst")
vst_selected_long <- left_join(vst_selected_long , coldata %>% select(sample,new_sample_type))

vst_selected_mean_vst_per_sample_type <- vst_selected_long %>% group_by(Geneid,new_sample_type) %>%
  summarise(mean_vst=mean(vst))

vst_selected_wide <- pivot_wider(vst_selected_mean_vst_per_sample_type,id_cols = 1,names_from = 2,values_from = 3)
gene_name=vst_selected_wide$Geneid
vst_selected_wide <- vst_selected_wide[,-1]
vst_selected_wide <- vst_selected_wide[,match(unique(coldata$new_sample_type),
                                              colnames(vst_selected_wide))]

vst_selected_wide <- as.matrix(vst_selected_wide)
rownames(vst_selected_wide) <- gene_name


# pheatmap expression ----
library(pheatmap)
outdir="outputs/plots/"
p=pheatmap(vst_selected_wide,show_rownames = F)
p
save_tiff_svg(p,outdir = outdir,filename = "heatmap_114_selected_candidates_new_sample_types")


# add annotation ----
gene_level_info_sel <- gene_level_info_sel %>% mutate(in_Delas_Luo=!is.na(Delas_best_gene)|!is.na(best_Luo_name))


gene_level_info_sel <- gene_level_info_sel %>% mutate(Annotation=ifelse(biotype=="potNovel"&in_Delas_Luo,
                                                                        "also_Delas/Luo",biotype))


selected_features <- gene_level_info_sel %>%dplyr::select(gene_name,
                                                      exonic_type,
                                                      mean_corr_signature,
                                                      evidence_level,
                                                      Enhancer,
                                                      Annotation,
                                                      classif)


selected_features <- selected_features %>% mutate(overlap_Enhancer=ifelse(Enhancer,"yes","no"))
rownames(selected_features)=selected_features$gene_name
selected_features <- selected_features[,-1]
selected_features <- selected_features[,colnames(selected_features)!="Enhancer"]

# annot colors
annotation_colors = list(
  exonic_type = c("monoexonic" = "#BC3C29FF", "multiexonic"="#0072B5FF"),

 evidence_level=c("none"="white","both"="#E18727FF", "DNA"="#0072B5FF"),
 Annotation=c("potNovel"="#0072B5FF","also_Delas/Luo"="#E18727FF", "lncRNA"= "#BC3C29FF","pseudogene"="#20854EFF",
              "TEC"= "#7876B1FF"),
 classif=c("intergenic"="#20854EFF","divergent"="#0072B5FF","antisense"="#BC3C29FF","sense"="#7876B1FF","intronic_sense"="#E18727FF",
           "intronic_antisense"="#FFDC91FF","convergent"="#6F99ADFF"),
 overlap_Enhancer=c("yes"="#E18727FF","no"="white"))

#"#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF" "#FFDC91FF" "#EE4C97FF"

p=pheatmap(vst_selected_wide,
           show_rownames = F,fontsize_col =12,
           annotation_row = selected_features,scale = "none",
           annotation_colors = annotation_colors)
p
save_tiff_svg(p,outdir = outdir,filename = "heatmap_114_selected_candidates_new_sample_types_with_annot_features")

