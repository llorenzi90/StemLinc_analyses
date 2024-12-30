#
# Purpose of script: corss-link different types of evidence in lncRNAs
#
# Author: Lucia Lorenzi
#
# Date Created: 2024-12-20
#
# Email: lucialorenzi90@gmail.com
#
#  Notes ---------------------------
#
#
#
#  Setup ---------------------------

library(tidyverse)
library(DBI)
library(RSQLite)
library(UpSetR)

# connect to the database
db <- dbConnect(RSQLite::SQLite(), "outputs/dbs/StemLinc.db")

# Retrieve gene-level tables from the metadata
gene_level_tables <- dbGetQuery(
  db,
  "SELECT table_name FROM metadata WHERE level = 'gene';"
)

# View the table names
print(gene_level_tables$table_name)

# Exclude specific tables (example: exclude "metadata")
excluded_tables <- c("gene_overlap_RepeatMasker_antisense_orientation", "gene_overlap_RepeatMasker_sense_orientation")
filtered_tables <- setdiff(gene_level_tables$table_name, excluded_tables)

# Load only the filtered tables
gene_level_data <- list()
for (table in filtered_tables) {
  gene_level_data[[table]] <- dbReadTable(db, table)
}

# Start with the first table
merged_gene_table <- gene_level_data[[1]]

# Merge all subsequent tables
for (table_name in names(gene_level_data)[-1]) {
  merged_gene_table <- full_join(
    merged_gene_table,
    gene_level_data[[table_name]],
    by = "gene_name"
  )
}

# View the merged table
head(merged_gene_table)
merged_gene_table <- merged_gene_table%>% filter(pass_filter==1)
merged_gene_table <- merged_gene_table %>% filter(biotype%in%c("protein_coding",
                                                               "lncRNA","potNovel",
                                                               "TEC","pseudogene"))

library(GGally)

colnames(merged_gene_table)
continuous_features <- merged_gene_table %>% dplyr::select(max_mean_tpm_LSK,
                                                           distance_closest_PCG,
                                                           phastCons_most_conserved_exon,
                                                           phastCons_most_conserved_intron,
                                                           HSPC_vs_diff_logFC,
                                                           HSPC_vs_prog_logFC,
                                                           prog_vs_diff_logFC,
                                                           mean_corr_signature,
                                                           max_corr_signature,
                                                           max_corr,
                                                           N_TFs,
                                                           maxABCscore_asTarget,
                                                           maxABC_score_asEnh,
                                                           max_ss5_MaxEntScan,
                                                           max_ss3_MaxEntScan,
                                                           fraction_covered_by_repeats,
                                                           total_exonic_length)

continuous_features <- merged_gene_table %>% dplyr::select(max_mean_tpm_LSK,
                                                           distance_closest_PCG,
                                                           phastCons_most_conserved_exon,
                                                           mean_corr_signature,
                                                           N_TFs,
                                                           maxABCscore_asTarget,
                                                           maxABC_score_asEnh,
                                                           max_ss5_MaxEntScan,
                                                           max_ss3_MaxEntScan,
                                                           fraction_covered_by_repeats,
                                                           width)
apply(continuous_features, 2, function(x){
  print(table(is.na(x)))
  print(sum(table(is.na(x))))})

ggpairs(continuous_features)
# convert some to log2

merged_gene_table$harborsEnhancer[is.na(merged_gene_table$harborsEnhancer)] <- 0
merged_gene_table$Enhancer_target[is.na(merged_gene_table$Enhancer_target)] <- 0
# Categorical features ----
categorical_features <- merged_gene_table %>% dplyr::select(exonic_type, biotype)

table(merged_gene_table$evidence_level,merged_gene_table$harborsEnhancer,
      merged_gene_table$biotype, merged_gene_table$Enhancer_target)

table(merged_gene_table$harborsEnhancer, merged_gene_table$Enhancer_target,
      merged_gene_table$biotype)


cat_evidence <- merged_gene_table %>%
  dplyr::mutate(
    multiexonic = exonic_type == "multiexonic",
    evidence_marks = evidence_level != "none",
    harborsEnhancer = harborsEnhancer == 1,
    Enhancer_target = Enhancer_target == 1,
    high_cons0.3 = phastCons_most_conserved_exon > 0.3,
    high_Exp0.78 = max_mean_tpm_LSK > 0.78,
    low_repeat_fraction = fraction_covered_by_repeats < 0.3,
    TFs_5 = N_TFs > 4,
    high_corr_sign0.2 = mean_corr_signature > 0.2,
    HSPC_enriched = HSPC_vs_diff == "UP" | HSPC_vs_prog == "UP"
  ) %>% dplyr::select(gene_name,biotype,multiexonic,
                      evidence_marks,harborsEnhancer,
                      Enhancer_target,high_cons0.3,
                      high_Exp0.78, low_repeat_fraction,
                      TFs_5, high_corr_sign0.2, HSPC_enriched)
cat_evidence[is.na(cat_evidence)] <- FALSE

table(rowSums(cat_evidence[,3:ncol(cat_evidence)]))
rowSums(cat_evidence[,3:ncol(cat_evidence)])[cat_evidence$gene_name=="XLOC_079903"]
table(cat_evidence$biotype,rowSums(cat_evidence[,3:ncol(cat_evidence)]))

for (feat in colnames(cat_evidence)[3:ncol(cat_evidence)]) {
  print(feat)
  print(table(cat_evidence$biotype,cat_evidence[,feat]))

  }

table(cat_evidence$biotype,rowSums(cat_evidence[,4:ncol(cat_evidence)])>3,cat_evidence$multiexonic)
table(cat_evidence$biotype,rowSums(cat_evidence[,4:ncol(cat_evidence)])>4,cat_evidence$multiexonic)

table(cat_evidence$biotype,rowSums(cat_evidence[,4:ncol(cat_evidence)])>5,cat_evidence$multiexonic)
table(cat_evidence$biotype,rowSums(cat_evidence[,4:ncol(cat_evidence)])>6,cat_evidence$multiexonic)
table(cat_evidence$biotype,rowSums(cat_evidence[,4:ncol(cat_evidence)])>7,cat_evidence$multiexonic)

cat_evidence_refined <- cat_evidence %>% mutate(Enh_reg = harborsEnhancer|Enhancer_target) %>%
  dplyr::select(evidence_marks,Enh_reg,high_cons0.3,high_Exp0.78,TFs_5,high_corr_sign0.2,HSPC_enriched)

table(cat_evidence$biotype, rowSums(cat_evidence_refined))



# Convert logical columns to a binary matrix (TRUE/FALSE -> 1/0)
upset_data <- cat_evidence_refined
upset_data <- data.frame(lapply(upset_data, as.numeric))
# Create the Upset plot
for (biot in unique(cat_evidence$biotype)) {

  grid::grid.newpage()
  print(UpSetR::upset(upset_data[cat_evidence$biotype==biot,],
                nsets = ncol(upset_data),
                order.by = "freq",
                sets.bar.color = "#56B4E9",
                main.bar.color = "#0072B2",
                text.scale = 1.5,
                empty.intersections = "on"
                ) )
  # Add a title
  grid::grid.text(
    paste0("UpSet Plot of ",biot," Features"),
    x = 0.5,
    y = 0.9,
    gp = grid::gpar(fontsize = 16)
  )

}


for (biot in unique(cat_evidence$biotype)) {
  for (et in c("multiexonic", "monoexonic")) {
    grid::grid.newpage()
    print(UpSetR::upset(upset_data[merged_gene_table$biotype==biot&merged_gene_table$exonic_type==et,],
                        nsets = ncol(upset_data),
                        order.by = "freq",
                        sets.bar.color = "#56B4E9",
                        main.bar.color = "#0072B2",
                        text.scale = 1.5
    ) )
    # Add a title
    grid::grid.text(
      paste("UpSet Plot of",et,biot,"Features"),
      x = 0.5,
      y = 0.9,
      gp = grid::gpar(fontsize = 16)
    )
  }



}
# make upset plot but including more features like DNA level RNA level and enhancer target or harbor
cat_evidence <- merged_gene_table %>%
  dplyr::mutate(
    DNA = evidence_level == "DNA" | evidence_level=="both",
    RNA = evidence_level == "RNA" | evidence_level=="both",
    harborsEnhancer = harborsEnhancer == 1,
    Enhancer_target = Enhancer_target == 1,
    high_cons0.3 = phastCons_most_conserved_exon > 0.3,
    high_Exp0.78 = max_mean_tpm_LSK > 0.78,
    low_repeat_fraction0.3 = fraction_covered_by_repeats < 0.3,
    TFs_5 = N_TFs > 4,
    high_corr_sign0.2 = mean_corr_signature>0.2,
    HSPC_enriched = HSPC_vs_diff == "UP" | HSPC_vs_prog == "UP"  ) %>%
  dplyr::select(DNA, RNA,harborsEnhancer,
                      Enhancer_target,high_cons0.3,
                      high_Exp0.78, low_repeat_fraction0.3,
                      TFs_5, high_corr_sign0.2, HSPC_enriched)
cat_evidence[is.na(cat_evidence)] <- FALSE

upset_data <- cat_evidence
upset_data <- data.frame(lapply(upset_data, as.numeric))
# Create the Upset plot
for (biot in unique(merged_gene_table$biotype)) {

  grid::grid.newpage()
  print(UpSetR::upset(upset_data[merged_gene_table$biotype==biot,],
                      nsets = ncol(upset_data),
                      order.by = "freq",
                      sets.bar.color = "#56B4E9",
                      main.bar.color = "#0072B2",
                      text.scale = 1.5,
                      empty.intersections = "on"
  ) )
  # Add a title
  grid::grid.text(
    paste0("UpSet Plot of ",biot," Features"),
    x = 0.5,
    y = 0.9,
    gp = grid::gpar(fontsize = 16)
  )

}


for (biot in unique(merged_gene_table$biotype)) {
  for (et in c("multiexonic", "monoexonic")) {
    grid::grid.newpage()
    print(UpSetR::upset(upset_data[merged_gene_table$biotype==biot&merged_gene_table$exonic_type==et,],
                        nsets = ncol(upset_data),
                        order.by = "freq",
                        sets.bar.color = "#56B4E9",
                        main.bar.color = "#0072B2",
                        text.scale = 1.5
    ) )
    # Add a title
    grid::grid.text(
      paste("UpSet Plot of",et,biot,"Features"),
      x = 0.5,
      y = 0.9,
      gp = grid::gpar(fontsize = 16)
    )
  }



}
# idea: make category "overlap other genes" something like that

cat_evidence <- merged_gene_table %>%
  dplyr::mutate(
    DNA = evidence_level == "DNA" | evidence_level=="both",
    RNA = evidence_level == "RNA" | evidence_level=="both",
    harborsEnhancer = harborsEnhancer == 1,
    Enhancer_target = Enhancer_target == 1,
    overlaps_PCG = grepl("antisense|overlap|intronic",merged_gene_table$best_classif_to_PCG),
    high_cons0.3 = phastCons_most_conserved_exon > 0.3,
    high_Exp0.78 = max_mean_tpm_LSK > 0.78,
    low_repeat_fraction0.3 = fraction_covered_by_repeats < 0.3,
    TFs_5 = N_TFs > 4,
    high_corr_sign0.2 = mean_corr_signature>0.2,
    HSPC_enriched = HSPC_vs_diff == "UP" | HSPC_vs_prog == "UP"  ) %>%
  dplyr::select(DNA, RNA,harborsEnhancer,
                Enhancer_target,overlaps_PCG,high_cons0.3,
                high_Exp0.78, low_repeat_fraction0.3,
                TFs_5, high_corr_sign0.2, HSPC_enriched)
cat_evidence[is.na(cat_evidence)] <- FALSE

upset_data <- cat_evidence
upset_data <- data.frame(lapply(upset_data, as.numeric))
# Create the Upset plot
for (biot in unique(merged_gene_table$biotype)) {

  grid::grid.newpage()
  print(UpSetR::upset(upset_data[merged_gene_table$biotype==biot,],
                      nsets = ncol(upset_data),
                      order.by = "freq",
                      sets.bar.color = "#56B4E9",
                      main.bar.color = "#0072B2",
                      text.scale = 1.5
  ) )
  # Add a title
  grid::grid.text(
    paste0("UpSet Plot of ",biot," Features"),
    x = 0.5,
    y = 0.9,
    gp = grid::gpar(fontsize = 16)
  )

}


for (biot in unique(merged_gene_table$biotype)) {
  for (et in c("multiexonic", "monoexonic")) {
    grid::grid.newpage()
    print(UpSetR::upset(upset_data[merged_gene_table$biotype==biot&merged_gene_table$exonic_type==et,],
                        nsets = ncol(upset_data),
                        order.by = "freq",
                        sets.bar.color = "#56B4E9",
                        main.bar.color = "#0072B2",
                        text.scale = 1.5
    ) )
    # Add a title
    grid::grid.text(
      paste("UpSet Plot of",et,biot,"Features"),
      x = 0.5,
      y = 0.9,
      gp = grid::gpar(fontsize = 16)
    )
  }
}

# make upset for chromatin marks alone ----
chromatin_marks <- merged_gene_table %>% dplyr::select(CAGE_within_100bp,polyAsite_within_100bp,Enhancer,H3K27ac_at_promoter,H3K36me3_at_geneBody,H3K4me1_at_promoter,H3K4me3_at_promoter,H3K27me3_at_promoter)

chromatin_marks[is.na(chromatin_marks)] <- FALSE

upset_data <- chromatin_marks
upset_data <- data.frame(lapply(upset_data, as.numeric))
# Create the Upset plot
for (biot in unique(merged_gene_table$biotype)) {
  tmp_data <- upset_data[merged_gene_table$biotype==biot,]
  grid::grid.newpage()
  print(UpSetR::upset(tmp_data,
                      nsets = ncol(upset_data),
                      order.by = "freq",
                      sets.bar.color = "#56B4E9",
                      main.bar.color = "#0072B2",
                      text.scale = 1.5,
                      empty.intersections = "on"
  ) )
  # Add a title
  grid::grid.text(
    paste0("UpSet Plot of ",biot," Features"),
    x = 0.5,
    y = 0.9,
    gp = grid::gpar(fontsize = 16)
  )

  # Calculate pairwise overlaps
  pairwise_overlap <- t(as.matrix(tmp_data)) %*% as.matrix(tmp_data)

  # Convert to data frame for inspection (optional)
  pairwise_overlap_df <- as.data.frame(as.table(pairwise_overlap))

  # View the pairwise overlaps
  print(pairwise_overlap_df)
  colnames(pairwise_overlap_df) <- c("Feature1", "Feature2", "OverlapCount")

  # Create a heatmap for pairwise overlaps
  g=ggplot(pairwise_overlap_df, aes(x = Feature1, y = Feature2, fill = OverlapCount)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    theme_minimal() +
    ggtitle("Pairwise Overlap Between Features")
  print(g)
}




for (biot in unique(merged_gene_table$biotype)) {
  for (et in c("multiexonic", "monoexonic")) {
    grid::grid.newpage()
    print(UpSetR::upset(upset_data[merged_gene_table$biotype==biot&merged_gene_table$exonic_type==et,],
                        nsets = ncol(upset_data),
                        order.by = "freq",
                        sets.bar.color = "#56B4E9",
                        main.bar.color = "#0072B2",
                        text.scale = 1.5
    ) )
    # Add a title
    grid::grid.text(
      paste("UpSet Plot of",et,biot,"Features"),
      x = 0.5,
      y = 0.9,
      gp = grid::gpar(fontsize = 16)
    )
  }
}
# upset plot with different enhancer definitions ----

enhancer_defs <- merged_gene_table %>% dplyr::select(Enhancer_Atlas_at_promoter,FANTOM_enhancers_at_promoter,
                                                     ABC_Enhancer = harborsEnhancer,Enhancer_from_marks_at_promoter)


enhancer_defs[is.na(enhancer_defs)] <- FALSE

upset_data <- enhancer_defs
upset_data <- data.frame(lapply(upset_data, as.numeric))
# Create the Upset plot
for (biot in unique(merged_gene_table$biotype)) {
  tmp_data <- upset_data[merged_gene_table$biotype==biot,]
  grid::grid.newpage()
  print(UpSetR::upset(tmp_data,
                      nsets = ncol(upset_data),
                      order.by = "freq",
                      sets.bar.color = "#56B4E9",
                      main.bar.color = "#0072B2",
                      text.scale = 1.5,
                      empty.intersections = "on"
  ) )
  # Add a title
  grid::grid.text(
    paste0("UpSet Plot of ",biot," Features"),
    x = 0.5,
    y = 0.9,
    gp = grid::gpar(fontsize = 16)
  )

  # Calculate pairwise overlaps
  pairwise_overlap <- t(as.matrix(tmp_data)) %*% as.matrix(tmp_data)

  # Convert to data frame for inspection (optional)
  pairwise_overlap_df <- as.data.frame(as.table(pairwise_overlap))

  # View the pairwise overlaps
  print(pairwise_overlap_df)
  colnames(pairwise_overlap_df) <- c("Feature1", "Feature2", "OverlapCount")

  # Create a heatmap for pairwise overlaps
  g=ggplot(pairwise_overlap_df, aes(x = Feature1, y = Feature2, fill = OverlapCount)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    theme_minimal() +
    ggtitle("Pairwise Overlap Between Features")
  print(g)
}




for (biot in unique(merged_gene_table$biotype)) {
  for (et in c("multiexonic", "monoexonic")) {
    grid::grid.newpage()
    print(UpSetR::upset(upset_data[merged_gene_table$biotype==biot&merged_gene_table$exonic_type==et,],
                        nsets = ncol(upset_data),
                        order.by = "freq",
                        sets.bar.color = "#56B4E9",
                        main.bar.color = "#0072B2",
                        text.scale = 1.5
    ) )
    # Add a title
    grid::grid.text(
      paste("UpSet Plot of",et,biot,"Features"),
      x = 0.5,
      y = 0.9,
      gp = grid::gpar(fontsize = 16)
    )
  }
}
# make a plot including "outliers" for some features, e.g,
#  genes with very high values compared to the distribution


# define a simple score:

basic_evidence_data <- upset_data %>% dplyr::select(DNA,RNA,harborsEnhancer,Enhancer_target,TFs_5)
table(merged_gene_table$biotype, rowSums(basic_evidence_data)>2,merged_gene_table$exonic_type)

basic_evidence_data <- basic_evidence_data %>% mutate(DNA_RNA=DNA|RNA,
                                                      Enh_reg=Enhancer_target|harborsEnhancer,
                                                      TFs_5) %>% select(DNA_RNA,Enh_reg,TFs_5)
table(merged_gene_table$biotype,
      rowSums(basic_evidence_data),
      merged_gene_table$exonic_type)

# define a more complex scoring system: think about it but not too much to start with
# We can always refine it later...
# I want to have a set of lncRNAs of which I am confident they are:
# 1. expressed consistently: for this, I've required reproducibility in our samples,
# I might check expression across blood cells: if not reproducible
# then their expression should be high enough (a reason why they are not found in other
# samples could be that they are non-polyadenylated), they have CAGE signals
# 2. likely regulated (by presence
# of TFs at their promoters, they have chromatin marks ar promoters, they are targets of enhancers (select for those with high ABC scores and
# in two replicates??), they have ATAC-seq signal(to be added),
# they have predicted binding to RNA, RBPs
# 3. likely regulators: predicted binding to miRNAs, close to important genes (TAD analysis),
# expression correlation with signature genes
# 4. transcribed from strong enhancers, prioritize genes that regulate HSC-related genes,
# where are the genes that these potNovel enhancers regulate? in their TADs?
# 5. Develop a simple scoring system with all data I have so far.
# 6. Try to classify genes into subclasses
# 7. Be careful with genes that overlap other genes, especially intronic ones, maybe focus on
# intronic antisense
