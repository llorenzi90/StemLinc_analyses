gene_level_path <- "outputs/transcriptome_characterization/jan25/LSK_StemLinc.combined_annotated_tracking.filtered.sel_biots.20250127_112131.gene_level.tsv"
source("scripts/functions/nejm_palette.R")
library(DBI)
library(RSQLite)
library(dbplyr)
library(tidyverse)
library(rtracklayer)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
require(ggvenn)
library(DESeq2)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)


db <- dbConnect(RSQLite::SQLite(), "/home/llorenzi/Rprojects/StemLinc_analyses/outputs/dbs/StemLinc.db")

dbListTables(db)

gene_core_data <- dbReadTable(db,"gene_core_data")

gene_level_data <- read.delim(gene_level_path)

table(gene_level_data$gene_name%in%gene_core_data$gene_name)
gene_level_data[!gene_level_data$gene_name%in%gene_core_data$gene_name,]

# Check characteristics of monoexons ----
# I know monoexons have in general less "evidence" that other lncRNAs
# But in all features I plot I often see there are outliers
# Are outliers in one feature enriched for other features?
# Should I use a combination of evidence features to retain only confident monoexons
# Compare to expressed monoexonic annotated lncRNAs in our assembly:
#   Are there annotated multiexonic isoforms in the same locus that are not retrieved?
#   How do evidence look like for these?

# Genomic classification ----
gene_classif <- dbReadTable(db,"gene_classification_to_PCG")

gene_level_data <- left_join(gene_level_data,gene_classif%>%dplyr::select(gene_name,best_classif_to_PCG,distance_closest_PCG))

table(gene_level_data$best_classif_to_PCG,gene_level_data$biotype,gene_level_data$exonic_type)

ggplot(gene_level_data %>% filter(biotype!="protein_coding"), aes(x=biotype,fill=best_classif_to_PCG)) + geom_bar() +
  facet_wrap(~exonic_type) + theme_minimal() +
  scale_fill_manual(values = colorRampPalette(nejm_pal)(length(unique(gene_level_data$best_classif_to_PCG))))

# Length ----
ggplot(gene_level_data, aes(x=biotype, y=width)) + geom_boxplot() +
  facet_wrap(~exonic_type) + theme_minimal() +
  scale_y_log10()

gene_level_data %>% group_by(biotype,exonic_type) %>%
  reframe(mean(width),
            median(width))
# Expression ----

ggplot(gene_level_data , aes(x=max_mean_tpm,col=biotype)) +
  geom_density(linewidth = 3) + theme_minimal() + facet_wrap(~exonic_type)+
  scale_x_log10()

ggplot(gene_level_data , aes(x=biotype,y=max_mean_tpm)) +
  geom_boxplot() + theme_minimal() + facet_wrap(~exonic_type)+
  scale_y_log10()

gene_level_data %>% group_by(biotype,exonic_type) %>%
  reframe(mean(max_mean_tpm),
          median(max_mean_tpm))

gene_level_data %>% group_by(biotype,exonic_type) %>%
  reframe(min(max_mean_tpm),
          max(max_mean_tpm))

# what if I filter by 0.5 TPM:
gene_level_data_0.5TPM <- gene_level_data %>% filter(max_mean_tpm>=0.5)
table(gene_level_data_0.5TPM$biotype,
      gene_level_data_0.5TPM$exonic_type)
ggplot(gene_level_data_0.5TPM , aes(x=biotype,y=max_mean_tpm)) +
  geom_boxplot() + theme_minimal() + facet_wrap(~exonic_type)+
  scale_y_log10() + ggtitle("Gene expression filtered at 0.5 TPM")

# what if I filter by 1 TPM:
gene_level_data_1TPM <- gene_level_data %>% filter(max_mean_tpm>=1)
table(gene_level_data_1TPM$biotype,
      gene_level_data_1TPM$exonic_type)
ggplot(gene_level_data_1TPM , aes(x=biotype,y=max_mean_tpm)) +
  geom_boxplot() + theme_minimal() + facet_wrap(~exonic_type)+
  scale_y_log10() + ggtitle("Gene expression filtered at 1 TPM")
# Conservation ----
gene_conservation_phastCons35way <- dbReadTable(db,"gene_conservation_phastCons35way")

gene_level_data <- left_join(gene_level_data,
                             gene_conservation_phastCons35way%>%
                               dplyr::select(gene_name,phastCons_mean_merged_exons,phastCons_most_conserved_promoter,phastCons_mean_introns))


ggplot(gene_level_data , aes(x=biotype,y=phastCons_mean_merged_exons)) +
  geom_boxplot() + theme_minimal() + facet_wrap(~exonic_type) +
  ggtitle("Mean phastCons exonic regions")

ggplot(gene_level_data , aes(x=biotype,y=phastCons_most_conserved_promoter)) +
  geom_boxplot() + theme_minimal() + facet_wrap(~exonic_type) +
  ggtitle("Mean phastCons most conserved promoter region")

gene_level_data$phastCons_promoter_vs_merged_Exons <-
  gene_level_data$phastCons_most_conserved_promoter/gene_level_data$phastCons_mean_merged_exons

ggplot(gene_level_data , aes(x=biotype,y=phastCons_promoter_vs_merged_Exons)) +
  geom_boxplot() + theme_minimal() + facet_wrap(~exonic_type) +
  ggtitle("Mean phastCons promoter/exons") +
  scale_y_log10()

# summary exon conservation
gene_level_data %>% filter(!is.na(phastCons_mean_merged_exons)) %>%
  group_by(biotype, exonic_type) %>% summarise(mean=mean(phastCons_mean_merged_exons),
                                               fraction_higher_promoter=sum(phastCons_mean_merged_exons>mean)/n(),
                                               median=median(phastCons_mean_merged_exons))

#summary ratio promoter/exon conservation

gene_level_data %>% filter(!is.na(phastCons_promoter_vs_merged_Exons)&is.finite(phastCons_promoter_vs_merged_Exons)) %>%
  group_by(biotype, exonic_type) %>% summarise(mean=mean(phastCons_promoter_vs_merged_Exons),
                                               fraction_higher_promoter=sum(phastCons_promoter_vs_merged_Exons>1)/n(),
                                               median=median(phastCons_promoter_vs_merged_Exons))


# Plot both together ----
phastcons_data <- gene_level_data %>% dplyr::select(gene_name,biotype,exonic_type,
                                             phastCons_mean_merged_exons,phastCons_most_conserved_promoter)%>%
  pivot_longer(cols = 4:5,names_to = "level",values_to = "phastCons")

ggplot(phastcons_data,aes(x=biotype, y=phastCons,fill=level)) +
  geom_boxplot() + facet_wrap(~exonic_type) +theme_minimal() +
  scale_fill_manual(values=nejm_pal)

# correlation expression - conservation

ggplot(gene_level_data, aes(x=max_mean_tpm,y=phastCons_most_conserved_promoter))+
  geom_point() + theme_minimal() + scale_x_log10() + geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add linear trend line
  stat_cor(method = "pearson", position = "jitter")
cor(gene_level_data$max_mean_tpm[!is.na(gene_level_data$max_mean_tpm)&!is.na(gene_level_data$phastCons_most_conserved_promoter)],
    gene_level_data$phastCons_most_conserved_promoter[!is.na(gene_level_data$max_mean_tpm)&!is.na(gene_level_data$phastCons_most_conserved_promoter)])


#[1] 0.008468279
summary(gene_level_data$phastCons_most_conserved_promoter)


# conservation density distribution

# exonic level
ggplot(gene_level_data, aes(x=phastCons_mean_merged_exons,col=biotype)) + geom_density(linewidth=2) +
  theme_minimal() + facet_wrap(~exonic_type)

# promoter level
ggplot(gene_level_data, aes(x=phastCons_most_conserved_promoter,col=biotype)) + geom_density(linewidth=2) +
  theme_minimal() + facet_wrap(~exonic_type)

# summary different cutoffs ----
for (co in c(0.1,0.2,0.3)) {
  gl_s <- gene_level_data %>% group_by(biotype,exonic_type) %>%
    summarise(sum(phastCons_mean_merged_exons>co,na.rm = T),
              sum(phastCons_mean_merged_exons>co,na.rm = T)/n())
  print(co)
  print(gl_s)
  gl_s <- gene_level_data %>% group_by(biotype,exonic_type) %>%
    summarise(sum(phastCons_most_conserved_promoter>co,na.rm = T),
              sum(phastCons_most_conserved_promoter>co,na.rm = T)/n())
  print("promoter")
  print(gl_s)

}

# correlation conservation promoter exon
ggplot(gene_level_data, aes(x=phastCons_mean_merged_exons,y=phastCons_most_conserved_promoter))+
  geom_point() + theme_minimal() + geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add linear trend line
  stat_cor(method = "pearson", position = "jitter")

# How is the expression of genes with high conservation?
co=0.2
gene_level_data_high_cons <- gene_level_data %>% filter(!is.na(phastCons_mean_merged_exons)&
                                                          phastCons_mean_merged_exons>=co)

ggplot(gene_level_data_high_cons,aes(x=biotype,y=max_mean_tpm)) +geom_boxplot() +
  theme_minimal() +facet_wrap(~exonic_type) + scale_y_log10() + ggtitle("High conservation genes (phastCons exons >=0.2)")

gene_level_data <- gene_level_data %>% mutate(highcons=!is.na(phastCons_mean_merged_exons)&phastCons_mean_merged_exons>=co)
ggplot(gene_level_data%>%filter(!is.na(phastCons_mean_merged_exons)),aes(x=biotype, fill=highcons,
                                                                                       y=max_mean_tpm)) +geom_boxplot() +
  theme_minimal() +facet_wrap(~exonic_type) + scale_y_log10() + ggtitle("High conservation genes (phastCons exons >=0.2)")

gene_level_data <- gene_level_data %>% mutate(highcons_promoter=!is.na(phastCons_most_conserved_promoter)&
                                                phastCons_most_conserved_promoter>=co)
ggplot(gene_level_data%>%filter(!is.na(phastCons_most_conserved_promoter)),aes(x=biotype, fill=highcons_promoter,
                                                                         y=max_mean_tpm)) +geom_boxplot() +
  theme_minimal() +facet_wrap(~exonic_type) + scale_y_log10() + ggtitle("High conservation genes (phastCons promoter >=0.2)")

# CAGE ----
cmarks <- dbReadTable(db,"gene_overlap_with_chromatin_marks")
gene_level_data <- left_join(gene_level_data,
                             cmarks %>%
                               dplyr::select(gene_name,evidence_level,
                                             Enhancer,CAGE_within_100bp,polyAsite_within_100bp,
                                            H3K4me3_at_promoter,H3K27ac_at_promoter,
                                            H3K4me1_at_promoter,H3K36me3_at_geneBody,H3K27me3_at_promoter))


ggplot(gene_level_data %>% filter(!is.na(CAGE_within_100bp)),
       aes(x=biotype,fill=CAGE_within_100bp==1)) + geom_bar(position = "fill")+
  theme_minimal() +facet_wrap(~exonic_type) + scale_fill_manual(values = nejm_pal)

gene_level_data_0.5TPM <- gene_level_data %>% filter(max_mean_tpm>=0.5)
ggplot(gene_level_data_0.5TPM %>% filter(!is.na(CAGE_within_100bp)),
       aes(x=biotype,fill=CAGE_within_100bp==1)) + geom_bar(position = "fill")+
  theme_minimal() +facet_wrap(~exonic_type) + scale_fill_manual(values = nejm_pal)+
  ggtitle(">=0.5TPM")

table(gene_level_data$biotype,
      gene_level_data$exonic_type,
      gene_level_data$CAGE_within_100bp==1)
table(gene_level_data_0.5TPM$biotype,
      gene_level_data_0.5TPM$exonic_type,
      gene_level_data_0.5TPM$CAGE_within_100bp==1)
  # expression with and without cage
ggplot(gene_level_data %>% filter(!is.na(CAGE_within_100bp)),
       aes(x=biotype,fill=CAGE_within_100bp==1,y=max_mean_tpm)) + geom_boxplot()+
  theme_minimal() +facet_wrap(~exonic_type) + scale_fill_manual(values = nejm_pal) +
  scale_y_log10()
# Chromatin marks ----
colnames(gene_level_data)

for (cn in colnames(gene_level_data)[21:27]) {
  p <- ggplot(gene_level_data %>% filter(!is.na(.data[[cn]])),
              aes(x = biotype, fill = .data[[cn]] == 1)) +
    geom_bar(position = "fill") +
    theme_minimal() +
    facet_wrap(~exonic_type) +
    scale_fill_manual(values = nejm_pal) +
    ggtitle(paste0("Fraction of genes with ",cn))

  print(p)
}

# Enhancer form marks, Enhancer Atlas and FANTOM ----
gene_level_data <- left_join(gene_level_data, cmarks %>% dplyr::select(gene_name,Enhancer_from_marks_at_promoter,Enhancer_Atlas_at_promoter,
                                                                FANTOM_enhancers_at_promoter,any_Enhancer=Enhancer))

colnames(gene_level_data)
for (cn in colnames(gene_level_data)[28:31]) {
  p <- ggplot(gene_level_data %>% filter(!is.na(.data[[cn]])),
              aes(x = biotype, fill = .data[[cn]] == 1)) +
    geom_bar(position = "fill") +
    theme_minimal() +
    facet_wrap(~exonic_type) +
    scale_fill_manual(values = nejm_pal) +
    ggtitle(paste0("Fraction of genes with ",cn))

  print(p)
}

# Venn diagram overlap enhancers ----
plot_three_venn <- function(traGL,
                            samp1_var,
                            samp2_var,
                            samp3_var,
                            feat_var,
                            filter_cond=rep(T,nrow(traGL)), # defaults the entire data
                            sampnames=c(sample_names,"Ref"),
                            title=""){
  traGL <- traGL %>% filter(filter_cond)
  if(is.logical(traGL%>%pull({{samp1_var}}))){
    vennlist <- list(traGL %>% filter({{samp1_var}}) %>% pull({{feat_var}}) ,
                     traGL %>% filter({{samp2_var}}) %>% pull({{feat_var}}),
                     traGL %>% filter({{samp3_var}}) %>% pull({{feat_var}}))
  }else{
    vennlist <- list(traGL %>% filter({{samp1_var}}!="-") %>% pull({{feat_var}}) ,
                     traGL %>% filter({{samp2_var}}!="-") %>% pull({{feat_var}}),
                     traGL %>% filter({{samp3_var}}!="-") %>% pull({{feat_var}}))
  }

  names(vennlist)=sampnames

  g=ggvenn(
    vennlist,
    fill_color = c("#0072B5FF","#FFDC91FF","#20854EFF"),
    stroke_size = 0.5, set_name_size = 4
  ) +ggtitle(title)
  print(g)
  #return(list(vennlist,g))
} # makes venn diagram of three sets of transcripts

enh_data <- gene_level_data %>% dplyr::select(biotype, gene_name,Enhancer_Atlas_at_promoter,Enhancer_from_marks_at_promoter,
                                       FANTOM_enhancers_at_promoter)

enh_data <- enh_data %>% drop_na()

enh_data[,3:ncol(enh_data)] <- enh_data[,3:ncol(enh_data)] == 1

for (bio in unique(enh_data$biotype)) {
  plot_three_venn(traGL = enh_data,samp1_var = Enhancer_Atlas_at_promoter,
                  samp2_var = Enhancer_from_marks_at_promoter,samp3_var = FANTOM_enhancers_at_promoter,
                  filter_cond = enh_data$biotype%in%bio,
                  feat_var = gene_name, sampnames = colnames(enh_data)[3:5],
                  title = bio)

}
# Enhancer from ABC ----
ABC_predictions <- dbReadTable(db,"gene_ABC_predictions_enhancers_and_or_targets")
ABC_predictions$biotype <- gene_level_data$biotype[match(ABC_predictions$gene_name,
                                                         gene_level_data$gene_name)]

table(ABC_predictions$biotype,ABC_predictions$harborsEnhancer)

ABC_harbors_enhancer <- ABC_predictions$gene_name[ABC_predictions$harborsEnhancer==1]

gene_level_data$Enhancer_from_ABC <- gene_level_data$gene_name%in%ABC_harbors_enhancer

table(gene_level_data$Enhancer_from_ABC,gene_level_data$biotype,gene_level_data$exonic_type)

cn="Enhancer_from_ABC"
p <- ggplot(gene_level_data %>% filter(!is.na(.data[[cn]])),
            aes(x = biotype, fill = .data[[cn]] == 1)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  facet_wrap(~exonic_type) +
  scale_fill_manual(values = nejm_pal) +
  ggtitle(paste0("Fraction of genes with ",cn))

print(p)

enh_data <- left_join(enh_data,gene_level_data%>%dplyr::select(gene_name,Enhancer_from_ABC))
for (bio in unique(enh_data$biotype)) {
  plot_three_venn(traGL = enh_data,samp1_var = Enhancer_Atlas_at_promoter,
                  samp2_var = Enhancer_from_ABC,samp3_var = FANTOM_enhancers_at_promoter,
                  filter_cond = enh_data$biotype%in%bio,
                  feat_var = gene_name, sampnames = colnames(enh_data)[c(3,6,5)],
                  title = bio)

}

# Enhancer in all 4 ----
gene_level_data <- gene_level_data %>% mutate(Enhancer_all4=Enhancer_from_marks_at_promoter==1&
                                                Enhancer_from_ABC&Enhancer_Atlas_at_promoter==1&FANTOM_enhancers_at_promoter==1)
table(gene_level_data$Enhancer_all4,
      gene_level_data$biotype,
      gene_level_data$exonic_type)

table(gene_level_data$FANTOM_enhancers_at_promoter,
      gene_level_data$biotype,
      gene_level_data$exonic_type)

table(gene_level_data$Enhancer_Atlas_at_promoter,
      gene_level_data$biotype,
      gene_level_data$exonic_type)

table(gene_level_data$Enhancer_from_marks_at_promoter,
      gene_level_data$biotype,
      gene_level_data$exonic_type)

table(gene_level_data$Enhancer_from_ABC,
      gene_level_data$biotype,
      gene_level_data$exonic_type)

# how many have CAGE and polyA
table(gene_level_data$CAGE_within_100bp==1&gene_level_data$polyAsite_within_100bp==1,
      gene_level_data$biotype,
      gene_level_data$exonic_type)
# how many have CAGE and H3K27ac?
table(gene_level_data$CAGE_within_100bp==1&gene_level_data$H3K27ac_at_promoter==1,
      gene_level_data$biotype,
      gene_level_data$exonic_type)


# how many have any mark?
table(gene_level_data$evidence_level,
      gene_level_data$biotype,
      gene_level_data$exonic_type)

table(gene_level_data$evidence_level,
      gene_level_data$Enhancer_all4,
      gene_level_data$biotype,
      gene_level_data$exonic_type)

# evidence level ----
gene_level_data <- gene_level_data %>% mutate(evidence_level=factor(evidence_level,levels=c("none","both","RNA","DNA")))
ggplot(gene_level_data%>%filter(!is.na(max_mean_tpm)&!is.na(evidence_level)),aes(x=biotype,fill=evidence_level)) +
  geom_bar(color="black") + theme_minimal() +
  scale_fill_manual(values = c("white",nejm_pal[c(3,2,4)])) + theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggplot(gene_level_data%>%filter(!is.na(max_mean_tpm)&!is.na(evidence_level)),aes(x=biotype,fill=evidence_level)) +
  geom_bar(color="black",position = "fill") + theme_minimal() +
  scale_fill_manual(values = c("white",nejm_pal[c(3,2,4)])) + theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggplot(gene_level_data%>%filter(!is.na(max_mean_tpm)&!is.na(evidence_level)),aes(x=biotype,fill=evidence_level)) +
  geom_bar(color="black",position = "fill") + theme_minimal() + facet_wrap(~ exonic_type) +
  scale_fill_manual(values = c("white",nejm_pal[c(3,2,4)])) + theme(axis.text.x = element_text(angle = 45,hjust = 1))

# expression per evidence level ----

ggplot(gene_level_data%>%filter(!is.na(max_mean_tpm)&!is.na(evidence_level)),aes(x=biotype,fill=evidence_level,y=max_mean_tpm)) +
  geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("white",nejm_pal[c(3,2,4)]))+
  scale_y_log10()

ggplot(gene_level_data%>%filter(!is.na(max_mean_tpm)&!is.na(evidence_level)),aes(x=biotype,fill=evidence_level,y=max_mean_tpm)) +
  geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("white",nejm_pal[c(3,2,4)]))+
  scale_y_log10() + facet_wrap(~exonic_type)

# Summary of evidence level ----
gene_level_data %>% filter(!is.na(evidence_level)) %>% group_by(biotype,exonic_type) %>%
  summarise(N_any_evidence=sum(evidence_level!="none"),
            fraction_any_evidence=sum(evidence_level!="none")/n())


# Correlation ----
gene_correlation_with_PCG_signature <- dbReadTable(db,"gene_correlation_with_PCG_signature")
gene_level_data <- left_join(gene_level_data,gene_correlation_with_PCG_signature)
summary(gene_correlation_with_PCG_signature$max_corr[gene_correlation_with_PCG_signature$gene_name%in%gene_level_data$gene_name[gene_level_data$biotype=="potNovel"]])
summary(gene_correlation_with_PCG_signature$max_corr)
table(gene_level_data$mean_corr_signature>0.5,
      gene_level_data$biotype,
      gene_level_data$exonic_type)

ggplot(gene_level_data, aes(x=biotype,y=max_corr_signature)) + geom_boxplot() +
  facet_wrap(~exonic_type)+ theme_minimal() + ggtitle("max correlation with HSC signature")
table(gene_level_data$mean_corr_signature>0.6,
      gene_level_data$biotype,
      gene_level_data$exonic_type)

# TF binding ----
gene_overlap_TFs_cistromes_rellevA <- dbReadTable(db,"gene_overlap_TFs_cistromes_rellevA")
gene_overlap_TFs_cistromes_rellevA$N_TFs
gene_level_data <- left_join(gene_level_data,gene_overlap_TFs_cistromes_rellevA)

ggplot(gene_level_data, aes(x=biotype,y=N_TFs)) + geom_boxplot() +
  facet_wrap(~exonic_type)+ theme_minimal() + ggtitle("TF binding at promoters")
table(gene_level_data$N_TFs>5,
      gene_level_data$biotype,
      gene_level_data$exonic_type)

# Targets of enhancer ----
gene_level_data <- mutate(gene_level_data,enhancer_target=gene_name%in%ABC_predictions$gene_name[ABC_predictions$Enhancer_target==1])
table(gene_level_data$enhancer_target,
      gene_level_data$biotype,
      gene_level_data$exonic_type)


# Develop an "scoring" system ----
# keep a most stringent set and a more inclusive set.
### high priority features (require 2): ----
# presence of CAGE
# presence of chromatin mark... which?? H3K4me3
# very high expression >=2 TMP
# very high conservation >= 0.5
# very high number of TFs >=40
# very high correlation with signature >=0.6

# rationale: CAGE is a strong independent indicator of transcription
# H3K4me3 is also a strong indicator of active transcription,
# prioritize this before H3K27ac bc the last is also indicator of enhancer
# and we do not want to include enhancers per-se unless they have other evidence
# of being important
# then the "very high" features are here to rescue "outliers" for
# different features of potential interest
#
# possible missings: genes that do not have a cap (how frequent does it happen?)
#
# should I include the other chromatin marks?
#
# fears: this semi-scoring system does not necessarily keep functional genes
# should I further filter by expression?? for differential analysis
# ideas: keep it less strict and look for functional evidence:
# another option is to try different filters for different analyses
# another option is to sort by score
# differential expression, should I filter before this??
# how many of the low priority are differentially expressed?
# check how many of the high corr are also differentially expressed
#
# check how many known lncRNAs enriched in HSC pass the filters
#

high_prior_feats <- gene_level_data %>% transmute(gene_name, biotype, exonic_type,
                                                  CAGE=CAGE_within_100bp==1,
                                                  H3K4me3_at_promoter=H3K4me3_at_promoter==1,

                                                  highExpression=max_mean_tpm>=2,
                                                  highConservation=phastCons_mean_merged_exons>=0.5,
                                                  highTFbinding=N_TFs>=40,
                                                  highCorrSign=max_corr_signature>=0.6)

table(high_prior_feats$biotype,
      high_prior_feats$exonic_type,
      rowSums(high_prior_feats[,4:ncol(high_prior_feats)])>1)

gene_level_data$highPrior_feats=rowSums(high_prior_feats[,4:ncol(high_prior_feats)])


### midium prority features (require 2): ----
# H3K27ac #maybe move above
# H3K36me3_at_geneBody
# 3 enhancers
# high expression >=1 TPM

N_enhancers <- rowSums(enh_data[,3:6])
gene_level_data$N_enh_evidence <- N_enhancers[match(gene_level_data$gene_name,
                                                    enh_data$gene_name)]

table(gene_level_data$N_enh_evidence,
      gene_level_data$biotype,
      gene_level_data$exonic_type)


mid_prior_feats <- gene_level_data %>% transmute(gene_name, biotype, exonic_type,
                                                 H3K27ac_at_promoter=H3K27ac_at_promoter==1,
                                                 H3K36me3_at_geneBody=H3K36me3_at_geneBody==1,
                                                  midhighExpression=max_mean_tpm>=1,
                                                  N_enh_evidence>=3)


table(mid_prior_feats$biotype,
      mid_prior_feats$exonic_type,
      rowSums(mid_prior_feats[,4:ncol(mid_prior_feats)])>1)

gene_level_data$midPrior_feats=rowSums(mid_prior_feats[,4:ncol(mid_prior_feats)])

table(gene_level_data$highPrior_feats,gene_level_data$midPrior_feats,
      gene_level_data$biotype)

gene_level_data <- gene_level_data %>% mutate(Tier1=highPrior_feats>1,
                                              Tier2=midPrior_feats>1|highPrior_feats>0)

table(gene_level_data$Tier1,gene_level_data$biotype)
table(gene_level_data$biotype, gene_level_data$exonic_type,gene_level_data$Tier1)
table(gene_level_data$Tier2,gene_level_data$biotype)
table(gene_level_data$biotype, gene_level_data$exonic_type,gene_level_data$Tier2)


min(gene_level_data$width[gene_level_data$biotype=="potNovel"&gene_level_data$Tier1])
min(gene_level_data$max_mean_tpm[gene_level_data$biotype=="potNovel"&gene_level_data$Tier1])
View(gene_level_data %>% arrange(desc(highPrior_feats),desc(max_mean_tpm),desc(midPrior_feats)) %>% filter(biotype=="potNovel"))
View(gene_level_data %>% arrange(desc(highPrior_feats),desc(midPrior_feats),desc(max_mean_tpm)) %>% filter(biotype=="potNovel"))

# I like the idea of sorting first by high priority features,
# then by expression and then by mid priority?


# More refined scoring system
# 1. Define Features and Assign Weights
# Assign weights to each feature based on its biological relevance and strength of evidence for functionality. Here's an example of how you might weight the features:
#
# High-Priority Features (Strong Evidence)
# CAGE peak within 100 bp of TSS: 3 points
#
# Expression (TPM): Min-max normalized × 3
#
# Correlation with signature genes (≥ 0.6): 2 points
#
# Conservation at exonic regions: Min-max normalized × 2
#
# TF binding at promoters: Min-max normalized × 2
#
# Medium-Priority Features
# Chromatin marks (H3K4me3, H3K27ac, H3K36me3): 1 point each
#
# Enhancer overlaps (Enhancer_Atlas, FANTOM): 1 point each
#
# Enhancer_from_marks: 0.5 points
#
# Conservation at promoter regions: Min-max normalized × 1
#
# Low-Priority Features
# Length > 500 bp: 0.5 points
#
# Min-max normalization function
min_max_normalize <- function(x) {
  (x - min(x,na.rm = T)) / (max(x,na.rm = T) - min(x,na.rm = T))
}

# Apply min-max normalization to continuous features
potNovel_genes <- gene_level_data%>%filter(biotype=="potNovel")

potNovel_genes <- potNovel_genes %>%
  mutate(
    Expression_Norm = min_max_normalize(log10(max_mean_tpm+1)),
    Exonic_Conservation_Norm = min_max_normalize(phastCons_mean_merged_exons),
    Promoter_Conservation_Norm = min_max_normalize(phastCons_most_conserved_promoter),
    TF_Binding_Norm = min_max_normalize(N_TFs)
  )

# Calculate total score
potNovel_genes <- potNovel_genes %>%
  mutate(Total_Score= CAGE_within_100bp*3 + Expression_Norm*3 +
         ifelse(max_corr_signature >= 0.6,2,0) +
          Exonic_Conservation_Norm * 2 +
           TF_Binding_Norm * 2 +
           H3K4me3_at_promoter * 1 + H3K27ac_at_promoter * 1 +
           H3K36me3_at_geneBody * 1 +
           Enhancer_Atlas_at_promoter * 1 + FANTOM_enhancers_at_promoter  * 1 +
           Enhancer_from_marks_at_promoter * 0.5 + Enhancer_from_ABC *1 +
           Promoter_Conservation_Norm * 1 +
           (ifelse(width > 500, 0.5, 0)) )


View(potNovel_genes %>% dplyr::select(gene_name,exonic_type,CAGE_within_100bp,Expression_Norm,
                               max_corr_signature,Exonic_Conservation_Norm,TF_Binding_Norm,
                               H3K4me3_at_promoter,H3K27ac_at_promoter,H3K36me3_at_geneBody,
                               Enhancer_Atlas_at_promoter,
                               FANTOM_enhancers_at_promoter,
                               Enhancer_from_ABC,Promoter_Conservation_Norm,
                               width,chr,Total_Score))

potNovel_genes <- potNovel_genes %>%
  mutate(Total_Score= CAGE_within_100bp*3 + Expression_Norm*3 +
           ifelse(max_corr_signature >= 0.6,2,0) +

           TF_Binding_Norm * 1 +
           H3K4me3_at_promoter * 1 + H3K27ac_at_promoter * 1 +
           H3K36me3_at_geneBody * 0.5 +
           Enhancer_Atlas_at_promoter * 0.5 + FANTOM_enhancers_at_promoter  * 0.5 +
           Enhancer_from_marks_at_promoter * 0.5 + Enhancer_from_ABC *1 +
           Promoter_Conservation_Norm * 1 +
           Exonic_Conservation_Norm * 2 +
           (ifelse(width > 500, 0.5, 0)) )

View(potNovel_genes %>% dplyr::select(gene_name,exonic_type,CAGE_within_100bp,Expression_Norm,
                               max_corr_signature,Exonic_Conservation_Norm,TF_Binding_Norm,
                               H3K4me3_at_promoter,H3K27ac_at_promoter,H3K36me3_at_geneBody,
                               Enhancer_Atlas_at_promoter,
                               FANTOM_enhancers_at_promoter,
                               Enhancer_from_ABC,Promoter_Conservation_Norm,
                               width,chr,Total_Score) %>% arrange(desc(Total_Score)))


# Rank genes by score
potNovel_genes <- potNovel_genes %>%
  arrange(desc(Total_Score))

# View the ranked genes
print(potNovel_genes)

# Score for all genes ----
# calculate min-max values per biotype
gene_level_data_score <- gene_level_data %>% group_by(biotype) %>%
  mutate(
    Expression_Norm = min_max_normalize(log10(max_mean_tpm+1)),
    Exonic_Conservation_Norm = min_max_normalize(phastCons_mean_merged_exons),
    Promoter_Conservation_Norm = min_max_normalize(phastCons_most_conserved_promoter),
    TF_Binding_Norm = min_max_normalize(N_TFs)
  )

gene_level_data_score <- gene_level_data_score %>%
  mutate(Total_Score= CAGE_within_100bp*3 + Expression_Norm*3 +
           ifelse(max_corr_signature >= 0.6,2,0) +

           TF_Binding_Norm * 1 +
           H3K4me3_at_promoter * 1 + H3K27ac_at_promoter * 1 +
           H3K36me3_at_geneBody * 0.5 +
           Enhancer_Atlas_at_promoter * 0.5 + FANTOM_enhancers_at_promoter  * 0.5 +
           Enhancer_from_marks_at_promoter * 0.5 + Enhancer_from_ABC *1 +
           Promoter_Conservation_Norm * 1 +
           Exonic_Conservation_Norm * 2 +
           (ifelse(width > 500, 0.5, 0)) )

gene_level_data <- left_join(gene_level_data,gene_level_data_score%>% dplyr::select(gene_name, Total_Score))
ggplot(gene_level_data,aes(exonic_type, fill = biotype,y=Total_Score)) +
  geom_boxplot() + theme_minimal() + ggtitle("Total Score distribution")

table(gene_level_data$biotype,gene_level_data$exonic_type,gene_level_data$Total_Score>=5)

# other ideas: just filter more by expression, minimum length of 500...
# minimum TF binding of 5??
# just derive a confident list by the end of the study...

# sort by
# high priority, expression, midpriority

# perform kallisto aggregation by new transcriptome def ----
kallisto_tl_path <- "outputs/expression_data/kallisto/filtered_assembly_20250127_112131_StemLinc_samples_kallisto.counts.tsv"
kallisto_tl_path <- "outputs/expression_data/kallisto/LSK_StemLinc.combined.filtered.20250127_112131.gene_name.subtracted_reps.counts.tsv"
kallisto_tl_counts <- read.delim(kallisto_tl_path)

filtered_annot <- read.delim("outputs/transcriptome_characterization/jan25/LSK_StemLinc.combined_annotated_tracking.filtered.sel_biots.20250127_112131.tsv")

kallisto_tl_counts <- kallisto_tl_counts %>% filter(transcript_id%in%filtered_annot$V1)
kallisto_tl_counts$gene_name <- filtered_annot$gene_name[match(kallisto_tl_counts$transcript_id,
                                                               filtered_annot$V1)]

# aggregate at gene level
kallisto_counts <- aggregate(kallisto_tl_counts[,2:(ncol(kallisto_tl_counts)-1)],
                             by = list(gene_name=kallisto_tl_counts$gene_name),sum)

# perform DGEA ----
# DGEA ----

# Load the count matrix and metadata
gene_names <- kallisto_counts$gene_name
kallisto_counts <- kallisto_counts[,-1]
colnames(kallisto_counts) <- c(paste0("LSK_StL_",c(1:3)),
                               paste0("macrophage_StL_",c(1:3)),
                               paste0("T-cell_StL_",c(1:3)))
coldata <- data.frame(sample_id=colnames(kallisto_counts),
                      sample=rep(c("LSK","macrophage","T-cell"),each=3),
                      replicate=rep(1:3,3))
rownames(coldata) <- coldata$sample_id
# Ensure sample names in metadata match column names in count matrix
all(colnames(kallisto_counts) == rownames(coldata))

# convert counts to integer
kallisto_counts <- apply(kallisto_counts,2,as.integer)
rownames(kallisto_counts) <- gene_names
# generate DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = kallisto_counts, colData = coldata, design = ~ sample)

table(rowSums(counts(dds)) > 10)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
resultsNames(dds)
rowData(dds)$biotype <- gene_level_data$biotype[match(rownames(dds),
                                                      gene_level_data$gene_name)]

vsd <- vst(dds, blind = FALSE)

# Perform PCA and visualize:

pca_data <- plotPCA(vsd, intgroup = "sample", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = sample)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_minimal()

# PCA per biotype
biots=c("protein_coding","lncRNA","potNovel","TEC","pseudogene")

for (biot in biots) {
  pca_tmp <- plotPCA(vsd[rowData(dds)$biotype == biot,],
                     intgroup = "sample", returnData = TRUE)
  percent_var <- round(100 * attr(pca_tmp, "percentVar"))

  g = ggplot(pca_tmp, aes(PC1, PC2, color = sample)) +
    geom_point(size = 5) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    theme_minimal() + ggtitle(paste("Biotype :",biot))

  print(g)

}

# Sort by significance and plot heatmap for each biotype:
results_list=list()
for(biot in biots){

  # T-cell
  res_T.cellvsLSK <- results(dds[rowData(dds)$biotype == biot,], contrast = c("sample",
                                                                              "T.cell", "LSK"))
  res_name=paste0(biot,"T.cell")
  res_T.cellvsLSK <- res_T.cellvsLSK[order(res_T.cellvsLSK$padj), ]
  summary(res_T.cellvsLSK)
  res_T.cellvsLSK <- as.data.frame(res_T.cellvsLSK)
  res_T.cellvsLSK <- res_T.cellvsLSK %>% mutate(Diff=ifelse(padj<0.05,ifelse(log2FoldChange<0,"DOWN","UP"),"NS"))
  res_T.cellvsLSK <- res_T.cellvsLSK %>% filter(!is.na(padj))
  results_list[[res_name]] <- as.data.frame(res_T.cellvsLSK)
  #write.csv(as.data.frame(res), "differential_genes.csv")
  g <- ggplot(res_T.cellvsLSK, aes(x = log2FoldChange, y = -log10(padj), col=Diff)) +
    geom_point() +
    #geom_hline(y = -log10(sig_threshold), linetype = "dashed") +
    #geom_vline(x = c(-fc_threshold, fc_threshold), linetype = "dashed") +
    xlab("Log2 Fold Change") +
    ylab("-log10 Adjusted p-value") +
    ggtitle(paste("Volcano Plot T-cell vs LSK",biot)) +
    scale_color_manual(values = nejm_pal[c(2,6,1)]) +
    theme_bw()
  print(g)
  # macrophage
  res_macrophagevsLSK <- results(dds[rowData(dds)$biotype == biot,], contrast = c("sample",
                                                                                  "macrophage", "LSK"))
  res_name=paste0(biot,"macrophage")
  res_macrophagevsLSK <- res_macrophagevsLSK[order(res_macrophagevsLSK$padj), ]
  summary(res_macrophagevsLSK)
  res_macrophagevsLSK <- as.data.frame(res_macrophagevsLSK)
  res_macrophagevsLSK <- res_macrophagevsLSK %>% mutate(Diff=ifelse(padj<0.05,ifelse(log2FoldChange<0,"DOWN","UP"),"NS"))
  res_macrophagevsLSK <- res_macrophagevsLSK %>% filter(!is.na(padj))
  results_list[[res_name]] <- as.data.frame(res_macrophagevsLSK)
  #write.csv(as.data.frame(res), "differential_genes.csv")
  g <- ggplot(res_macrophagevsLSK, aes(x = log2FoldChange, y = -log10(padj), col=Diff)) +
    geom_point() +
    #geom_hline(y = -log10(sig_threshold), linetype = "dashed") +
    #geom_vline(x = c(-fc_threshold, fc_threshold), linetype = "dashed") +
    xlab("Log2 Fold Change") +
    ylab("-log10 Adjusted p-value") +
    ggtitle(paste("Volcano Plot macrophage vs LSK",biot)) +
    scale_color_manual(values = nejm_pal[c(2,6,1)]) +
    theme_bw()
  print(g)
  #Filter significant genes (adjusted p-value < 0.05):

  sig_genes1 <- res_T.cellvsLSK[which(res_T.cellvsLSK$padj < 0.05), ]
  sig_genes2 <- res_macrophagevsLSK[which(res_macrophagevsLSK$padj<0.05),]

  table(rownames(sig_genes1)%in%rownames(sig_genes2))
  sig_genes <- unique(rownames(sig_genes1),rownames(sig_genes2))

  #sig_genes1 <- as.data.frame(sig_genes1)
  # Step 5: Generate a Heatmap
  # Extract normalized counts for significant genes:


  normalized_counts <- counts(dds, normalized = TRUE)
  sig_gene_counts <- normalized_counts[sig_genes, ]

  #Scale data and plot heatmap:

  pheatmap(sig_gene_counts,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = FALSE,
           treeheight_row = 0,
           treeheight_col = 0,
           scale = "row",
           #annotation_col = coldata,
           main = paste("Heatmap of",biot, "Differentially Expressed Genes"))

}

res_LSKvsT.cell <- as.data.frame(results(dds, contrast = c("sample", "LSK","T.cell")))
res_LSKvsT.cell <- res_LSKvsT.cell%>%filter(!is.na(padj))
res_LSKvsT.cell <- res_LSKvsT.cell %>% mutate(Diff=ifelse(padj<0.05,
                                                          ifelse(log2FoldChange<0,
                                                                 "DOWN","UP"),"NS"))

res_LSKvsmacrophage <- as.data.frame(results(dds, contrast = c("sample", "LSK","macrophage")))
res_LSKvsmacrophage <- res_LSKvsmacrophage%>%filter(!is.na(padj))
res_LSKvsmacrophage <- res_LSKvsmacrophage %>% mutate(Diff=ifelse(padj<0.05,
                                                                  ifelse(log2FoldChange<0,
                                                                         "DOWN","UP"),"NS"))

gene_level_data$LSK_vs_T.cell <- res_LSKvsT.cell$Diff[match(gene_level_data$gene_name,
                                                            rownames(res_LSKvsT.cell))]

gene_level_data$LSK_vs_macrophage <- res_LSKvsmacrophage$Diff[match(gene_level_data$gene_name,
                                                                    rownames(res_LSKvsmacrophage))]

table(gene_level_data$LSK_vs_T.cell=="UP"&gene_level_data$LSK_vs_macrophage=="UP",
      gene_level_data$biotype,gene_level_data$exonic_type)

table(gene_level_data$LSK_vs_T.cell=="UP"&gene_level_data$LSK_vs_macrophage=="UP",
      gene_level_data$biotype,gene_level_data$Tier2)

# Analyse DE genes ----
gene_level_data <- gene_level_data %>% mutate(UP_vs=ifelse(LSK_vs_T.cell=="UP",
                                        ifelse(LSK_vs_macrophage=="UP",
                                               "both","T-cell"),
                                        ifelse(LSK_vs_macrophage=="UP","macrophage","none")))

gene_level_data$UP_vs[is.na(gene_level_data$UP_vs)] <- "none"
# Divide in: UP in both, UP in one,
# analyse:
### CAGE ----
ggplot(gene_level_data%>% filter(!is.na(CAGE_within_100bp)),
       aes(x=UP_vs,fill=CAGE_within_100bp==1)) + geom_bar(position = "fill")+
  theme_minimal() + facet_wrap(~biotype) +scale_fill_manual(values = nejm_pal) +
  ggtitle("Fraction of genes with CAGE within 100 bp from TSS")
### Tiers ----
ggplot(gene_level_data%>% filter(!is.na(Tier1)),
       aes(x=UP_vs,fill=Tier1)) + geom_bar(position = "fill")+
  theme_minimal() + facet_wrap(~biotype) +scale_fill_manual(values = nejm_pal) +
  ggtitle("Fraction of DE genes in Tier1")

ggplot(gene_level_data%>% filter(!is.na(Tier2)),
       aes(x=UP_vs,fill=Tier2)) + geom_bar(position = "fill")+
  theme_minimal() + facet_wrap(~biotype) +scale_fill_manual(values = nejm_pal) +
  ggtitle("Fraction of DE genes in Tier1")
### Score ----
ggplot(gene_level_data%>% filter(!is.na(Total_Score)),
       aes(x=UP_vs,y=Total_Score)) + geom_boxplot()+
  theme_minimal() + facet_wrap(~biotype) +scale_fill_manual(values = nejm_pal) +
  ggtitle("Score distribution")
### Expression ----
vst_df <- as.data.frame(assay(vsd))
mean_vst_LSK=rowMeans(vst_df[,1:3])
gene_level_data$mean_vst_LSK <- mean_vst_LSK[match(gene_level_data$gene_name,
                                                   names(mean_vst_LSK))]
ggplot(gene_level_data%>% filter(!is.na(mean_vst_LSK)),
       aes(x=UP_vs,y=mean_vst_LSK)) + geom_boxplot()+
  theme_minimal() + facet_wrap(~biotype) +scale_fill_manual(values = nejm_pal) +
  ggtitle("Mean LSK vst ")

### genomic distribution ----
ggplot(gene_level_data%>% filter(biotype!="protein_coding"&!is.na(best_classif_to_PCG)),
       aes(x=UP_vs,fill=best_classif_to_PCG)) + geom_bar(position = "fill")+
  theme_minimal() + facet_wrap(~biotype) +
  scale_fill_manual(values = colorRampPalette(nejm_pal)(length(unique(gene_level_data$best_classif_to_PCG)))) +
  ggtitle("Genomic position of DE genes")
### evidence level ----
ggplot(gene_level_data%>%filter(!is.na(max_mean_tpm)&!is.na(evidence_level)),
       aes(x=UP_vs,fill=evidence_level)) +
  geom_bar(color="black",position = "fill") + theme_minimal() +
  scale_fill_manual(values = c("white",nejm_pal[c(3,2,4)])) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) + facet_wrap(~biotype)



### correlation with signature ----
ggplot(gene_level_data%>% filter(!is.na(max_corr_signature)),
       aes(x=UP_vs,y=max_corr_signature)) + geom_boxplot()+
  theme_minimal() + facet_wrap(~biotype) +scale_fill_manual(values = nejm_pal) +
  ggtitle("Mean correlation with Signature genes ")
### distance to closest PCG ----
ggplot(gene_level_data%>% filter(!is.na(distance_closest_PCG)),
       aes(x=UP_vs,y=abs(distance_closest_PCG))) + geom_boxplot()+
  theme_minimal() + facet_wrap(~biotype) +scale_fill_manual(values = nejm_pal) +
  ggtitle("Distance to closest gene ") + scale_y_log10()

###
### NTFs ----
ggplot(gene_level_data%>% filter(!is.na(N_TFs)),
       aes(x=UP_vs,y=N_TFs)) + geom_boxplot()+
  theme_minimal() + facet_wrap(~biotype) +scale_fill_manual(values = nejm_pal) +
  ggtitle("# TFs bound at promoter")
### enhancer or target ----
ggplot(gene_level_data %>% filter(!is.na(Enhancer_from_ABC)),
       aes(x=UP_vs,fill=Enhancer_from_ABC)) +
  geom_bar(position = "fill")+
  theme_minimal() + facet_wrap(~biotype) +scale_fill_manual(values = nejm_pal) +
  ggtitle("Fraction of genes identified as harboring enhancer by ABC")

## harbors enhancer form ABC ----
ggplot(gene_level_data %>% filter(!is.na(enhancer_target)),aes(x=UP_vs,fill=enhancer_target)) +
  geom_bar(position = "fill")+
  theme_minimal() + facet_wrap(~biotype) +scale_fill_manual(values = nejm_pal) +
  ggtitle("Fraction of genes identified as enhancer targets by ABC")
## GO enrichment of closest genes ----
gene_level_data <- left_join(gene_level_data,gene_classif%>% dplyr::select(gene_name,best_closest_PCG))

potNovel_diff <- gene_level_data%>% filter(biotype=="potNovel"&UP_vs!="none")
potNovel_diff_closest_genes <- potNovel_diff$best_closest_PCG

length(unique(potNovel_diff_closest_genes))
go_enrichment <- enrichGO(gene = unique(potNovel_diff_closest_genes),
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)
# View the top GO terms
head(go_enrichment)

# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) +
  ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")


gene_level_data <- left_join(gene_level_data,gene_classif%>%
                               dplyr::select(gene_name,best_closest_PCG))

lncRNA_diff <- gene_level_data%>% filter(biotype=="lncRNA"&UP_vs!="none")
lncRNA_diff_closest_genes <- lncRNA_diff$best_closest_PCG

length(unique(lncRNA_diff_closest_genes))
go_enrichment <- enrichGO(gene = unique(lncRNA_diff_closest_genes),
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)
# View the top GO terms
head(go_enrichment)

# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) +
  ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")


## potNovel "none" category
potNovel_non_diff <- gene_level_data%>% filter(biotype=="potNovel"&UP_vs=="none")
potNovel_non_diff_closest_genes <- potNovel_non_diff$best_closest_PCG
length(unique(potNovel_non_diff_closest_genes))
go_enrichment <- enrichGO(gene = unique(potNovel_non_diff_closest_genes),
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)
# View the top GO terms
head(go_enrichment)

# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) +
  ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")



lncRNA_non_diff <- gene_level_data%>% filter(biotype=="lncRNA"&UP_vs=="none")
lncRNA_non_diff_closest_genes <- lncRNA_non_diff$best_closest_PCG
length(unique(lncRNA_non_diff_closest_genes))
go_enrichment <- enrichGO(gene = unique(lncRNA_non_diff_closest_genes),
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)
# View the top GO terms
head(go_enrichment)

# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) +
  ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")


# Non differentially expressed genes ----
# nor up nor down
non_diff_genes <- gene_level_data %>% filter(LSK_vs_T.cell=="NS"&LSK_vs_macrophage=="NS")

table(non_diff_genes$biotype)

closestPCG_to_non_diff_potNovel <- non_diff_genes%>%filter(biotype=="potNovel")%>%
  pull(best_closest_PCG)
length(unique(closestPCG_to_non_diff_potNovel))
go_enrichment <- enrichGO(gene = unique(closestPCG_to_non_diff_potNovel),
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)
# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) +
  ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")


# GO enrichment of DE genes most correlated with these
potNovel_diff_most_corr_genes <- potNovel_diff$maxcorrPCG

length(unique(potNovel_diff_most_corr_genes))
go_enrichment <- enrichGO(gene = unique(potNovel_diff_most_corr_genes),
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)
# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) +
  ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")

# GO enrichment of closeby genes
# later: GO enrichment of DE only in macroph, same T-cells

# perform correlation kallisto feature counts or kallisto vs StringTie ----
# genes with no evidence from marks:
genes_none <- gene_level_data %>% filter(evidence_level=="none")

potNovel_none <- genes_none %>% filter(biotype=="potNovel")
table(potNovel_none$max_mean_tpm>=0.5)
table(potNovel_none$H3K27me3_at_promoter)
table(potNovel_none$max_mean_tpm>=0.5,
      potNovel_none$N_TFs>5)
table(potNovel_none$enhancer_target,
      potNovel_none$highcons)

table(potNovel_none$highcons)
table(potNovel_none$highcons,
      potNovel_none$H3K27me3_at_promoter)
table(potNovel_none$highcons,
      potNovel_none$max_mean_tpm>0.5)

summary(potNovel_none$max_corr)
table(potNovel_none$Enhancer_from_ABC)
table(potNovel_none$best_classif_to_PCG)
summary(potNovel_none$max_mean_tpm)
table(potNovel_none$max_corr_signature>0.5)

# expression of enhancers harboring ----
# Transcription factors ----
# make combinations ----

# Repeat overlap ----

# expression by biotype and class ----


