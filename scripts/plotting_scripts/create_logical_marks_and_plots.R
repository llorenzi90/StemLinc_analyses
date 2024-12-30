# collect info on overlaps ----
library(tidyverse)
gene_level_info_path="outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.biotypes_of_interest.tsv"
CAGE_polyA <- "outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.biotypes_of_interest.last.CAGE_polyAsite.tsv"
histone_marks <- "outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.biotypes_of_interest.last.histone_marks.tsv"
enhancers <- "outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.biotypes_of_interest.last.enhancers.tsv"

CAGE_polyA=read.table(CAGE_polyA,header = T)
histone_marks=read.table(histone_marks,header = T)
enhancers=read.table(enhancers,header = T)

gene_level_info <- left_join(CAGE_polyA,histone_marks)
gene_level_info <- left_join(gene_level_info, enhancers)

# generate a logic table of presence or not of each feature ----

# criteria to define overlap with each mark ----
# CAGE and polyA : peak within 100 bp
RNA_marks_co=100
gene_level_info$CAGE_within_100bp <- abs(gene_level_info$closest_CAGE)<=RNA_marks_co
gene_level_info$polyAsite_within_100bp <- abs(gene_level_info$closest_polyA)<=RNA_marks_co

# H3K4me3 at promoters -----
# H3K4me3 is tipically found at promoter regions
# I will keep the genes with any overlap ate their promoters
# gene_level_info$H3K4me3_promoter_cov>0

gene_level_info$H3K4me3_at_promoter <- gene_level_info$H3K4me3_promoter_cov>0
table(gene_level_info$H3K4me3_at_promoter,gene_level_info$biotype)
table(gene_level_info$CAGE_within_100bp,gene_level_info$biotype,gene_level_info$H3K4me3_at_promoter)

# H3K27ac and H3Kme1 are tipically located at TSS of eRNAs,
# so I will use the same criteria that I used for H3Kme3

# H3K27ac at promoters ----
# this mark is associated both to promoters and to active enhancers
gene_level_info$H3K27ac_at_promoter <- gene_level_info$H3K27ac_promoter_cov>0

# H3K4me1 at promoters ----
# this mark is associated to active enhancers
gene_level_info$H3K4me1_at_promoter <- gene_level_info$H3K4me1_promoter_cov>0


# H3K36me3 at gene bodies ----
# H3K36me3 is tipically found at gene bodies
# therefore I will select for this mark if it covers 1/4 of gene body (not too strict)
gene_level_info$H3K36me3_at_geneBody <- gene_level_info$H3K36me3_total_cov >0.25

# Enhancer Atlas at promoters ----
gene_level_info$Enhancer_Atlas_at_promoter <- gene_level_info$Enhancer_Atlas_promoter_cov>0


# FANTOM enhancer at promoteres ----
gene_level_info$FANTOM_enhancers_at_promoter <- gene_level_info$FANTOM_enhancers_promoter_cov>0


# Enhancer defined by marks ----
# add enhancer definition form Luo et al.
# "we found that LncHSC-1 and LncHSC-2 are transcribed from enhancer regions, marked by histone H3 Lys27
# acetylation (H3K27ac) or histone H3 Lys4 monomethylation (H3K4me1)
# but not H3K4me3 and H3K27me3"

# before I had used this definition:

# ol_marks <- ol_marks%>%mutate(Enhancer_from_marks= ifelse((H3K27ac_geneBody>cutoff|
#                                                              H3K4me3_geneBody>cutoff)&
#                                                             !(H3K4me3_geneBody>0|H3K27me3_geneBody>0),T,F))
# ol_marks$Enhancer=ol_marks$Enhancer_from_marks|ol_marks$EnhancerAtlas|ol_marks$EnhancerFANTOM

# but now I will use the promoter region bc transcribed from enhancer regions sounds more like promoter...
# and because I have defined overlap for these marks at promoter level
gene_level_info$H3K27me3_at_promoter <- gene_level_info$H3K27me3_promoter_cov>0
gene_level_info <- gene_level_info %>% mutate(Enhancer_from_marks_at_promoter = ifelse((H3K27ac_at_promoter|H3K4me1_at_promoter)&
                                                                                         !(H3K4me3_at_promoter|H3K27me3_at_promoter),T,F))



# remove repressive mark K27me3
source("scripts/source_all_functions.R")


# for each feature make fraction plot
features=colnames(gene_level_info_logic_marks)[-1]
outdir="outputs/plots/logic_marks"
dir.create(outdir)
for (feat in features) {
  feat_var <- sym(feat)
  tmp_data <- gene_level_info %>% group_by(biotype) %>%
    summarise(nmark=sum({{feat_var}}),
              fraction=nmark/n())
  g=ggplot(tmp_data, aes(x=biotype, y=fraction, fill=biotype)) +
    geom_col( ) +
    ggtitle(feat) + scale_fill_manual(values = nejm_pal)
  print(g)
  save_tiff_svg(g,outdir = outdir,
                  filename = paste0("fraction_of_",feat,"_per_biotype"))

}

# add exonic type
for (feat in features) {
  feat_var <- sym(feat)
  tmp_data <- gene_level_info %>% group_by(biotype, exonic_type) %>%
    summarise(nmark=sum({{feat_var}}),
              fraction=nmark/n())
  g=ggplot(tmp_data, aes(x=biotype, y=fraction, fill=biotype)) +
    geom_col( ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle(feat) + scale_fill_manual(values = nejm_pal) +facet_wrap(~exonic_type)
  print(g)
  save_tiff_svg(g,outdir = outdir,
                filename = paste0("fraction_of_",feat,"_per_biotype_exonic_type"))

}


# combination of any mark at RNA and DNA level:
gene_level_info <- gene_level_info %>%
  mutate(RNA_level=CAGE_within_100bp|polyAsite_within_100bp,
         ,
         Enhancer=Enhancer_Atlas_at_promoter|FANTOM_enhancers_at_promoter|
                      Enhancer_from_marks_at_promoter,
         DNA_level=H3K4me3_at_promoter|H3K27ac_at_promoter|
           H3K4me1_at_promoter|H3K36me3_at_geneBody|Enhancer,
         both=RNA_level&(DNA_level|Enhancer),
         RNA=RNA_level&!both,
         DNA=DNA_level&!both,
         none=!(RNA|DNA|both)
         )

data2plot <- pivot_longer(gene_level_info,cols = c(RNA, DNA,both,none),
                           names_to = "evidence_level",
                          values_to = "evidence") %>% filter(evidence) %>%
  mutate(evidence_level=factor(evidence_level,levels=c("none","both","RNA","DNA")))

g=ggplot(data2plot,aes(biotype,fill=evidence_level)) +
  geom_bar(col="black")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("white","#E18727FF","#0072B5FF", "#20854EFF"))

g

save_tiff_svg(g,outdir = outdir,filename = "level_evidence_perbiotype")

g=ggplot(data2plot,aes(biotype,fill=evidence_level)) +
  geom_bar(col="black", position = "fill")+
  theme_minimal() + ylab("Fraction")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("white","#E18727FF","#0072B5FF", "#20854EFF"))

g

save_tiff_svg(g,outdir = outdir,filename = "level_evidence_perbiotype_fraction")




g=ggplot(data2plot %>%filter(evidence_level!="none"),aes(Enhancer,fill=evidence_level)) +
  geom_bar(col="black")+
  theme_minimal() + ylab("Fraction")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ biotype           )+
  ggtitle("Genes with evidence vs overlap Enhancer")+
  scale_fill_manual(values = c("#E18727FF","#0072B5FF", "#20854EFF"))

g
save_tiff_svg(g,outdir = outdir,filename = "level_evidence_perbiotype_andEnhancer")


g=ggplot(data2plot %>%filter(evidence_level!="none"),aes(Enhancer,fill=evidence_level)) +
  geom_bar(col="black",position = "fill")+
  theme_minimal() + ylab("Fraction")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ biotype           )+
  ggtitle("Genes with evidence vs overlap Enhancer")+
  scale_fill_manual(values = c("#E18727FF","#0072B5FF", "#20854EFF"))

g
save_tiff_svg(g,outdir = outdir,filename = "level_evidence_perbiotype_andEnhancer_fraction")


g=ggplot(data2plot %>%filter(evidence_level!="none"),aes(biotype,fill=Enhancer)) +
  geom_bar(col="black",position = "fill")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("#E18727FF","#0072B5FF"))

g

save_tiff_svg(g,outdir = outdir,filename = "genes_with_any_evidence_fraction_enhancer")



gene_level_info_logic_marks <- data2plot %>%select(c(1,88:102))
timest=get_timestamp(gene_level_info_path)
out_file <- paste0("outputs/overlap_marks/gene_level_info.",timest,"logicmarks.tsv")


# write data ----
write.table(gene_level_info_logic_marks,out_file,quote = F, row.names = F,sep = "\t")


