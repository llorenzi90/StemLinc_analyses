options(scipen = 999)
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
library(ggplot2)
library(rtracklayer)
library(ggpubr)
source("scripts/source_all_functions.R")

outdir="outputs/plots/CAGE_polyAsite"
dir.create(outdir)
gene_level_data_path <- "outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.biotypes_of_interest.last.tsv"
check_file_header(gene_level_data_path)
gene_level_info <- read.table(gene_level_data_path,header = T)
#gene_level_bed <- extract_bed(gene_level_info)

# define sets of biotypes of interest ----
biots_of_interest=c("protein_coding","lncRNA","TEC","potNovel","pseudogene")
ncbiots=c("lncRNA","TEC","potNovel","pseudogene")
lncRNA_potNovel_TEC=c("lncRNA","TEC","potNovel")
potNovel_TEC=c("TEC","potNovel")
biotypes2select=list(biots_of_interest,
                     ncbiots,
                     lncRNA_potNovel_TEC,
                     potNovel_TEC)

classes=sort(unique(gene_level_info$best_classif_to_PCG))

simpl_classes=c("antisense","convergent",
                "convergent","divergent","divergent",
                "intergenic","intronic_antisense",
                "intronic_sense","sense","sense","sense","sense",
                "sense")

#gene_level_info <- gene_level_info %>% mutate(classif=simpl_classes[match(best_classif_to_PCG,
                                                                  #classes)])
#write_info_table(gene_level_info,gene_level_data_path,suffix = "last")

#add exonic type ----
transcript_level_data_path="outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.gene_classif.tsv"
transcript_level_info <- read.table(transcript_level_data_path,header = T)
transcript_level_info <- transcript_level_info %>% filter(gene_name%in%gene_level_info$gene_name)
gene_exonic_type <- transcript_level_info %>% group_by(gene_name) %>%
  summarise(max_exons=max(Nexons),
            exonic_type=ifelse(any(Nexons>1),"multiexonic","monoexonic"))
#gene_level_info <- left_join(gene_level_info,gene_exonic_type)

# write_info_table(gene_level_info,
#                  gene_level_data_path,modify_path = F)

# functions ----
plot_closest_mark <- function(data=gene_level_info,
                              mark,
                              outdir="outputs/plots/",
                              outfile="plot",
                              filter_var= biotype,
                              col_var = {{filter_var}},
                              filter_list=biotypes2select,
                              xlims=c(-500,500),
                              color_pal=nejm_pal,
                              h=5,
                              w=8,
                              save_plot=T){
  for (biots in filter_list) {
    tmp_outfile=paste0(outfile,".",paste0(biots,collapse = "_"),"_",paste0(xlims,collapse = "_"))
    g=ggplot(data%>%filter({{filter_var}}%in%biots)
             ,
             aes({{mark}},col={{col_var}})) + geom_density(linewidth = 1) +
      xlim(xlims) +
      theme_minimal() +scale_color_manual(values = color_pal)
    print(g)
    if(save_plot) save_tiff_svg(g,outdir = outdir,filename = tmp_outfile, h=h,w=w)

  }
}


# CAGE ----
CAGE_path="outputs/overlap_marks/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.closest_CAGE.tsv"
CAGE=read.table(CAGE_path,header = T)

gene_level_info <- left_join(gene_level_info, CAGE)


# distance closest CAGE ----
gene_level_info %>% plot_closest_mark(mark = closest_CAGE,
                                      outfile = "closest_CAGE",
                                      outdir = outdir)


# by gene classif

gene_level_info %>% filter(biotype%in%ncbiots) %>%
  plot_closest_mark(mark = closest_CAGE,
                    outdir = outdir,
                    outfile = "CAGE_non_coding_biotypes",
                    filter_var = classif,
                    filter_list = list(unique(gene_level_info$classif)))

# Now plot fraction of genes with CAGE within X bp

dist_co=500
frac_cage_plot=gene_level_info %>%
  filter(biotype %in% biots_of_interest) %>%
  group_by(biotype) %>%
  summarise(nclose=sum(abs(closest_CAGE)<=dist_co),
            fraction=nclose/n()) %>%
  ggplot(aes(x=biotype, y=fraction, fill=biotype)) +
  geom_col() +
  scale_fill_manual(values = nejm_pal) +
  ggtitle(paste0("Fraction of genes with CAGE within ",dist_co," bp")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

frac_cage_plot

save_tiff_svg(frac_cage_plot,outdir = outdir,
              filename = paste0("Fraction_CAGE_within_",dist_co,"_TSS" ))

gene_level_info <- gene_level_info %>% mutate(CAGE_within_500 = abs(closest_CAGE) <= dist_co )


## polyA ----

polyA_path="outputs/overlap_marks/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.closest_polyA.tsv"
polyA=read.table(polyA_path,header = T)

gene_level_info <- left_join(gene_level_info, polyA)


# distance closest polyA ----
gene_level_info %>% plot_closest_mark(mark = closest_polyA,
                                      outdir = outdir,
                                      outfile = "closest_polyA")

# by gene classif

gene_level_info %>% filter(biotype%in%ncbiots) %>%
  plot_closest_mark(mark = closest_polyA,
                    outdir = outdir,
                    outfile = "polyA_non_coding_biotypes",
                    filter_var = classif,
                    filter_list = list(unique(gene_level_info$classif)))


# Now plot fraction of genes with polyA within X bp

dist_co=500
frac_polyA_plot=gene_level_info %>%
  filter(biotype %in% biots_of_interest) %>%
  group_by(biotype) %>%
  summarise(nclose=sum(abs(closest_polyA)<=dist_co),
            fraction=nclose/n()) %>%
  ggplot(aes(x=biotype, y=fraction, fill=biotype)) +
  geom_col() +
  scale_fill_manual(values = nejm_pal) +
  ggtitle(paste0("Fraction of genes with polyA within ",dist_co," bp"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

frac_polyA_plot

save_tiff_svg(frac_polyA_plot,
              outdir = outdir,
              filename = paste0("Fraction_polyA_within_",dist_co,"_TES" ))

gene_level_info <- gene_level_info %>% mutate(polyA_within_500 = abs(closest_polyA) <= dist_co )

# write extended gene level info ----
write_info_table(gene_level_info,data_path = gene_level_data_path,modify_path = T,
                 suffix = "CAGE_polyAsite")
