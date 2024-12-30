############################################
##
## Purpose of script: plot gene classes
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-09-30
##
## Email: lucialorenzi90@gmail.com
##
#  Notes ---------------------------
##
##
##
#  Setup ---------------------------

options(scipen = 999)
library(RColorBrewer)
# Load packages---------------------------
# load up the packages we will need:  (uncomment as required)
library(ggplot2)
library(tidyverse)
source("scripts/source_all_functions.R")
# Load data---------------------------
gene_level_data_path <- "outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.tsv"
gene_level_data <- read.table(gene_level_data_path, header = T)

biots2keep <- c("lncRNA","TEC","potNovel","pseudogene")

plot_data <- gene_level_data %>% filter(biotype %in% biots2keep)

classes=sort(unique(plot_data$best_classif_to_PCG))

simpl_classes=c("antisense","convergent","convergent","divergent","divergent",
                "intergenic","intronic_antisense","intronic_sense","sense","sense","sense","sense",
                "sense")

plot_data <- plot_data %>% mutate(simpl_class=simpl_classes[match(best_classif_to_PCG,
                                                                  classes)])
# plots ----
g1 <- ggplot(plot_data, aes(x=biotype, fill = best_classif_to_PCG)) +
  geom_bar(position = "dodge") + scale_fill_manual(values = colorRampPalette(colors = nejm_pal)(length(classes)))
g1

g2 <- ggplot(plot_data, aes(x=biotype, fill = simpl_class)) +
  geom_bar(position = "dodge") + scale_fill_manual(values = nejm_pal)
g2

outdir="outputs/plots/"
save_tiff_svg(g1,outdir = outdir, filename = "non_coding_gene_classification.all", h = 8, w = 10)
save_tiff_svg(g2,outdir = outdir, filename = "non_coding_gene_classification.simple", h = 8, w = 10)
