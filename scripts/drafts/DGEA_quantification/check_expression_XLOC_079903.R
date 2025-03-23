exp=read.table("~/Rprojects/StemLinc_analyses/outputs/expression_data/featureCounts_20240930_173216.ordered_samples.selected_biotypes.vst.txt",header = T)

library(ggplot2)
library(tidyverse)
xloc79903=exp %>% filter(Geneid=="XLOC_079903")
xloc79903=pivot_longer(xloc79903,cols = 2:ncol(xloc79903),names_to = "sample",values_to = "vst")
xloc79903$sample=factor(xloc79903$sample,levels = unique(xloc79903$sample))
ggplot(xloc79903, aes(sample,vst)) + geom_point() + theme(axis.text.x = element_text(angle = 90))


exp=read.table("~/Rprojects/StemLinc_analyses/outputs/expression_data/featureCounts_20240930_173216.ordered_samples.selected_biotypes.TPM.txt",header = T)


xloc79903=exp %>% filter(Geneid=="XLOC_079903")
xloc79903=pivot_longer(xloc79903,cols = 2:ncol(xloc79903),names_to = "sample",values_to = "TPM")
xloc79903$sample=factor(xloc79903$sample,levels = unique(xloc79903$sample))
ggplot(xloc79903, aes(sample,TPM)) + geom_point() +
  theme(axis.text.x = element_text(angle = 90)) + scale_y_log10()
