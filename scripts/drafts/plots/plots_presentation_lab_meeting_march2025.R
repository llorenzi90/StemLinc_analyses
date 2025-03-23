options(scipen = 999)
## Load packages---------------------------
require(tidyverse)
require(data.table)
require(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)

# Check if an input dir is provided
if (length(args) == 0) {
   stop("Please provide the working directory as a command-line argument.")
}

# Setup ----
dir_path <- args[1] # path to wdir
dir_path="outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation"
setwd(dir_path) #

# Load overlap with other catalogs ----
ol_prev_catalogs <- read.delim("overlap_other_catalogs/gene_level_data_lncRNAs_overlap.tsv")

# Load gene level info for expressed genes ----
gene_annot <- read.delim("expressed_merged_genes_annotation.sel_biotypes.assembled_col.tsv")

gene_annot <- left_join(gene_annot,ol_prev_catalogs%>%select(gene_id,Delas,Klimmeck,Luo,Qian,Nother_datasets))
gene_annot <- gene_annot%>%filter(assembled)
table(gene_annot$biotype,gene_annot$Nother_datasets)
gene_annot <- gene_annot %>% mutate(class=ifelse(biotype=="potNovel",
                                                 ifelse(Nother_datasets>0,"PreRep","StemLinc_only"),
                                                 biotype))
table(gene_annot$class)

source("functions/nejm_palette.R")

# Plot overlap with other catalogs ----
pie(table(gene_annot$class),col=nejm_pal)
pie(table(gene_annot$class[gene_annot$biotype!="protein_coding"]),col=nejm_pal)
pdf("overlap_other_catalogs/non_coding_genes_pie.pdf")
pie(table(gene_annot$class[gene_annot$biotype!="protein_coding"]),col=nejm_pal)
dev.off()

nc_genes <- gene_annot%>%filter(biotype!="protein_coding")
table(nc_genes$class,nc_genes$exonic_type)
ggplot(nc_genes,aes(x=class,fill=exonic_type)) + geom_bar() +theme_minimal() +
  scale_fill_manual(values = nejm_pal) + theme(text = element_text(size = 18))

# Load gene classification ----
classif <- read.delim("gene_classif/all_merged_genes.sel_biotypes.gene_classif.tsv")
table(nc_genes$gene_name%in%classif$gene_name)
gene_annot <- left_join(gene_annot,classif%>%select(gene_name,best_classif_to_PCG,best_closest_PCG,dist_closest_PCG))
nc_genes <- left_join(nc_genes,classif%>%select(gene_name,best_classif_to_PCG,best_closest_PCG,dist_closest_PCG))

table(nc_genes$best_classif_to_PCG)
simpl_classif=c("antisense","convergent","convergent","divergent","divergent","intergenic","intronic_antisense",
                "intronic_sense","sense_close","sense_overlap","sense_overlap","sense_overlap","sense_close")

names(simpl_classif) <- names(table(nc_genes$best_classif_to_PCG))
length(unique(simpl_classif))

nc_genes$classif_to_PCG <- simpl_classif[match(nc_genes$best_classif_to_PCG,
                                               names(simpl_classif))]

gene_annot$classif_to_PCG <- simpl_classif[match(gene_annot$best_classif_to_PCG,
                                                 names(simpl_classif))]
# Plot gene classification ----
ggplot(nc_genes,aes(x=class,fill=classif_to_PCG)) + geom_bar() +theme_minimal() + scale_fill_manual(values = nejm_pal) +
   theme(text = element_text(size = 18))

ggplot(nc_genes,aes(x=class,fill=classif_to_PCG)) + geom_bar(position = "fill") +theme_minimal() + scale_fill_manual(values = nejm_pal) +
  theme(text = element_text(size = 18)) + ylab("fraction")


# Plot distance to closest PCG ----
ggplot(nc_genes, aes(x=dist_closest_PCG+1, col=class)) + geom_density()+theme_minimal() +
  scale_x_log10()

table(nc_genes$dist_closest_PCG<=100000,nc_genes$class)

# polyA site
polyA_site <- read.delim("overlap_marks/all_merged_genes_annotation.sel_biotypes.closest_polyA.tsv")
table(gene_annot$gene_name%in%polyA_site$gene_name)
gene_annot <- left_join(gene_annot,polyA_site)

nc_genes <- gene_annot%>%filter(biotype!="protein_coding")

ggplot(nc_genes, aes(x=closest_polyA, col=class)) + geom_density()+ theme_minimal() +
  scale_x_log10()

ggplot(gene_annot, aes(x=closest_polyA, col=class)) + geom_density(linewidth = 2)+ theme_minimal() +
  scale_x_log10() + scale_color_manual(values = nejm_pal) + theme(text = element_text(size=18)) +
xlab("distance closest polyA site to TES")

ggplot(gene_annot, aes(y=closest_polyA, fill=class)) + geom_boxplot()+ theme_minimal() +
  scale_y_log10() + scale_fill_manual(values = nejm_pal) + theme(text = element_text(size=18)) +
  xlab("distance closest polyA site to TES")

gene_annot$polyA_site_within_100bp <- abs(gene_annot$closest_polyA)<=100

table(gene_annot$polyA_site_within_100bp,gene_annot$class)

ggplot(gene_annot,aes(x=class,fill=polyA_site_within_100bp)) +
  geom_bar(position = "fill") + theme_minimal() +scale_fill_manual(values = nejm_pal) +
  theme(text = element_text(size=18))
