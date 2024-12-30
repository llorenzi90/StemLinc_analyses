source("~/Rprojects/StemLinc_analyses/scripts/correlation_with_signature.R")

coldata <- read.table(coldata_path,header = T)

vst_long <- pivot_longer_from_matrix(vst,val2 = "vst",nam2 = "sample")

vst_long <- left_join(vst_long, coldata)

vst_long <- left_join(vst_long, gene_level_info%>%select(gene_name, exonic_type,biotype,classif))

vst_long <- vst_long %>% mutate(cell_class=factor(cell_class,levels=c("HSPC","progenitor","differentiated")))


## DGEA ----
library(DESeq2)
counts=read.table(counts_path,header = T)
all(counts$Geneid==gene_level_info$gene_name) # this has to be TRUE!

# remove gene id column
rownames(counts) <- counts$Geneid
counts <- counts[,-1]

rownames(coldata)=coldata$sample


dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design= ~ study + cell_class)

#dds <- DESeq(dds)
tstamp=get_timestamp(counts_path)
dir.create("outputs/DESeq2")

#saveRDS(dds,paste0("outputs/DESeq2/dds_DESeq2_all_blood_cells.",
 #                  tstamp,"RDS"))
dds <- readRDS(paste0("outputs/DESeq2/dds_DESeq2_all_blood_cells.",
                      tstamp,"RDS"))
resultsNames(dds)

contrast1=c("HSPC","differentiated")
contrast2=c("HSPC","progenitor")
contrast3=c("progenitor","differentiated")

# volcano plots ----
for (cont in list(contrast1,contrast2,contrast3)) {
  res <- results(dds,contrast = c("cell_class",
                                  cont))

  ### volcano plot signature ----
  res <- as.data.frame(res) %>%
    mutate(diff=ifelse(padj<0.05&!is.na(padj),
                       ifelse(log2FoldChange>0,"UP","DOWN"),"NS"))
  table(res$diff)

  res$gene_class=res$diff
  res$gene_name <- rownames(res)
  res$gene_class[res$gene_name%in%HSC_signature$V1]="fingerprint"

  # volcano plot with fingerprint ----
  vp=ggplot(res, aes(x = log2FoldChange , y = -log2(padj))) +
    geom_point(data = res %>%filter(gene_class!="fingerprint"),
               aes(colour = gene_class),size=3) +
    geom_point(data = res %>%filter(gene_class=="fingerprint"),
               aes(colour = gene_class),size=3) +
    scale_colour_manual(values = mycolors[c(2,4,6,1)])+
    #geom_hline(yintercept = -log2(0.05), linetype = "dashed") +
    #geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    labs(x = "Log2 Fold Change", y = "-log2(padj)")+
    theme_classic()+
    ggtitle(paste0(cont,collapse = " vs "))+
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA)
    )
  vp

  save_tiff_svg(plo = vp,outdir = outdir,filename = paste0("volcano_plot_fingerprint_",
                                                           paste0(cont,collapse = "_vs_")))

  res <- left_join(res, gene_level_info%>%select(gene_name,biotype,classif))
  res <- res%>% mutate(biotype = ifelse(diff=="NS",diff,biotype))

  # plot non coding biotypes
  vp=ggplot(res%>%filter(biotype!="protein_coding"), aes(x = log2FoldChange , y = -log2(padj))) +
    geom_point(aes(color = biotype),size=2) +
    scale_color_manual(values = mycolors[c(1,6,2,3,4)])+
    #geom_hline(yintercept = -log2(0.05), linetype = "dashed") +
    #geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    labs(x = "Log2 Fold Change", y = "-log2(padj)")+
    theme_classic()+
    ggtitle(paste0(paste0(cont,collapse = " vs ")," non-coding biotypes"))+
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA)
    )


  vp


  save_tiff_svg(plo = vp,outdir = outdir,filename = paste0("volcano_plot_nonCoding_biotypes_",
                                                           paste0(cont,collapse = "_vs_")))

  res <- res %>% mutate(classif=ifelse(diff=="NS",diff,classif))
  res$biotype=gene_level_info$biotype
  # plot non coding biotypes colouring by classif
  vp=ggplot(res%>%filter(biotype!="protein_coding"), aes(x = log2FoldChange , y = -log2(padj))) +
    geom_point(aes(color = classif, shape = biotype),size=2) +
    scale_color_manual(values = mycolors)+
    #geom_hline(yintercept = -log2(0.05), linetype = "dashed") +
    #geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    labs(x = "Log2 Fold Change", y = "-log2(padj)")+
    theme_classic()+
    ggtitle(paste0(paste0(cont,collapse = " vs ")," non-coding biotypes"))+
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA)
    )


  vp


  save_tiff_svg(plo = vp,outdir = outdir,filename = paste0("volcano_plot_nonCoding_biotypes_shapeclassif_",
                                                           paste0(cont,collapse = "_vs_")))

}


HSPC_vs_differentiated <- results(dds,contrast = c("cell_class",
                                                          contrast1))

res <- HSPC_vs_differentiated
res <- as.data.frame(res) %>%
  mutate(diff=ifelse(padj<0.05&!is.na(padj),
                     ifelse(log2FoldChange>0,"UP","DOWN"),"NS"))
res$gene_class=res$diff
res$gene_name <- rownames(res)
res$gene_class[res$gene_name%in%HSC_signature$V1]="fingerprint"
res <- left_join(res, gene_level_info%>%select(gene_name,biotype,classif,mean_corr_signature,max_exons,exonic_type))

# boxplot expression all genes across the 3 sample classes ----

### boxplot per cell class and biotype ----
g=ggplot(vst_long, aes(x=biotype, fill=cell_class, y=vst)) + geom_boxplot() +
  scale_fill_manual(values = mycolors[c(1,2,3)]) +theme_minimal()
g

save_tiff_svg(g, outdir = outdir,
              filename = "boxplot_expression_across_sample_classes_per_biotype")
### boxplot per cell class biotype and exonic type ----
g=ggplot(vst_long, aes(x=biotype, fill=cell_class, y=vst)) + geom_boxplot() +
  scale_fill_manual(values = mycolors[c(1,2,3)]) +theme_minimal() + facet_wrap(~exonic_type)
g
save_tiff_svg(g, outdir = outdir,
              filename = "boxplot_expression_across_sample_classes_per_biotype_exonic_type")

vst_long <- left_join(vst_long, res%>% select(gene_name,log2FoldChange))
colnames(vst_long)[ncol(vst_long)] <- "HSPC_vs_differentiated_log2FC"

ggplot(vst_long %>% filter(!duplicated(gene_name)), aes(x=biotype, y=HSPC_vs_differentiated_log2FC)) + geom_boxplot() +
  theme_minimal()

# heatmaps ----
vst_mean_expression <- vst_long %>% group_by(cell_class, gene_name, biotype) %>%
  summarise(mean_vst=mean(vst))

# function that we will need to cluster per group:
# arrange genes according to cluster within each biotype
order_based_on_gene_cluster <- function(data){
  data <- as.data.frame(data)
  ord <- hclust( dist(data[,-1], method = "euclidean"), method = "ward.D" )$order
  return(data[,1][ord])
}

### top expressed genes in HSPC ----
for (x in c(100,500,1000)) {
  exp_data <- vst_long
  genes2keep <- vst_mean_expression %>% filter(cell_class=="HSPC") %>%
    arrange(-mean_vst) %>% group_by(biotype) %>% slice_head(n = x)
  exp_data <- exp_data%>% filter(gene_name%in%genes2keep$gene_name)

  gene_order=c()
  for (biot in unique(exp_data$biotype)) {
    wide_data <- pivot_wider(exp_data%>% filter(biotype==biot),id_cols = gene_name,
                             names_from = sample,
                             values_from = vst)

    gene_order=c(gene_order,order_based_on_gene_cluster(wide_data))
  }
  exp_data <- exp_data %>% mutate(gene_name=factor(gene_name,levels=gene_order))
  g=ggplot(exp_data, aes(x = sample, y = gene_name, fill = vst)) +
    geom_tile() +
    facet_wrap(~ biotype, scales = "free_y") +  # Separate heatmaps for each biotype
    scale_fill_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = median(exp_data$vst, na.rm = TRUE)) +
    theme_minimal() +
    labs(x = "Samples", y = "Genes", fill = "vst") +
    theme(axis.text.x = element_text(angle = 90,
                                     vjust = 0.5 ,
                                     size =2)) +
    theme(axis.text.y = element_blank(),  # Remove y-axis text
          axis.ticks.y = element_blank())

  g
  save_tiff_svg(g,outdir = outdir,filename = paste0("heatmap_all_samples_",x,"topHSPC_per_biotype.vst"))
}

### top differential genes between different contrasts----
for (cont in list(contrast1,contrast2,contrast3)) {
  res <- results(dds,contrast = c("cell_class",
                                  cont))

  ###
  res <- as.data.frame(res) %>%
    mutate(diff=ifelse(padj<0.05&!is.na(padj),
                       ifelse(log2FoldChange>0,"UP","DOWN"),"NS"))
  table(res$diff)

  res$gene_class=res$diff
  res$gene_name <- rownames(res)
  res <- left_join(res,gene_level_info %>% select(gene_name,biotype))

  for (x in c(10,20,50)) {
    exp_data <- vst_long
    genes2keep <- res  %>%
      arrange(-log2FoldChange) %>% group_by(biotype) %>% slice_head(n = x)
    exp_data <- exp_data%>% filter(gene_name%in%genes2keep$gene_name)

    gene_order=c()
    for (biot in unique(exp_data$biotype)) {
      wide_data <- pivot_wider(exp_data%>% filter(biotype==biot),id_cols = gene_name,
                               names_from = sample,
                               values_from = vst)

      gene_order=c(gene_order,order_based_on_gene_cluster(wide_data))
    }
    exp_data <- exp_data %>% mutate(gene_name=factor(gene_name,levels=gene_order))
    g=ggplot(exp_data, aes(x = sample, y = gene_name, fill = vst)) +
      geom_tile() +
      facet_wrap(~ biotype, scales = "free_y") +  # Separate heatmaps for each biotype
      scale_fill_gradient2(low = "blue",
                           mid = "white",
                           high = "red",
                           midpoint = median(exp_data$vst, na.rm = TRUE)) +
      theme_minimal() +
      labs(x = "Samples", y = "Genes", fill = "vst") +
      theme(axis.text.x = element_text(angle = 90,
                                       vjust = 0.5 ,
                                       size =2)) +
      theme(axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank())

    g
    save_tiff_svg(g,outdir = outdir,filename = paste0("heatmap_all_samples_",x,"topFC",paste0(cont,collapse = "_vs_"),"_per_biotype.vst"))

    ### pheatmap version
    exp_data_wide <- pivot_wider(exp_data, id_cols = gene_name,names_from = sample,values_from = vst)
    exp_data_wide=as.data.frame(exp_data_wide)
    rownames(exp_data_wide)=exp_data_wide$gene_name

    annot_row=exp_data %>% filter(!duplicated(gene_name)) %>% select(gene_name, biotype)
    rnames=annot_row$gene_name
    annot_row <- data.frame(biotype=factor(annot_row$biotype,levels = unique(annot_row$biotype)))
    rownames(annot_row)=rnames

    annot_col=data.frame(cell_class=factor(coldata$cell_class,levels = unique(coldata$cell_class)))
    rownames(annot_col)=colnames(exp_data_wide)[-1]

    vst_colors=colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(100)
    ann_colors = list(

      cell_class = c(HSPC=mycolors[1],progenitor=mycolors[2],differentiated=mycolors[3])
      ,
      biotype = c(potNovel=mycolors[4],
                  lncRNA=mycolors[5],
                  TEC=mycolors[6],
                  pseudogene=mycolors[7],
                  protein_coding=mycolors[8])
    )
    p=pheatmap::pheatmap(exp_data_wide[,-1],
                         cluster_cols = F, cluster_rows = F ,
                         annotation_row = annot_row,show_rownames = F, scale = "none",
                         show_colnames = F, annotation_col = annot_col ,
                         annotation_colors = ann_colors,
                         color = vst_colors
    )
    save_tiff_svg(p,outdir = outdir,filename = paste0("heatmap_all_samples_",x,"topFC",paste0(cont,collapse = "_vs_"),"_per_biotype.pheatmap.vst"))
  }

}


# heatmap with mean across cell class

### top differential per cell class ----
for (cont in list(contrast1,contrast2,contrast3)) {
  res <- results(dds,contrast = c("cell_class",
                                  cont))

  ###
  res <- as.data.frame(res) %>%
    mutate(diff=ifelse(padj<0.05&!is.na(padj),
                       ifelse(log2FoldChange>0,"UP","DOWN"),"NS"))
  table(res$diff)

  res$gene_class=res$diff
  res$gene_name <- rownames(res)
  res <- left_join(res,gene_level_info %>% select(gene_name,biotype))

  for (x in c(10,20,50)) {
    exp_data <- vst_mean_expression
    genes2keep <- res  %>%
      arrange(-log2FoldChange) %>% group_by(biotype) %>% slice_head(n = x)

    exp_data <- exp_data%>% filter(gene_name%in%genes2keep$gene_name)


    gene_order=c()
    for (biot in unique(exp_data$biotype)) {
      wide_data <- pivot_wider(exp_data%>%
                                 filter(biotype==biot),id_cols = gene_name,
                               names_from = cell_class,
                               values_from = mean_vst)

      gene_order=c(gene_order,order_based_on_gene_cluster(wide_data))
    }

    exp_data <- exp_data %>% mutate(gene_name=factor(gene_name,levels=gene_order))
    g=ggplot(exp_data, aes(x = cell_class, y = gene_name, fill = mean_vst)) +
      geom_tile() +
      facet_wrap(~ biotype, scales = "free_y") +  # Separate heatmaps for each biotype
      scale_fill_gradient2(low = "blue",
                           mid = "white",
                           high = "red",
                           midpoint = median(exp_data$mean_vst, na.rm = TRUE)) +
      theme_minimal() +
      labs(x = "Samples", y = "Genes", fill = "mean_vst") +
      theme(axis.text.x = element_text(angle = 90,
                                       vjust = 0.5 ,
                                       size =8)) +
      theme(axis.text.y = element_blank(),  # Remove y-axis text
            axis.ticks.y = element_blank())

    g
    save_tiff_svg(g,outdir = outdir,filename = paste0("heatmap_cell_class_",x,"topFC",paste0(cont,collapse = "_vs_"),"_per_biotype.mean_vst"))

    ### pheatmap version
    exp_data_wide <- pivot_wider(exp_data, id_cols = gene_name,
                                 names_from = cell_class,values_from = mean_vst)
    exp_data_wide=as.data.frame(exp_data_wide)
    rownames(exp_data_wide)=exp_data_wide$gene_name

    annot_row=exp_data %>% filter(!duplicated(gene_name)) %>%
      select(gene_name, biotype)
    annot_row=as.data.frame(annot_row)
    annot_row=annot_row[!duplicated(annot_row$gene_name),]
    rnames=annot_row$gene_name
    annot_row <- data.frame(biotype=factor(annot_row$biotype,
                                           levels = unique(annot_row$biotype)))
    rownames(annot_row)=rnames


    vst_colors=colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(100)
    ann_colors = list(


      biotype = c(potNovel=mycolors[4],
                  lncRNA=mycolors[5],
                  TEC=mycolors[6],
                  pseudogene=mycolors[7],
                  protein_coding=mycolors[8])
    )
    exp_data_wide <- exp_data_wide %>% arrange(match(gene_name,gene_order))
    p=pheatmap::pheatmap(exp_data_wide[,-1],
                         cluster_cols = F, cluster_rows = F ,
                         annotation_row = annot_row,show_rownames = F, scale = "none",
                         show_colnames = F ,
                         annotation_colors = ann_colors,
                         color = vst_colors
    )
    save_tiff_svg(p,outdir = outdir,filename = paste0("heatmap_cell_class_",x,"topFC",paste0(cont,collapse = "_vs_"),"_per_biotype.pheatmap.mean_vst"))
  }

}


# genes highly expressed in LSK but lowly in other samples ----
# first define what is highly expressed in LSK

vst_StemLinc <- vst_long %>% filter(grepl("StL",sample))
vst_StemLinc <- left_join(vst_StemLinc,coldata)
g <- ggplot(vst_StemLinc,aes(x=sample_type, y = vst , fill=biotype)) +
  geom_boxplot() + theme_minimal() +scale_fill_manual(values = mycolors)
g


save_tiff_svg(plo = g,outdir = outdir,
                filename = "boxplot_StemLinc_samples_per_biotype.vst")

sample_types=c("LT.HSC","LSK",unique(grep("MPP",coldata$sample_type,value = T)))
vst_hspcs <- vst_long %>% filter(sample_type %in%sample_types)
vst_hspcs$sample_type[grepl("MPP",vst_hspcs$sample_type)] ="MPPs"
vst_hspcs$study=ifelse(vst_hspcs$study=="StL","StemLinc","other")
vst_hspcs <- vst_hspcs %>% group_by(gene_name,sample_type,study,biotype) %>%
  summarise(mean_vst=mean(vst))

g <- ggplot(vst_hspcs,aes(x=sample_type, y = mean_vst , fill=study)) +
  geom_boxplot() + theme_minimal() +scale_fill_manual(values = mycolors)+
  facet_wrap(~biotype)
g

save_tiff_svg(plo = g,outdir = outdir,
              filename = "boxplot_StemLinc_and_other_studies_hspcs.vst")

g <- ggplot(vst_StemLinc,aes(x=biotype, y = vst , fill=sample_type)) +
  geom_boxplot() + theme_minimal() +scale_fill_manual(values = mycolors)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g

save_tiff_svg(plo = g,outdir = outdir,
              filename = "boxplot_StemLinc_samples_per_biotype2.vst")


vst_StemLinc <- vst_long %>% filter(grepl("StL",sample))
vst_StemLinc <- left_join(vst_StemLinc,coldata)
g <- ggplot(vst_StemLinc,aes(x=vst, col=sample_type)) +
  geom_density() + theme_minimal() +scale_color_manual(values = mycolors) +
  facet_wrap(~biotype)

g

save_tiff_svg(g, outdir ,filename = "expression_of_LSK_genes_in_other_StemLinc_samples")



g <- ggplot(vst_StemLinc,aes(x=sample_type, y = vst ,fill=sample_type)) +
  geom_boxplot() + theme_minimal() +scale_fill_manual(values = mycolors) +
  facet_wrap(~biotype)
g

save_tiff_svg(plo = g,outdir = outdir,
              filename = "boxplot_StemLinc_samples_facet_per_biotype.vst")


### add other non-polyA samples
#vst_long <- left_join(vst_long, coldata%>% select(sample,sample_type,study))

non_polyA_studies=c("StL","Cua")
vst_non_polyA <- vst_long %>% filter(study%in%non_polyA_studies)

g <- ggplot(vst_non_polyA,aes(x=sample_type, y = vst , fill=biotype)) +
  geom_boxplot() + theme_minimal() +scale_fill_manual(values = mycolors)
g

save_tiff_svg(plo = g,outdir = outdir,
              filename = "boxplot_non_polyA_samples_per_biotype.vst")


g <- ggplot(vst_non_polyA,aes(x=vst, col=sample_type)) +
  geom_density() + theme_minimal() +scale_color_manual(values = mycolors) +
  facet_wrap(~biotype)
g

save_tiff_svg(g, outdir ,filename = "expression_of_LSK_genes_in_other_non_polyA_samples")


g <- ggplot(vst_non_polyA,aes(x=sample_type, y = vst ,fill=sample_type)) +
  geom_boxplot() + theme_minimal() +scale_fill_manual(values = mycolors) +
  facet_wrap(~biotype)
g

save_tiff_svg(plo = g,outdir = outdir,
              filename = "boxplot_non_polyA_samples_facet_per_biotype.vst")


# add other LSK samples
vst_current <- vst_long %>%filter(study%in%non_polyA_studies|sample_type=="LSK")

g <- ggplot(vst_current,aes(x=interaction(sample_type, study), y = vst , fill=biotype)) +
  geom_boxplot() + theme_minimal() +scale_fill_manual(values = mycolors)
g

save_tiff_svg(plo = g,outdir = outdir,
              filename = "boxplot_non_polyA_and_otherLSK_samples_per_biotype.vst")


g <- ggplot(vst_current,aes(x=vst, col=interaction(sample_type,study))) +
  geom_density() + theme_minimal() +scale_color_manual(values = mycolors) +
  facet_wrap(~biotype)
g

save_tiff_svg(g, outdir ,filename = "expression_of_LSK_genes_in_non_polyA_and_otherLSK_samples")


g <- ggplot(vst_current,aes(x=interaction(sample_type, study), y = vst ,
                            fill=interaction(sample_type, study))) +
  geom_boxplot() + theme_minimal() +scale_fill_manual(values = mycolors) +
  facet_wrap(~biotype)
g

save_tiff_svg(plo = g,outdir = outdir,
              filename = "boxplot_LSK_genes_in_non_polyA_and_otherLSK_samples_facet_biotype.vst")


## add other LT-HSC samples

vst_current <- vst_long %>%filter(study%in%non_polyA_studies|
                                    sample_type%in%c("LSK","LT.HSC"))

g <- ggplot(vst_current,aes(x=interaction(sample_type, study), y = vst , fill=biotype)) +
  geom_boxplot() + theme_minimal() +scale_fill_manual(values = mycolors)
g

#save_tiff_svg(plo = g,outdir = outdir,
 #             filename = "boxplot_non_polyA_and_otherLSK_LT.HSCsamples_per_biotype.vst")


g <- ggplot(vst_current,aes(x=vst, col=interaction(sample_type,study))) +
  geom_density() + theme_minimal()  +
  facet_wrap(~biotype)
g

#save_tiff_svg(g, outdir ,filename = "expression_of_LSK_genes_in_non_polyA_and_otherLSK_LT.HSCsamples")


g <- ggplot(vst_current,aes(x=interaction(sample_type, study), y = vst ,
                            fill=interaction(sample_type, study))) +
  geom_boxplot() + theme_minimal() +
  facet_wrap(~biotype) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g

save_tiff_svg(plo = g,outdir = outdir,
              filename = "boxplot_LSK_genes_in_non_polyA_and_otherLSK_LT.HSCsamples_facet_biotype.vst")


# select genes that are highly expressed in all of our samples
# and compare to mean expression across other cells
# mean_expression LSK_StL, StL, LSK_polyA, LT.HSC_polyA
# max_expression other sample/study

mean_expression_LSK_StL <- vst_long %>% filter(sample_type=="LSK"&study=="StL") %>%
  group_by(gene_name) %>% summarise(mean_vst=mean(vst))
mean_expression_LSK_StL$group="LSK_StL"

mean_expression_StL <- vst_long %>% filter(study=="StL") %>%
  group_by(gene_name) %>% summarise(mean_vst=mean(vst))
mean_expression_StL$group="StL"

mean_expression_otherStL <- vst_long %>% filter(study=="StL"&sample_type!="LSK") %>%
  group_by(gene_name) %>% summarise(mean_vst=mean(vst))
mean_expression_otherStL$group="otherStL"

mean_expression_per_group <- rbind(mean_expression_LSK_StL,
                                   mean_expression_otherStL)


mean_expression_per_group_wide <- pivot_wider(mean_expression_per_group,
                                              id_cols = gene_name,
                                              names_from = group,
                                              values_from = mean_vst)

mean_expression_per_group_wide <- left_join(mean_expression_per_group_wide,
                                            gene_level_info %>% select(gene_name,biotype))
g=ggplot(mean_expression_per_group_wide%>%filter(biotype!="protein_coding"),aes(x=LSK_StL,y=otherStL ,col=biotype)) +
  geom_point() + theme_minimal() + scale_color_manual(values = mycolors)

g
save_tiff_svg(g,outdir = outdir,filename = "correlation_mean_vst_ncGenes_LSK_StemLinc_vs_other_StemLinc")

mean_expression_otherLSK_LT.HSC <- vst_long %>% filter(study!="StL"&sample_type%in%c("LSK","LT.HSC") )%>%
  group_by(gene_name,study,sample_type) %>% summarise(mean_vst=mean(vst))
mean_expression_otherStL$group="otherStL"

mean_expression_LSK_Kli <- vst_long %>% filter(study=="Kli"&sample_type=="LSK")  %>%
  group_by(gene_name) %>% summarise(mean_vst=mean(vst))
mean_expression_LSK_Kli$group="LSK_Kli"



mean_expression_per_group <- rbind(mean_expression_LSK_StL,
                                   mean_expression_LSK_Kli)


mean_expression_per_group_wide <- pivot_wider(mean_expression_per_group,
                                              id_cols = gene_name,
                                              names_from = group,
                                              values_from = mean_vst)

mean_expression_per_group_wide <- left_join(mean_expression_per_group_wide,
                                            gene_level_info %>% select(gene_name,biotype))
g=ggplot(mean_expression_per_group_wide%>%
           filter(biotype!="protein_coding"),aes(x=LSK_StL,y=LSK_Kli ,col=biotype)) +
  geom_point() + theme_minimal() + scale_color_manual(values = mycolors)

g
save_tiff_svg(g,outdir = outdir,filename = "correlation_mean_vst_ncGenes_LSK_StemLinc_vs_other_StemLinc")

# select genes above 8 vst in StemLinc and less than 5 in Klimmeck
mean_expression_per_group_wide_highStemLinc_lowKlimmeck <-
  mean_expression_per_group_wide %>% filter(LSK_StL>8&LSK_Kli<5)

# histone genes
histones_H1 <- paste0("H1f",1:5)
histones_H2 <- grep("-ps",
                    grep("^H2ac|H2bc",
                         gene_level_info$gene_name,value = T),
                    invert = T,value = T)
histones_H3 <- grep("^H3c",gene_level_info$gene_name,value = T)
histones_H3 <- c(histones_H3,"H3f3a","H3f3b")
histones_H4 <- grep("^H4c",gene_level_info$gene_name,value = T)

histones <- c(histones_H1,histones_H2,histones_H3,histones_H4)
histones
#write.table(histones,"outputs/histone_genes_5Oct.txt",quote = F,col.names = F,row.names = F)
View(mean_expression_per_group_wide%>% filter(gene_name%in%histones))
mean_expression_per_group_wide <- mean_expression_per_group_wide %>%
  mutate(FC=LSK_StL/LSK_Kli)
mean_expression_per_group_wide <- mean_expression_per_group_wide %>%
  mutate(biotype=ifelse(gene_name%in%histones,"histone",biotype))

g=ggplot(mean_expression_per_group_wide,aes(biotype,FC)) +
  geom_boxplot() + theme_minimal() +
  ggtitle("FC vst expression LSK StemLinc vs LSK Klimmeck") +
  ylab("FC totalRNA - polyA")

g

save_tiff_svg(g,outdir,"boxplot_FC_SL_Kli_per_biotype_histones")

median(mean_expression_per_group_wide%>%filter(biotype=="histone")%>%
         pull(FC))

median_histone_FC <- median(mean_expression_per_group_wide%>%
         filter(biotype=="histone")%>%
         pull(FC))

mean_expression_per_group_wide %>% group_by(biotype) %>%
  summarise(sum(FC>=median_histone_FC))

# mean_expression_per_group_wide %>% group_by(biotype) %>%
# +   summarise(sum(FC>=median_histone_FC))
# # A tibble: 6 × 2
# biotype        `sum(FC >= median_histone_FC)`
# <chr>                                   <int>
#   1 TEC                                        78
# 2 histone                                    33
# 3 lncRNA                                    101
# 4 potNovel                                  539
# 5 protein_coding                            282
# 6 pseudogene                                 21


# potential novel genes are indeed enriched for high FC!

## add other groups
vst_long <- vst_long %>% mutate(sample_study=
                                  paste(sample_type,
                                        study,sep = "_"))
mean_expression_per_study_sample_type <- vst_long %>%
  group_by(gene_name,sample_study) %>% summarise(mean_vst=mean(vst))

mean_expression_per_study_sample_type_wide <-
  pivot_wider(mean_expression_per_study_sample_type,id_cols = gene_name,
              names_from = sample_study,
              values_from = mean_vst)

table(mean_expression_LSK_StL$gene_name==
        mean_expression_per_study_sample_type_wide$gene_name)
mean_expression_per_study_sample_type_wide_FC <- apply(
  mean_expression_per_study_sample_type_wide[,-1],2,
  function(x)mean_expression_LSK_StL$mean_vst/x)

gene_biot=gene_level_info$biotype[match(mean_expression_per_study_sample_type_wide$gene_name,
                                  gene_level_info$gene_name)]
gene_biot[mean_expression_per_study_sample_type_wide$gene_name%in%histones]="histone"

#order columns
mean_expression_per_study_sample_type_wide_FC <- mean_expression_per_study_sample_type_wide_FC[,match(unique(vst_long$sample_study),
  colnames(mean_expression_per_study_sample_type_wide_FC))]
#p=pheatmap(mean_expression_per_study_sample_type_wide_FC[gene_biot=="potNovel",],
 #        cluster_cols = F)
# save_tiff_svg(p,outdir , "heatmap_FC_SL_LSK_vs_all_other_sample_study.potNovel")
# p=pheatmap(mean_expression_per_study_sample_type_wide_FC[gene_biot=="lncRNA",],
#          cluser_cols = F)
# save_tiff_svg(p,outdir , "heatmap_FC_SL_LSK_vs_all_other_sample_study.lncRNA")
#
# p=pheatmap(mean_expression_per_study_sample_type_wide_FC[gene_biot=="histone",],
#          cluster_cols = F)

#save_tiff_svg(p,outdir , "heatmap_FC_SL_LSK_vs_all_other_sample_study.histones")

set.seed(1234)
rand_proteins = sample(seq_len(sum(gene_biot=="protein_coding")),500,replace = F)
p=pheatmap(mean_expression_per_study_sample_type_wide_FC[which(gene_biot=="protein_coding")[rand_proteins],],
         cluster_cols = F)
save_tiff_svg(p,outdir , "heatmap_FC_SL_LSK_vs_all_other_sample_study.500_random_proteins")

# finally, for each sample calculate the mean FC for histone genes

mean_histones <- apply(mean_expression_per_study_sample_type_wide_FC[gene_biot=="histone",],2,
                       mean)
median_histones <- apply(mean_expression_per_study_sample_type_wide_FC[gene_biot=="histone",],2,
                       median)
# select genes that have FC higher than the histone median
higher_than_hist <- t(apply(mean_expression_per_study_sample_type_wide_FC,1,
                          function(x)x>=median_histones))

samps2remove <- unique(vst_long %>% filter(study%in%non_polyA_studies) %>% pull(sample_study))

higher_than_hist <- higher_than_hist[,!colnames(higher_than_hist)%in%samps2remove]
rownames(higher_than_hist) <- mean_expression_per_study_sample_type_wide$gene_name

Nhigher_than_hist=rowSums(higher_than_hist)
gene_level_info$Nhigher_than_hist <- Nhigher_than_hist[match(gene_level_info$gene_name,
                                                             names(Nhigher_than_hist))]

table(gene_level_info$Nhigher_than_hist>19,gene_level_info$biotype)
#write_info_table(gene_level_info,gene_level_data_path,modify_path = F)

# FC relative to HSPCs
HSPC_samples_polyA <- unique(vst_long%>%filter(cell_class=="HSPC"&!study%in%non_polyA_studies) %>% pull(sample_study))
mean_expression_HSPC_polyA <- mean_expression_per_study_sample_type_wide[,colnames(mean_expression_per_study_sample_type_wide)%in%HSPC_samples_polyA]
mean_expression_HSPC_polyA <- as.data.frame(mean_expression_HSPC_polyA)
rownames(mean_expression_HSPC_polyA) <- mean_expression_per_study_sample_type_wide$gene_name

mean_vst_HSPC_polyA <- apply(mean_expression_HSPC_polyA,1, mean)
mean_expression_LSK_StL
mean_vst_HSPC_polyA <- data.frame(gene_name=names(mean_vst_HSPC_polyA),
                                  mean_vst=mean_vst_HSPC_polyA,
                                    group="HSPC_polyA"
                                  )

mean_vst_LSK_and_HSPC_polyA=rbind(mean_vst_HSPC_polyA,mean_expression_LSK_StL)
mean_vst_LSK_and_HSPC_polyA_wide <- pivot_wider(mean_vst_LSK_and_HSPC_polyA,
                                                id_cols = gene_name,names_from = group,
                                                values_from = mean_vst)
mean_vst_LSK_and_HSPC_polyA_wide <- left_join(mean_vst_LSK_and_HSPC_polyA_wide,
                                              gene_level_info%>%select(gene_name,biotype,classif,Nhigher_than_hist))
g=ggplot(mean_vst_LSK_and_HSPC_polyA_wide %>% filter(biotype!="protein_coding"),
       aes(LSK_StL,HSPC_polyA, col=biotype)) +  geom_point() +
  theme_minimal() + scale_color_manual(values = mycolors)

g

save_tiff_svg(g,outdir = outdir,filename = "scatter_mean_vst_LSK_vs_HSPC_polyA_ncGenes")

g=g +geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")

g

save_tiff_svg(g,outdir = outdir,filename = "scatter_mean_vst_LSK_vs_HSPC_polyA_ncGenes_with_dashed_line")

mean_vst_LSK_and_HSPC_polyA_wide <- mean_vst_LSK_and_HSPC_polyA_wide %>% mutate(FC=LSK_StL/HSPC_polyA)

g=ggplot(mean_vst_LSK_and_HSPC_polyA_wide %>% filter(biotype!="protein_coding"),
         aes(LSK_StL,FC, col=biotype)) +  geom_point() +
  theme_minimal() + scale_color_manual(values = mycolors) + ylab("vst FC StemLinc vs mean_HSPCs")

g
save_tiff_svg(g,outdir = outdir,
              filename = "scatter_mean_vst_LSK_vs_FC_HSPC_samples_ncGenes")

gene_level_info <- left_join(gene_level_info ,
                             mean_vst_LSK_and_HSPC_polyA_wide %>%
                               select(gene_name,LSK_StL,HSPC_polyA,FC))

colnames(gene_level_info)[22:24] <- c("mean_vst_LSK_StL",
                                      "mean_vst_HSPC_polyA",
                                      "FC_vst_LSK_StL_vs_HSPC_polyA")

table(gene_level_info$FC_vst_LSK_StL_vs_HSPC_polyA>1.7,
      gene_level_info$biotype)
g=g + geom_hline(yintercept = 1.7)
save_tiff_svg(g,outdir = outdir,filename = "scatter_mean_vst_LSK_vs_FC_HSPC_samples_ncGenes_FCcutoff_1.7")


g+geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")

### add info on polyA site ----
## polyA ----

polyA_path="outputs/overlap_marks/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.closest_polyA.tsv"
polyA=read.table(polyA_path,header = T)

gene_level_info <- left_join(gene_level_info, polyA)


# distance closest polyA ----
polyA_co=500
histones=read.table("outputs/histone_genes_5Oct.txt")
histones=histones$V1
gene_level_info <- gene_level_info %>%
  mutate(polyA_class = ifelse(gene_name%in%histones,
                              "histone",
                              ifelse(biotype=="potNovel",
                                     ifelse(FC_vst_LSK_StL_vs_HSPC_polyA>1.7,
                                            "potNovel_highFC",
                                            ifelse(FC_vst_LSK_StL_vs_HSPC_polyA<1.2,
                                                   "potNovel_lowFC",
                                                   "potNovel_intermediate")),
                                     biotype)))

gene_level_info <- gene_level_info %>%
  mutate(polyA_biotype = ifelse(gene_name%in%histones,
                              "histone",
                              biotype),
         FC_class = ifelse(FC_vst_LSK_StL_vs_HSPC_polyA>1.7,
                           "highFC",
                           ifelse(FC_vst_LSK_StL_vs_HSPC_polyA<1.2,
                                  "lowFC",
                                  "intermediateFC")))

gene_level_info %>% plot_closest_mark(mark = closest_polyA,filter_var = polyA_class,
                                        filter_list = list(all=unique(gene_level_info)),
                                        outfile = "closest_polyA",save_plot = F)


g=ggplot(gene_level_info,
         aes(closest_polyA,col=biotype)) + geom_density(linewidth = 1) +
  xlim(c(-20,20)) +
  theme_minimal() +scale_color_manual(values = nejm_pal)
print(g)

save_tiff_svg(g,outdir,"polyAsite2.0_within20TES_all_biotypes")

g=ggplot(gene_level_info%>%filter(biotype%in%ncbiots),
         aes(closest_polyA,col=biotype)) + geom_density(linewidth = 1) +
  xlim(c(-20,20)) +
  theme_minimal() +scale_color_manual(values = nejm_pal)
print(g)
save_tiff_svg(g,outdir,"polyAsite2.0_within20TES_non-coding_biotypes")


g=ggplot(gene_level_info%>%filter(biotype%in%c("potNovel")),
         aes(closest_polyA,col=FC_class)) + geom_density(linewidth = 1) +
  xlim(c(-20,20)) +
  theme_minimal() +scale_color_manual(values = nejm_pal)
print(g)
save_tiff_svg(g,outdir,"polyAsite2.0_within20TES_potNovel_FC_class")




dist_co=10

frac_polyA_plot=gene_level_info %>%

  group_by(polyA_class) %>%
  summarise(nclose=sum(abs(closest_polyA)<=dist_co),
            fraction=nclose/n()) %>%
  ggplot(aes(x=polyA_class, y=fraction, fill=polyA_class)) +
  geom_col() +
  scale_fill_manual(values = nejm_pal) +
  ggtitle(paste0("Fraction of genes with polyA within ",dist_co," bp")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

frac_polyA_plot

for (dist_co in c(5,10,100,500)) {
  frac_polyA_plot=gene_level_info %>%

    group_by(polyA_class) %>%
    summarise(nclose=sum(abs(closest_polyA)<=dist_co),
              fraction=nclose/n()) %>%
    ggplot(aes(x=polyA_class, y=fraction, fill=polyA_class)) +
    geom_col() +
    scale_fill_manual(values = nejm_pal) +
    ggtitle(paste0("Fraction of genes with polyA within ",dist_co," bp")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

  frac_polyA_plot

  save_tiff_svg(frac_polyA_plot,outdir = outdir,
                filename = paste0("Fraction_polyA_within_",dist_co,"bp" ))

}



frac_polyA_plot=gene_level_info %>%

  group_by(polyA_biotype,FC_class) %>%
  summarise(nclose=sum(abs(closest_polyA)<=dist_co),
            fraction=nclose/n()) %>%
  ggplot(aes(x=FC_class, y=fraction, fill=FC_class)) +
  geom_col() +
  scale_fill_manual(values = nejm_pal) +
  ggtitle(paste0("Fraction of genes with polyA within ",dist_co," bp")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~polyA_biotype)

frac_polyA_plot


for (dist_co in c(5,10,100,500)) {
  frac_polyA_plot=gene_level_info %>%

    group_by(polyA_biotype,FC_class) %>%
    summarise(nclose=sum(abs(closest_polyA)<=dist_co),
              fraction=nclose/n()) %>%
    ggplot(aes(x=FC_class, y=fraction, fill=FC_class)) +
    geom_col() +
    scale_fill_manual(values = nejm_pal) +
    ggtitle(paste0("Fraction of genes with polyA within ",dist_co," bp")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~polyA_biotype)

  frac_polyA_plot
  save_tiff_svg(frac_polyA_plot,outdir = outdir,
                filename = paste0("Fraction_polyA_within_",dist_co,"bp_perFC_class_biotype" ))
}


for (dist_co in c(5,10,100,500)) {
  frac_polyA_plot=gene_level_info %>%

    group_by(polyA_biotype,FC_class) %>%
    summarise(nclose_polyAsignal=sum(abs(closest_polyA)<=dist_co&!is.na(polyA_signal)),
              fraction=nclose_polyAsignal/n()) %>%
    ggplot(aes(x=FC_class, y=fraction, fill=FC_class)) +
    geom_col() +
    scale_fill_manual(values = nejm_pal) +
    ggtitle(paste0("Fraction of genes with polyA signal within ",dist_co," bp")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~polyA_biotype)

  frac_polyA_plot
  save_tiff_svg(frac_polyA_plot,outdir = outdir,
                filename = paste0("Fraction_polyA_signal_within_",dist_co,"bp_perFC_class_biotype" ))
}


#write_info_table(gene_level_info,gene_level_data_path,modify_path = F)
#how do they compare to genes that are both highly expressed in our
# samples and in others samples?

# should I select only intergenic genes?



