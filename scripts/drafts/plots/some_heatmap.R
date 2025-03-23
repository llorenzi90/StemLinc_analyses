tt=pivot_wider(exp_data%>% filter(biotype=="potNovel"),id_cols = gene_name,names_from = sample,values_from = vst)
data <- scale(t(tt[,-1]))
ord <- hclust( dist(tt[,-1], method = "euclidean"), method = "ward.D" )$order
tt$gene_name[ord]

order_based_on_gene_cluster <- function(data){
  data <- as.data.frame(data)
  ord <- hclust( dist(data[,-1], method = "euclidean"), method = "ward.D" )$order
  return(data[,1][ord])
}

gene_order=c()
for (biot in unique(exp_data$biotype)) {
  wide_data <- pivot_wider(exp_data%>% filter(biotype==biot),id_cols = gene_name,
                           names_from = sample,
                           values_from = vst)

  gene_order=c(gene_order,order_based_on_gene_cluster(wide_data))
}
gene_order
exp_data <- exp_data %>% arrange(match(gene_order,gene_name))

exp_data_wide <- pivot_wider(exp_data, id_cols = gene_name,names_from = sample,values_from = vst)

annot_row=exp_data %>% filter(!duplicated(gene_name)) %>% select(gene_name, biotype)
annot_row=as.data.frame(annot_row)
rnames=annot_row$gene_name
annot_row <- data.frame(biotype=factor(annot_row$biotype,levels = unique(annot_row$biotype)))
rownames(annot_row)=rnames
exp_data_wide=as.data.frame(exp_data_wide)
rownames(exp_data_wide)=exp_data_wide$gene_name

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
save_tiff_svg(p,outdir = outdir,filename = "test")

"#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF" "#FFDC91FF" "#EE4C97FF"
