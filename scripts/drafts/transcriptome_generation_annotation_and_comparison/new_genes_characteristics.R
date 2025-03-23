## Head -------------------------------------
##
##
## Purpose of script: Summarize characteristics potential novel candidates
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-04-20
##
## Email: lucialorenzi90@gmail.com
##
## Notes ---------------------------
##
##
##   
##
## Setup ---------------------------

options(scipen = 999) 
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
library(ggplot2)
library(knitr)
library(DT)
library(ggsci)
library(tiff)
library(grid)
library(pheatmap)
## Load data---------------------------

# correlation with signature

vst=read.table("/home/llorenzi/work_local/expression_files/featureCounts/PCGslncRNAsPseudogenes.110424.vst.txt")
gene_biot=read.table("/home/llorenzi/work_local/expression_files/featureCounts/PCGslncRNAsPseudogenes.110424.annot.txt")
gene_level_info=read.delim("/home/llorenzi/work_local/APRIL/gene_level_info_APRIL.txt")
intramodular_connectivity=read.delim("/home/llorenzi/work_local/APRIL/WGCNA/intramodular_connectivity.txt")
table(gene_biot$simpl_biotype[!gene_biot$Geneid%in%gene_level_info$gene_id])
gene_translation_file <- "/home/llorenzi/references/conversion_table_ENSEMBL_NCBI_ids.csv"
gene_translation <- fread(gene_translation_file)
gene_biot$gene_name=gene_translation$gene_name[match(gene_biot$Geneid,
                                                     gene_translation$gene_id)]
vst_lncRNAs=vst[rownames(vst)%in%gene_biot$Geneid[gene_biot$simpl_biotype!="PCG"],]

gene_signature=read.table("/home/llorenzi/work_local/HSC_signature.txt")
#gene_signature=read.table("HSC_signature_Sergi.txt")
gene_sign_ids=gene_biot$Geneid[match(gene_signature$V1,gene_biot$gene_name)]
vst_signature=vst[rownames(vst)%in%gene_sign_ids,]

corr=cor(t(vst_lncRNAs),t(vst_signature))

corr=as.data.frame(corr)
corr$lncRNA=rownames(corr)
corr=pivot_longer(corr,cols=1:(ncol(corr)-1),
                  values_to = "corr",names_to = "PCG_id")
corr$PCG_name=gene_translation$gene_name[match(corr$PCG_id,gene_translation$gene_id_vM31)]
corr$lncRNA_name=gene_level_info$gene_name[match(corr$lncRNA,gene_level_info$gene_id)]

corr_summ=corr %>% arrange(-corr)%>%group_by(lncRNA)%>% 
  filter(any(corr>=0.8)) %>%
  summarise(N_0.8=sum(corr>=0.8),
            N_0.65=sum(corr>=0.65),
            PCGs=paste(PCG_name,collapse = ","),
            corr=paste(round(corr,2),collapse = ","))

corr_summ=left_join(corr_summ,gene_level_info%>%dplyr::select(gene_id,closest_gene_PCG,distance_PCG),
                    by=c(lncRNA="gene_id"))

corr_summ$lncRNA_name=gene_level_info$gene_name[match(corr_summ$lncRNA,
                                                      gene_level_info$gene_id)]

corr_summ=corr_summ%>%arrange(-N_0.8,-N_0.65)


corr_summ=left_join(corr_summ,intramodular_connectivity %>%
                      dplyr::select(gene,is_primary_corr,blue_module_membership,
                                    assigned_module_membership,kWithin,kTotal,assigned_module),
                    by=c(lncRNA="gene"))

corr_summ_extended=left_join(corr_summ,gene_level_info%>%
                               dplyr::select(gene_id,exonic_type,
                                             Coding_prob,
                                             meanPhastCons35way,
                                             closest_CAGE,
                                             dist_closest_enhancer,
                                             Delas_overlap,
                                             Delas_enrichment,
                                             Luo_gene),by = c(lncRNA="gene_id"))

counts=read.table("/home/llorenzi/work_local/expression_files/featureCounts/PCGslncRNAsPseudogenes.110424.counts.txt")
mean_counts_LSK_StL=apply(counts[,grep("LSK_StL",colnames(counts))],1,mean)

corr_summ_extended$mean_counts_LSK_StL=mean_counts_LSK_StL[match(corr_summ_extended$lncRNA,
                                                                 names(mean_counts_LSK_StL))]
mean_corr_signature=corr%>%group_by(lncRNA)%>%summarise(mean_corr=mean(corr))
mean_corr_signature=mean_corr_signature%>%arrange(-mean_corr)

View(mean_corr_signature%>%filter(mean_corr>0.645))
corr_summ_extended$mean_corr=mean_corr_signature$mean_corr[match(corr_summ_extended$lncRNA,
                                                                 mean_corr_signature$lncRNA)]

###
make_longer <- function(expression_data,add_gene_id=T,cols=2:ncol(expression_data),meas="tpm"){
  if(add_gene_id){
    expression_data$gene_id=rownames(expression_data)
    expression_data <- expression_data[,c(ncol(expression_data),1:(ncol(expression_data) - 1))]
  }
  expression_data <- pivot_longer(expression_data,cols = cols,names_to = "samples",values_to = meas)
  return(expression_data)
}

make_summ <- function(long_expdata,group_var,var,suff="tpm"){
  name1=paste0("mean_",suff)
  name2=paste0("max_",suff)
  long_expdata %>% group_by({{group_var}}) %>% 
    summarise("{name1}" := mean({{var}}),
              "{name2}" :=max({{var}}))
}


vst=make_longer(vst,meas = "vst")
vst_summ <- vst %>% make_summ(gene_id,vst,suff = "vst")



# Select potnovel candidates ----
potnovel=unique(c(grep("XLOC",corr_summ_extended$lncRNA,value = T),
                  grep("XLOC",mean_corr_signature$lncRNA[mean_corr_signature$mean_corr>0.645],value = T),
                  grep("XLOC",intramodular_connectivity$gene[intramodular_connectivity$assigned_module=="blue"&intramodular_connectivity$blue_module_membership>0],value = T)))

gene_level_info_potnovel=gene_level_info%>% filter(gene_id%in%potnovel)

prev_18=read.table("/home/llorenzi/work_local/potnovel_lncRNAs_info/gene_level_info_18potNovel.txt",header = T,sep = "\t")
table(prev_18$lncRNA%in%gene_level_info_potnovel$gene_id)


## Heatmap of expression 18 potNovel lncRNAs  
vst_pnew= vst %>% filter(gene_id%in% potnovel)

vst_hm=pivot_wider(vst_pnew %>% dplyr::select(gene_id,samples,vst),names_from = samples, values_from = vst)
vst_hm <- as.data.frame(vst_hm)

rownames(vst_hm)=vst_hm$gene_id

vst_hm=vst_hm[,-1]



pheatmap(vst_hm,fontsize_col=8, 
         fontsize_row=8,cluster_cols = F, scale = "none")

dir.create("/home/llorenzi/work_local/APRIL/plots_summary_candidates")
setwd("/home/llorenzi/work_local/APRIL/plots_summary_candidates")
tiff("heatmap_potNew_APRIL.tiff",res=300, width = 10,height = 10,units = "in")
pheatmap(vst_hm,fontsize_col=8, 
         fontsize_row=8,cluster_cols = F, scale = "none")

dev.off()

# add the previously defined genes
potnovel=unique(c(potnovel,prev_18$lncRNA))
gene_level_info_potnovel=gene_level_info%>% filter(gene_id%in%potnovel)

write.table(gene_level_info_potnovel,"/home/llorenzi/work_local/APRIL/selected_genes_all.27.txt",quote = F,row.names = F,sep = "\t")
vst_pnew= vst %>% filter(gene_id%in% potnovel)

vst_hm=pivot_wider(vst_pnew %>% dplyr::select(gene_id,samples,vst),names_from = samples, values_from = vst)
vst_hm <- as.data.frame(vst_hm)

rownames(vst_hm)=vst_hm$gene_id

vst_hm=vst_hm[,-1]



pheatmap(vst_hm,fontsize_col=8, 
         fontsize_row=8,cluster_cols = F, scale = "none")

tiff("heatmap_potNew_APRIL.27candidates.tiff",res=300, width = 10,height = 10,units = "in")
pheatmap(vst_hm,fontsize_col=8, 
         fontsize_row=8,cluster_cols = F, scale = "none")

dev.off()
# Plots categorical features   



discrete_vars=gene_level_info_potnovel %>% group_by(gene_id) %>%
  summarise(is_multiexonic=exonic_type!="monoexonic",
            tr_level_support=tr_level_support==3,
            inBlue_module=assigned_module=="blue",
            is_intergenic=best_class=="intergenic",
            primary_corr_0.6=is_primary_corr>0.6,
            lowcoding_prob=(Coding_prob<0.44|is.na(Coding_prob)),
            CAGE=closest_CAGE==0,
            Enhancer=dist_closest_enhancer==0,
            repeat_0.5=repeat.fraction>0.5,
            mean_corr_0.65=mean_corr>0.65)

discrete_vars <- as.data.frame(discrete_vars)

discrete_vars[is.na(discrete_vars)] <- FALSE

discrete_vars_mat=matrix(as.numeric(as.matrix(discrete_vars[,-1])),nrow = nrow(discrete_vars),ncol = (ncol(discrete_vars) - 1))

colnames(discrete_vars_mat)=colnames(discrete_vars)[-1]
rownames(discrete_vars_mat)=discrete_vars$gene_id

pheatmap(t(discrete_vars_mat),scale = "none",
         cluster_rows = F,color = c("grey","darkgreen"))

tiff("categorical_heatmap.tiff",res=300, width = 10,height = 10,units = "in")
pheatmap(t(discrete_vars_mat),scale = "none",
         cluster_rows = F,color = c("grey","darkgreen"))
dev.off()

# VST heatmap with categorical features      


new_annorow=discrete_vars
rownames(new_annorow)=discrete_vars$gene_id
new_annorow=new_annorow[match(rownames(vst_hm),rownames(new_annorow)),]
new_annorow[new_annorow==T]="yes"
new_annorow[new_annorow==F]="no"

new_annorow <- left_join(new_annorow,gene_level_info_potnovel%>%
                           dplyr::select(gene_id,
                                         mean_corr_signature,
                                         mean_counts_LSK,
                                         is_primary_corr,
                                         meanPhastCons35way,
                                         repeat.fraction,
                                         Coding_prob,
                                         best_class))

rownames(new_annorow)=new_annorow$gene_id
new_annorow=new_annorow[,c(2:4,18,12:17)]
new_annorow[is.na(new_annorow)]=0

colnames(new_annorow)[5]="mean_corr_signature"

pheatmap(vst_hm,
         fontsize_col=4, 
         annotation_row = new_annorow,
         fontsize_row=8,cluster_cols = F, 
         scale = "none",legend = F)

gene_level_info_potnovel$mean_counts_LSK=mean_counts_LSK_StL[match(gene_level_info_potnovel$gene_id,
                                                                   names(mean_counts_LSK_StL))]

ordered_potnov=gene_level_info_potnovel%>% arrange(-mean_corr_signature,
                                    -mean_counts_LSK,
                                    -is_primary_corr,
                                    repeat.fraction,
                                    Coding_prob)
View(ordered_potnov%>%filter(mean_counts_LSK>40)%>%dplyr::select(gene_id,best_class,
                                    exonic_type,
                                    exonic_length,
                                    mean_corr_signature,
                                    mean_counts_LSK,
                                    closest_gene_PCG,
                                    distance_PCG,
                                    repeat.fraction,
                                    Coding_prob,
                                    meanPhastCons35way))

tiff("heatmap_vst_with_continuous_and_categorical_feats.tiff",
     res = 300,height = 10,width = 10,units = "in")
pheatmap(vst_hm[match(ordered_potnov$gene_id,rownames(vst_hm)),],
         fontsize_col=4, 
         annotation_row = new_annorow,
         fontsize_row=8,cluster_cols = F, 
         cluster_rows = F,
         scale = "none",legend = F)
dev.off()

cols2remove=c(3,7,8)


tiff("heatmap_vst_with_continuous_and_categorical_feats.7feats.tiff",
     res = 300,height = 10,width = 10,units = "in")
pheatmap(vst_hm[match(ordered_potnov$gene_id,rownames(vst_hm)),],
         fontsize_col=4, 
         annotation_row = new_annorow[,-cols2remove],
         fontsize_row=8,cluster_cols = F, 
         cluster_rows = F,
         scale = "none",legend = F)
dev.off()

tiff("heatmap_vst_with_continuous_and_categorical_feats.7feats.scalerow.tiff",
     res = 300,height = 10,width = 10,units = "in")
pheatmap(vst_hm[match(ordered_potnov$gene_id,rownames(vst_hm)),],
         fontsize_col=4, 
         annotation_row = new_annorow[,-cols2remove],
         fontsize_row=8,cluster_cols = F, 
         cluster_rows = F,
         scale = "row",legend = F)
dev.off()

# Add Spehd
potnovelANDSpehd=c(potnovel,"ENSMUSG00000120454.1")
gene_level_info_potnovelANDSpehd=gene_level_info%>% filter(gene_id%in%potnovelANDSpehd)

vst_pnew= vst %>% filter(gene_id%in% potnovelANDSpehd)

vst_hm=pivot_wider(vst_pnew %>% dplyr::select(gene_id,samples,vst),names_from = samples, values_from = vst)
vst_hm <- as.data.frame(vst_hm)

rownames(vst_hm)=vst_hm$gene_id

vst_hm=vst_hm[,-1]

#rownames(vst_hm)[28]="Spehd"

discrete_vars=gene_level_info_potnovelANDSpehd %>% group_by(gene_id) %>%
  summarise(is_multiexonic=exonic_type!="monoexonic",
            tr_level_support=tr_level_support==3,
            inBlue_module=assigned_module=="blue",
            is_intergenic=best_class=="intergenic",
            primary_corr_0.6=is_primary_corr>0.6,
            lowcoding_prob=(Coding_prob<0.44|is.na(Coding_prob)),
            CAGE=closest_CAGE==0,
            Enhancer=dist_closest_enhancer==0,
            repeat_0.5=repeat.fraction>0.5,
            mean_corr_0.65=mean_corr_signature>0.65)

discrete_vars <- as.data.frame(discrete_vars)

discrete_vars[is.na(discrete_vars)] <- FALSE

discrete_vars_mat=matrix(as.numeric(as.matrix(discrete_vars[,-1])),
                         nrow = nrow(discrete_vars),
                         ncol = (ncol(discrete_vars) - 1))

colnames(discrete_vars_mat)=colnames(discrete_vars)[-1]
rownames(discrete_vars_mat)=discrete_vars$gene_id

new_annorow=discrete_vars
rownames(new_annorow)=discrete_vars$gene_id
new_annorow=new_annorow[match(rownames(vst_hm),rownames(new_annorow)),]
new_annorow[new_annorow==T]="yes"
new_annorow[new_annorow==F]="no"
gene_level_info_potnovelANDSpehd$mean_counts_LSK=
  mean_counts_LSK_StL[match(gene_level_info_potnovelANDSpehd$gene_id,
                                                                   names(mean_counts_LSK_StL))]

new_annorow <- left_join(new_annorow,gene_level_info_potnovelANDSpehd%>%
                           dplyr::select(gene_id,
                                         mean_corr_signature,
                                         mean_counts_LSK,
                                         repeat.fraction,
                                         Coding_prob,
                                         best_class))

rownames(new_annorow)=new_annorow$gene_id
new_annorow=new_annorow[,c(2:3,16,12:15)]
new_annorow[is.na(new_annorow)]=0
new_annorow$tr_level_support[28]="yes"

rownames(vst_hm)[28]="Spehd"
rownames(new_annorow)[28]="Spehd"
ordered_potnovANDSpehd=gene_level_info_potnovelANDSpehd%>%arrange(-mean_corr_signature,
                                                                  -mean_counts_LSK,
                                                                  -is_primary_corr,
                                                                  repeat.fraction,
                                                                  Coding_prob)

ordered_potnovANDSpehd$gene_id[1]="Spehd"
tiff("heatmap_vst_with_continuous_and_categorical_feats.7feats.Spehd.tiff",
     res = 300,height = 10,width = 10,units = "in")
pheatmap(vst_hm[match(ordered_potnovANDSpehd$gene_id,rownames(vst_hm)),],
         fontsize_col=6, 
         annotation_row = new_annorow,
         fontsize_row=8,cluster_cols = F, 
         cluster_rows = F,
         scale = "none",legend = F)
dev.off()
