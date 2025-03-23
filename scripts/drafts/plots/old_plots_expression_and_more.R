# "#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF" "#FFDC91FF" "#EE4C97FF"
# plot density with previously validated genes ----
g=ggplot(mean_corr_ext,aes(x=mean_corr))+
  geom_density(linewidth=2)+theme_classic() +
  geom_vline(xintercept =
               mean_corr_ext$mean_corr[mean_corr_ext$lncRNA_name%in%c("Spehd",
                                                                      "4930519L02Rik/Lnc-HSC-2")],
             col=mycolors[1],linewidth=0.5)+
  theme(text = element_text(size = 30)) +xlab("mean correlation")

g

save_tiff_svg(g,"tmp_plots/","density_correlation_signature",
              h=6,w=10)


g=ggplot(mean_corr_ext,aes(x=mean_corr))+
  geom_density(linewidth=2)+theme_classic() +
  theme(text = element_text(size = 30)) +xlab("mean correlation")

g

save_tiff_svg(g,"tmp_plots/","density_correlation_signature",
              h=6,w=10)




g=ggplot(mean_corr_ext,aes(x=mean_corr))+
  geom_density(linewidth=2)+theme_classic() +
  geom_vline(xintercept = 0.65,
             col=mycolors[1],linewidth=2)+
  theme(text = element_text(size = 30)) +xlab("mean correlation")

g

save_tiff_svg(g,"tmp_plots/","density_correlation_signature_cutoff",
              h=6,w=10)

mean(mean_corr$mean_corr[mean_corr$mean_corr>0])
mean(abs(mean_corr$mean_corr))
# prev studies defined as enriched in progenitors ----
mean_corr_ext_inprevstudies=mean_corr_ext%>%filter(!is.na(Delas_enrichment)|!is.na(Luo_gene))

### manually correct delas enrichment terms----
mean_corr_ext_inprevstudies$Delas_enrichment[is.na(mean_corr_ext_inprevstudies$Delas_enrichment)]="HSC_enriched"
mean_corr_ext_inprevstudies$
  Delas_enrichment[mean_corr_ext_inprevstudies$
                     Delas_enrichment=="AML_enriched,NA"]="AML_enriched"

table(mean_corr_ext_inprevstudies$Delas_enrichment)


mean_corr_ext_inprevstudies$
  Delas_enrichment[mean_corr_ext_inprevstudies$
                     Delas_enrichment=="NA,NA"]="HSC_enriched"


mean_corr_ext_inprevstudies$
  Delas_enrichment[grepl("HSC_enriched",mean_corr_ext_inprevstudies$
                           Delas_enrichment)]="HSC_enriched"

mean_corr_ext_inprevstudies$
  Delas_enrichment[grepl("lymphoid",mean_corr_ext_inprevstudies$
                           Delas_enrichment)]="lymphoid_enriched"


mean_corr_ext_inprevstudies$
  Delas_enrichment[grepl("progenitor_vs_differentiated",mean_corr_ext_inprevstudies$
                           Delas_enrichment)]="progenitor_vs_differentiated"

### plots ----
validated_genes=c("Spehd","4930519L02Rik/Lnc-HSC-2","Meg3")
validated_genes=c("Spehd","4930519L02Rik/Lnc-HSC-2")

mean_corr_ext_inprevstudies$validated=mean_corr_ext_inprevstudies$lncRNA_name%in%validated_genes
table(mean_corr_ext_inprevstudies$Delas_enrichment)

mean_corr_ext_inprevstudies$Delas_enrichment[mean_corr_ext_inprevstudies$Delas_enrichment=="progenitor_vs_differentiated"]="HSPC_enriched"
mean_corr_ext_inprevstudies$Delas_enrichment[mean_corr_ext_inprevstudies$Delas_enrichment=="HSC_enriched"]="HSPC_enriched"


test_cols=rep("lightgrey",nrow(mean_corr_ext_inprevstudies))
mean_corr_ext_inprevstudies$gene_name[mean_corr_ext_inprevstudies$validated]
test_cols[mean_corr_ext_inprevstudies$validated]=mycolors[1]
#test_cols[mean_corr_ext_inprevstudies$validated]=c(mycolors[1],mycolors[6],mycolors[1])

mean_corr_ext_inprevstudies$Delas_enrichment=factor(mean_corr_ext_inprevstudies$Delas_enrichment,
                                                    levels = c("lymphoid_enriched","AML_enriched",
                                                               "HSPC_enriched"))
g=ggplot(mean_corr_ext_inprevstudies,
         aes(x=Delas_enrichment,y=mean_corr)) +geom_violin() +
  theme_classic() +theme(text = element_text(size = 20),
                         axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") + ylab("mean correlation") +
  geom_jitter(col=test_cols,size=3)

g
save_tiff_svg(g,"tmp_plots/",
              "correlation_per_enriched_group_in_prev_studies.3groups.violin",
              h=6,w=10)

g=ggplot(mean_corr_ext_inprevstudies,
         aes(x=Delas_enrichment,y=mean_corr)) +geom_boxplot(outlier.shape = NA) +
  theme_classic() +theme(text = element_text(size = 20),
                         axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") + ylab("mean correlation") +geom_jitter(col=test_cols,size=3)
g
save_tiff_svg(g,
              "tmp_plots/",
              "correlation_per_enriched_group_in_prev_studies.3groups.1.boxplot",
              h=6,w=10)

mean_corr_ext$Delas_enrichment=mean_corr_ext_inprevstudies$Delas_enrichment[match(mean_corr_ext$gene_name,
                                                                                  mean_corr_ext_inprevstudies$gene_name)]

mean_corr_ext$Delas_enrichment[is.na(mean_corr_ext$Delas_enrichment)&mean_corr_ext$mean_corr>0.65]="HighCorr"
mean_corr_ext$Delas_enrichment[is.na(mean_corr_ext$Delas_enrichment)]="Not_reported"

mean_corr_ext$validated=mean_corr_ext$lncRNA_name%in%validated_genes

test_cols=rep("lightgrey",nrow(mean_corr_ext))
test_cols[mean_corr_ext$validated]=mycolors[1]

mean_corr_ext$Delas_enrichment=factor(mean_corr_ext$Delas_enrichment,
                                      levels = c("lymphoid_enriched",
                                                 "AML_enriched",
                                                 "Not_reported","progenitor_vs_differentiated",
                                                 "HSC_enriched","HighCorr"))
g=ggplot(mean_corr_ext,
         aes(x=Delas_enrichment,y=mean_corr)) +geom_boxplot(outlier.shape = NA) +
  theme_classic() +theme(text = element_text(size = 20),
                         axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") + ylab("mean correlation")+geom_jitter(col=test_cols,size=3)
g

save_tiff_svg(g,
              "tmp_plots/",
              "correlation_per_enriched_group_in_prev_studies.boxplot.with_notreported_and_Highcorr_lncRNAs",
              h=6,w=10)


mean_corr_ext_reported_only=mean_corr_ext%>%filter(Delas_enrichment!="Not_reported")
test_cols=rep("lightgrey",nrow(mean_corr_ext_reported_only))
test_cols[mean_corr_ext_reported_only$validated]=mycolors[1]
g=ggplot(mean_corr_ext%>%filter(Delas_enrichment!="Not_reported"),
         aes(x=Delas_enrichment,y=mean_corr)) +geom_boxplot(outlier.shape = NA) +
  theme_classic() +theme(text = element_text(size = 20),
                         axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") + ylab("mean correlation")+geom_jitter(col=test_cols,size=3)
g

save_tiff_svg(g,
              "tmp_plots/",
              "correlation_per_enriched_group_in_prev_studies.boxplot.and_Highcorr_lncRNAs",
              h=6,w=10)



mean_corr_ext_reported_only=mean_corr_ext%>%filter(Delas_enrichment!="Not_reported")
test_cols=rep("lightgrey",nrow(mean_corr_ext_reported_only))
test_cols[mean_corr_ext_reported_only$validated]=mycolors[1]
g=ggplot(mean_corr_ext%>%filter(Delas_enrichment!="Not_reported"),
         aes(x=Delas_enrichment,y=mean_corr)) +geom_boxplot(outlier.shape = NA) +
  theme_classic() +theme(text = element_text(size = 20),
                         axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") + ylab("mean correlation")+geom_jitter(col=test_cols,size=3)
g

save_tiff_svg(g,
              "tmp_plots/",
              "correlation_per_enriched_group_in_prev_studies.boxplot.and_Highcorr_lncRNAs",
              h=6,w=10)






### rm overlapping lncRNAs in sense or antisense ----

overlapping_classes=c("antisense","intronic_sense","sense_overlap","sense_overlap_downstream",
                      "sense_overlap_upstream")

mean_corr_ext_non_overlapping=mean_corr_ext%>%filter(!best_class%in%overlapping_classes)

mean_corr_ext_inprevstudies=mean_corr_ext_inprevstudies%>%filter(!best_class%in%overlapping_classes)

test_cols=rep("lightgrey",nrow(mean_corr_ext_inprevstudies))
test_cols[mean_corr_ext_inprevstudies$validated]= mycolors[1]

g=ggplot(mean_corr_ext_inprevstudies,
         aes(x=Delas_enrichment,y=mean_corr)) +geom_violin() +
  theme_classic() +theme(text = element_text(size = 20),
                         axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") + ylab("mean correlation") +
  geom_jitter(col=test_cols,size=3)
g
save_tiff_svg(g,"tmp_plots/",
              "correlation_per_enriched_group_in_prev_studies.violin.nonOlclasses",
              h=6,w=10)

g=ggplot(mean_corr_ext_inprevstudies,
         aes(x=Delas_enrichment,y=mean_corr)) +geom_boxplot(outlier.shape = NA) +
  theme_classic() +theme(text = element_text(size = 20),
                         axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") + ylab("mean correlation") +geom_jitter(col=test_cols,size=3)
g
save_tiff_svg(g,
              "tmp_plots/",
              "correlation_per_enriched_group_in_prev_studies.boxplot.nonOlclasses",
              h=6,w=10)


table(mean_corr_ext_non_overlapping$mean_corr>0.6)


max_corr=corr%>%arrange(-corr)%>%
  group_by(lncRNA_name)%>%
  summarise(max_corr=corr[1],
            max_correlated_gene=PCG_name[1])


mean_corr_ext_non_overlapping <- left_join(mean_corr_ext_non_overlapping,
                                           max_corr)

table(mean_corr_ext_non_overlapping$mean_corr>=0.6|mean_corr_ext_non_overlapping$max_corr>=0.8)


new_selected_genes=mean_corr_ext_non_overlapping%>%filter(mean_corr>=0.6,
                                                          max_corr>=0.8)

new_selected_genes_info=prev_sel%>%filter(gene_name%in%new_selected_genes$gene_name)


## expression levels vs other genes ----

gene_info_lncRNAs=gene_info%>%filter(biotype!="PCG")
gene_info_lncRNAs$Selected_genes=gene_info_lncRNAs$gene_name%in%new_selected_genes_info$gene_name

ggplot(gene_info_lncRNAs,aes(y=mean_LSK_TPM,x=Selected_genes)) +
  geom_boxplot() +scale_y_log10()

ggplot(gene_info_lncRNAs,aes(y=mean_phastCons35way,x=Selected_genes)) +
  geom_boxplot()


ggplot(gene_info_lncRNAs,aes(y=Coding_prob,x=Selected_genes)) +
  geom_boxplot()
table(gene_info_lncRNAs$best_class[gene_info_lncRNAs$Selected_genes],
      gene_info_lncRNAs$biotype[gene_info_lncRNAs$Selected_genes])

ggplot(gene_info_lncRNAs,aes(y=exonic_length,x=Selected_genes)) +
  geom_boxplot() +scale_y_log10()



cor_top_FC=cor(t(vst[rownames(vst)%in%(mean_corr_ext%>%arrange(-log2FoldChange)%>%dplyr::select(gene_name))$gene_name[1:10],]),
               t(vst[rownames(vst)%in%gene_info$gene_name[gene_info$biotype=="PCG"],]))

cor_top_FC
cor_top_FC=as.data.frame(cor_top_FC)
cor_top_FC$lncRNA=rownames(cor_top_FC)
cor_top_FC=pivot_longer(cor_top_FC,cols=1:(ncol(cor_top_FC)-1),
                        values_to = "cor_top_FC",names_to = "PCG_name")

cor_top_FC$lncRNA_name=prev_sel$Gname[match(cor_top_FC$lncRNA,prev_sel$gene_name)]
cor_top_FC$lncRNA_name[is.na(cor_top_FC$lncRNA_name)]=cor_top_FC$lncRNA[is.na(cor_top_FC$lncRNA_name)]

resultsNames(dds) # lists the coefficients
res <- results(dds, name="cell_class_primary_vs_differentiated")
res=as.data.frame(res)
res$gene_name=rownames(res)
new_selected_genes_info=left_join(new_selected_genes_info,res)
View(res%>%filter(gene_name%in%HSC_signature$V1))

mean_corr_ext=left_join(mean_corr_ext,res)

ggplot(mean_corr_ext,aes(x=mean_corr,y=log2FoldChange)) +geom_point(size=2) +
  geom_hline(yintercept = 7,col="darkred") + geom_vline(xintercept = 0.65,col=mycolors[1])

g=ggplot(mean_corr_ext,aes(x=mean_corr,y=log2FoldChange)) +geom_point(size=2) +
  geom_vline(xintercept = 0.65,col=mycolors[1],size=2) +theme_classic() +
  theme(text = element_text(size=20)) + ylab("log2FC HSPC vs differentiated")+
  xlab("mean correlation HSC fingerprint")

save_tiff_svg(g,"tmp_plots/","FC_vs_mean_corr_cutoff")


table(mean_corr_ext$padj<0.05&mean_corr_ext$log2FoldChange>5)
table((mean_corr_ext$padj<0.05&mean_corr_ext$log2FoldChange>7)|mean_corr_ext$mean_corr>=0.6)

sel_with_FC=mean_corr_ext%>%filter((padj<0.05 & log2FoldChange>7)|mean_corr>=0.6)
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="cell_class_primary_vs_differentiated", type="apeglm")


# plot expression fingerprint genes ----
pheatmap::pheatmap(vst[rownames(vst)%in%HSC_signature$V1,],cluster_cols = F,scale = "row")
vst_signature_long=pivot_longer_from_matrix(vst_signature,val2 = "vst",nam2 = "sample")

vst_signature_long$sample_type=coldata$cell_class[match(vst_signature_long$sample,
                                                        coldata$sample)]
vst_signature_long$sample=factor(vst_signature_long$sample,
                                 levels=unique(vst_signature_long$sample))

vst_signature_long$sample_type[vst_signature_long$sample_type=="primary"] = "HSPC"
vst_signature_long$sample_type=factor(vst_signature_long$sample_type,levels = c("HSPC","progenitor","differentiated"))
g=ggplot(vst_signature_long,
         aes(x = sample,y=vst,col=sample_type)) + geom_point() +
  facet_wrap(~gene_name) + theme_minimal() +
  scale_color_manual("cell type",values = mycolors)+theme(axis.text.x = element_blank(),
                                                          text = element_text(size = 25))+
  xlab("blood cells")

g

save_tiff_svg(g,"tmp_plots/",filename = "expression_signature_vst",
              h=6,w=10)

## public studies study overview ----
coldata$blood_cell=sapply(strsplit(coldata$sample,split = "_"),function(x)x[1])
colnames(coldata)[2]="cell_type"
coldata$cell_type[coldata$cell_type=="primary"]="HSPC"
coldata$cell_type=factor(coldata$cell_type,levels = c("HSPC","progenitor","differentiated"))
coldata$blood_cell=factor(coldata$blood_cell,levels = rev(unique(coldata$blood_cell)))

g=ggplot(coldata, aes(x=batch,y=blood_cell,col=cell_type)) +
  geom_point(size=5) +xlab("Study ID")+ theme_minimal() +
  theme(text = element_text(size=20)) +ylab("blood cell") +
  scale_color_manual(values=mycolors)

g

save_tiff_svg(g,"tmp_plots/",filename = "public_studies_samples_new",
              h=6,w=10)

# select candidates ----
# "Meg3" is dispensable for hematopoiesis
mean_corr_ext_selected=mean_corr_ext%>%filter(mean_corr>0.65)
table(mean_corr_ext_selected$annot)
new_selected_genes_info=left_join(mean_corr_ext_selected,gene_info%>%
                                    dplyr::select(gene_name,exonic_type,exonic_length))
table(new_selected_genes_info$best_class)
table(new_selected_genes_info$exonic_type,new_selected_genes_info$best_class%in%overlapping_classes)
table(new_selected_genes_info$exonic_type,new_selected_genes_info$annot)
table(new_selected_genes_info$exonic_type,new_selected_genes_info$best_class)

g=ggplot(new_selected_genes_info,aes(x=best_class, fill=exonic_type)) +
  geom_bar() + theme_classic() +
  scale_fill_manual(values=mycolors[c(6,2)]) + theme(text = element_text(size=20),
                                                     axis.text.x = element_text(angle = 45,
                                                                                hjust = 1))+
  xlab("") + ylab("# candidate lncRNA")
g
save_tiff_svg(g,"tmp_plots/","Ncandidates_by_class_exonix_type",h=6,w=10)


# exonic length ----

gene_info_lncRNAs$is_candidate=gene_info_lncRNAs$gene_name%in%new_selected_genes_info$gene_name

tmp_gene_info=gene_info
tmp_gene_info$biotype[tmp_gene_info$biotype%in%c("lncRNA","TEC")]="lncRNA"
tmp_gene_info$biotype=factor(tmp_gene_info$biotype,levels=c("lncRNA","PotNovel","PCG"))

median_Exonic_length=tmp_gene_info%>%
  group_by(biotype)%>%
  summarise(mdian_Exonic_length=median(exonic_length/1000),
            ypos=median(exonic_length/1000))
median_Exonic_length$lab=c("1.7 Kb","2.5 Kb","4.2 Kb")


g=ggplot(tmp_gene_info,aes(x=exonic_length,col=biotype))+
  stat_ecdf() +
  theme_classic() +
  ggtitle("Gene length")+scale_x_log10()

g

g=ggplot(tmp_gene_info,
         aes(x=biotype,fill=biotype,y=exonic_length/1000))+
  geom_violin() +
  theme_classic() + geom_text(data = median_Exonic_length,
                              aes(label=lab,y=ypos))+
  ggtitle("Gene exonic length")+scale_y_log10()+
  ylab("Exonic length (Kb)")+
  xlab("")+
  scale_fill_manual(values = mycolors)+
  theme(text=element_text(size=20))
g
save_tiff_svg(g,"tmp_plots/","length_dist_potnovel_lncRNA_PCG")


# mean phastcons ----

# overlap with marks ----
CAGE=read.table("data/overlap_marks/genes_with_CAGE_at_promoters.txt",header = T)
gene_info_lncRNAs$CAGE=gene_info_lncRNAs$gene_name%in%CAGE$gene_name
H3K27ac=read.table("data/overlap_marks/H3K27ac_gene_overlap.txt",header = T)
gene_info_lncRNAs$H3K27ac_fraction=
  H3K27ac$total_fraction_covered[match(
    gene_info_lncRNAs$gene_name,
    H3K27ac$gene_name)]

EA=read.table("data/overlap_marks/Enhancer_Atlas_gene_body.txt",header = T)
EF=read.table("data/overlap_marks/FANTOM_Enhancers_gene_body.txt",header = T)
gene_info_lncRNAs$EnhancerAtlas=EA$total_fraction_covered[match(
  gene_info_lncRNAs$gene_name,
  EA$gene_name)]

H3K27me3=read.delim("data/overlap_marks/H3K27me3_gene_overlap.txt")
H3K36me3=read.delim("data/overlap_marks/H3K36me3_gene_overlap.txt")
H3K4me3=read.delim("data/overlap_marks/H3K4me3_TSS_overlap.txt")
H3K4me1=read.delim("data/overlap_marks/H3K4me1_gene_overlap.txt")
H3K4me3_gene_body=read.delim("data/overlap_marks/H3K4me3_gene_overlap.txt")
gene_info_lncRNAs$EnhancerFANTOM=EF$total_fraction_covered[match(
  gene_info_lncRNAs$gene_name,
  EF$gene_name)]

gene_info_lncRNAs$H3K27me3_fraction=H3K27me3$total_fraction_covered[match(
  gene_info_lncRNAs$gene_name,H3K27me3$gene_name
)]
gene_info_lncRNAs$H3K4me3_atTSS=H3K4me3$TSS_ol_H3K4me3[match(
  gene_info_lncRNAs$gene_name,H3K4me3$gene_name
)]

gene_info_lncRNAs$H3K36_fraction=H3K36me3$total_fraction_covered[match(
  gene_info_lncRNAs$gene_name,H3K36me3$gene_name
)]

gene_info_lncRNAs$H3K4me1=H3K4me1$total_fraction_covered[match(
  gene_info_lncRNAs$gene_name,H3K4me1$gene_name
)]

gene_info_lncRNAs$H3K4me3_gene_body=H3K4me3_gene_body$total_fraction_covered[match(
  gene_info_lncRNAs$gene_name,H3K4me3_gene_body$gene_name
)]

gene_info_lncRNAs=gene_info_lncRNAs%>%mutate(possible_enhancer=
                                               ifelse((H3K27ac_fraction>0.2|
                                                         H3K4me1>0.2)&!(H3K4me3_gene_body>0|H3K27me3_fraction>0),T,F))
gene_info_lncRNAs=gene_info_lncRNAs%>%
  mutate(activation_signal=H3K4me3_atTSS&H3K36_fraction>0.3)

gene_info_lncRNAs=gene_info_lncRNAs%>%
  mutate(Enhancer_signal=EnhancerAtlas>0.4|EnhancerFANTOM>0.4|possible_enhancer)

summary_is_cand=gene_info_lncRNAs%>%group_by(is_candidate)%>%
  summarise(mean(mean_phastCons35way,na.rm=T),
            mean(repeat.fraction),
            mean(dist_closest_enhancer,na.rm=T),
            mean(mean_LSK_TPM),mean(exonic_length),
            sum(CAGE)/n(),
            mean(H3K27ac_fraction),
            sum(H3K27ac_fraction>0)/n(),
            mean(EnhancerAtlas),
            mean(EnhancerFANTOM),
            sum(EnhancerAtlas>0)/n(),
            sum(EnhancerFANTOM>0)/n(),
            mean(H3K27me3_fraction),
            sum(H3K4me3_atTSS)/n(),
            mean(H3K36_fraction),
            mean(H3K4me1),
            sum(activation_signal)/n()
  )

write.csv(gene_info_lncRNAs,"data/gene_info_lncRNA_candidates_0.65corr.olmarks.csv",row.names = F)
gene_info_lncRNAs_cand=gene_info_lncRNAs%>%filter(gene_info_lncRNAs$is_candidate)
gene_info_lncRNAs=read.csv("data/gene_info_lncRNA_candidates_0.65corr.olmarks.csv")

table(gene_info_lncRNAs_cand$EnhancerAtlas>0.2|gene_info_lncRNAs_cand$EnhancerFANTOM>0.2|gene_info_lncRNAs_cand$possible_enhancer)

gene_info_lncRNAs$H3K27ac=gene_info_lncRNAs$H3K27ac_fraction>0.2
table(gene_info_lncRNAs$H3K4me3_atTSS)
gene_info_lncRNAs$H3K27me3=gene_info_lncRNAs$H3K27me3_fraction>0.2
tmp_data_ol=gene_info_lncRNAs%>%dplyr::select(gene_name,is_candidate,
                                              CAGE,
                                              H3K27ac,
                                              H3K4me3_atTSS



)
colnames(tmp_data_ol)[5]="H3K4me3"


tmp_data_ol$is_candidate=ifelse(tmp_data_ol$is_candidate,"candidate","rest")

tmp_data_ol=tmp_data_ol%>% pivot_longer(cols = 3:5,names_to = "mark",values_to = "overlap")
ggplot(tmp_data_ol,aes(x=is_candidate,fill=overlap))+geom_bar(position = "fill") +
  facet_wrap(~ mark) +theme_classic()

tmp_data_ol=tmp_data_ol%>% group_by(is_candidate,mark)%>%
  summarise(fraction=sum(overlap)/n())

g=ggplot(tmp_data_ol,aes(x=is_candidate,y=fraction,fill=is_candidate))+geom_col() +
  facet_wrap(~ mark) +theme_classic() +
  scale_fill_manual("lncRNA",values = mycolors[c(2,6)]) +xlab("")+
  theme(strip.background = element_blank(),
        text = element_text(size=20))


g

save_tiff_svg(g,"tmp_plots/","enrichment_marks_candidates_rest",
              h=6,w=10)


tmp_data_ol=gene_info_lncRNAs%>%filter(exonic_type=="monoexonic")%>%dplyr::select(gene_name,is_candidate,
                                                                                  CAGE,
                                                                                  H3K27ac,
                                                                                  H3K4me3_atTSS



)
colnames(tmp_data_ol)[5]="H3K4me3"


tmp_data_ol$is_candidate=ifelse(tmp_data_ol$is_candidate,"candidate","rest")

tmp_data_ol=tmp_data_ol%>% pivot_longer(cols = 3:5,names_to = "mark",values_to = "overlap")
ggplot(tmp_data_ol,aes(x=is_candidate,fill=overlap))+geom_bar(position = "fill") +
  facet_wrap(~ mark) +theme_classic()

tmp_data_ol=tmp_data_ol%>% group_by(is_candidate,mark)%>%
  summarise(fraction=sum(overlap)/n())

g=ggplot(tmp_data_ol,aes(x=is_candidate,y=fraction,fill=is_candidate))+geom_col() +
  facet_wrap(~ mark) +theme_classic() +
  scale_fill_manual("lncRNA",values = mycolors[c(2,6)]) +xlab("")+
  theme(strip.background = element_blank(),
        text = element_text(size=20))

g


table(gene_info_lncRNAs_cand$EnhancerAtlas>0|gene_info_lncRNAs_cand$EnhancerFANTOM>0)


enhancer_list=list(Enhancer_Atlas=gene_info_lncRNAs_cand$gene_name[gene_info_lncRNAs_cand$EnhancerAtlas>0&!
                                                                     gene_info_lncRNAs_cand$EnhancerFANTOM==0],
                   Enhancer_FANTOM=gene_info_lncRNAs_cand$gene_name[gene_info_lncRNAs_cand$EnhancerAtlas==0&
                                                                      gene_info_lncRNAs_cand$EnhancerFANTOM>0],
                   both=gene_info_lncRNAs_cand$gene_name[gene_info_lncRNAs_cand$EnhancerAtlas>0&gene_info_lncRNAs_cand$EnhancerFANTOM>0])

enhancer_regions=sapply(enhancer_list,length)
pie(enhancer_regions,col=mycolors[c(2,6,1)],labels = paste0(names(enhancer_regions),
                                                            " (",enhancer_regions,")"))
save_tiff_svg(pie(enhancer_regions,col=mycolors[c(2,6,1)],labels = paste0(names(enhancer_regions),
                                                                          " (",enhancer_regions,")")),
              "tmp_plots/","pie_enhancers",h=5,w=7)


# TF binding ----
TF_binding=read.table("data/overlap_marks/Promoters_overlap_ALL_TFS_cistrome.txt")
TF_binding=TF_binding%>%filter(!is.na(TF_binding$V14))
TF_binding=TF_binding%>%filter(V14=="A")
TF_binding=TF_binding%>%filter(V4%in%gene_info_lncRNAs$gene_name)


table(TF_binding$V10[TF_binding$is_candidate])

TF_binding=TF_binding%>%group_by(V4,V10)%>%
  summarise(n())

TF_binding$is_candidate=gene_info_lncRNAs$is_candidate[match(TF_binding$V4,
                                                             gene_info_lncRNAs$gene_name)]

candidates=sort(table(TF_binding$V10[TF_binding$is_candidate])/sum(gene_info_lncRNAs$is_candidate))

rest=sort(table(TF_binding$V10[!TF_binding$is_candidate])/sum(!gene_info_lncRNAs$is_candidate))

rest=rest[match(names(candidates),names(rest))]
candidates=as.data.frame(candidates)
colnames(candidates)[2]="candidate"
rest=as.data.frame(rest)

TF_bind_fractions=left_join(candidates,rest)
colnames(TF_bind_fractions)=c("TF","candidates","rest")
sort(table(TF_binding$V10[!TF_binding$is_candidate]))
TF_bind_fractions$TF=rownames(TF_bind_fractions)

TF_bind_fractions$FC=TF_bind_fractions$candidates - TF_bind_fractions$rest




# expresion of selected candidates ----
# i'll manually select some

gene_info_lncRNAs <- left_join(gene_info_lncRNAs,
                               mean_corr_ext%>%dplyr::select(
                                 gene_name,mean_corr,log2FoldChange,padj,diff
                               ))

gene_info_lncRNAs_cand=gene_info_lncRNAs%>%filter(is_candidate)

View(gene_info_lncRNAs_cand%>%filter(!best_class%in%overlapping_classes) %>%
       dplyr::select(gene_name,mean_LSK_TPM,mean_corr,
                     exonic_type,annot,best_class,CAGE,closest_gene_PCG,EnhancerFANTOM,EnhancerAtlas))

monoexonic_candidates=gene_info_lncRNAs_cand$gene_name[gene_info_lncRNAs_cand$exonic_type=="monoexonic"&
                                                         !gene_info_lncRNAs_cand$best_class%in%overlapping_classes]

vst_monoexonic=vst[rownames(vst)%in%monoexonic_candidates,]

vst_monoexonic_long=pivot_longer_from_matrix(vst_monoexonic,val2 = "vst",nam2 = "sample")

vst_monoexonic_long$sample_type=coldata$cell_type[match(vst_monoexonic_long$sample,
                                                        coldata$sample)]
vst_monoexonic_long$sample=factor(vst_monoexonic_long$sample,
                                  levels=unique(vst_monoexonic_long$sample))

vst_monoexonic_long$sample_type=factor(vst_monoexonic_long$sample_type,levels = c("HSPC","progenitor","differentiated"))
g=ggplot(vst_monoexonic_long,
         aes(x = sample,y=vst,col=sample_type)) + geom_point() +
  facet_wrap(~gene_name) + theme_minimal() +
  scale_color_manual("cell type",values = mycolors)+theme(axis.text.x = element_blank(),
                                                          text = element_text(size = 25))+
  xlab("blood cells")

g

save_tiff_svg(g,"tmp_plots/",filename = "expression_signature_vst",
              h=6,w=10)

####
### final summary of all features ----
# features to plot

gene_info_lncRNAs_cand$Gname=prev_sel$Gname[match(gene_info_lncRNAs_cand$gene_name,
                                                  prev_sel$gene_name)]

gene_info_lncRNAs_cand$Gname[is.na(gene_info_lncRNAs_cand$Gname)]=
  gene_info_lncRNAs_cand$gene_name[is.na(gene_info_lncRNAs_cand$Gname)]

features2plot=gene_info_lncRNAs_cand%>%dplyr::select(gene_name,Gname,exonic_type,annot,
                                                     mean_corr,best_class)


vst_selected=vst[rownames(vst)%in%features2plot$gene_name,]

rownames(vst_selected)=features2plot$Gname[match(rownames(vst_selected),
                                                 features2plot$gene_name)]

vst_selected_long=pivot_longer_from_matrix(vst_selected,val2 = "vst",nam2 = "sample")


vst_selected_long$sample_type=coldata$cell_type[match(vst_selected_long$sample,
                                                      coldata$sample)]
vst_selected_long$sample=factor(vst_selected_long$sample,
                                levels=unique(vst_selected_long$sample))

vst_selected_long$sample_type=factor(vst_selected_long$sample_type,levels = c("HSPC","progenitor","differentiated"))


top10=(gene_info_lncRNAs_cand %>% arrange(-mean_corr)%>%filter(!best_class%in%overlapping_classes))$Gname[1:10]

g=ggplot(vst_selected_long%>%filter(gene_name%in%top10
),
aes(x = sample,y=vst,col=sample_type)) + geom_point() +
  facet_wrap(~gene_name) + theme_minimal() +
  scale_color_manual("cell type",values = mycolors)+theme(axis.text.x = element_blank(),
                                                          text = element_text(size = 25))+
  xlab("blood cells")

g


vst_selected_long=vst_selected_long%>%group_by(sample_type,gene_name)%>%
  summarise(mean_vst=mean(vst))


# mean vst per sample_type
vst_selected_wide=pivot_wider(vst_selected_long,
                              names_from = "sample_type",
                              values_from = "mean_vst")
vst_selected_wide=as.data.frame(vst_selected_wide)
rownames(vst_selected_wide)=vst_selected_wide$gene_name

vst_selected_wide=vst_selected_wide[,-1]



library(pheatmap)
pheatmap(vst_selected_wide,cluster_cols = F)
anno_row=features2plot
anno_row=as.data.frame(anno_row)
rownames(features2plot)=features2plot$Gname
anno_row=anno_row[,3:ncol(anno_row)]
pheatmap(vst_selected_wide,cluster_cols = F,annotation_row = anno_row)
anno_row$class=anno_row$best_class
anno_row$class[grepl("antisense",anno_row$class)]="antisense"
anno_row$class[grepl("sense",anno_row$class)]="sense"
anno_row=anno_row[,c(1:3,5)]
anno_row$annot=ifelse(anno_row$annot%in%c("lncRNA","TEC"),"annotated","not-annotated")
anno_row=anno_row%>%arrange(-mean_corr)

vst_selected_wide=vst_selected_wide[match(rownames(anno_row),
                                          rownames(vst_selected_wide)),]


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
pheatmap(as.matrix(vst_selected_wide),
         cluster_cols = F,cluster_rows = F,
         annotation_row = anno_row,fontsize = 20,show_rownames = F,  angle_col="0" )

tiff("tmp_plots/heatmap_multiple_features.tiff",
     res = 300, width = 10,height = 8,units = "in")
pheatmap(as.matrix(vst_selected_wide),
         cluster_cols = F,cluster_rows = F,
         annotation_row = anno_row[,1:3],fontsize = 20,show_rownames = F,  angle_col="0" )
dev.off()

final_feats$tr_level_support[is.na(final_feats$tr_level_support)]=3
final_feats[is.na(final_feats)]=0
final_feats$tr_level_support=final_feats$tr_level_support==3
final_feats$mean.counts_LSKStL_avobe80=final_feats$mean_counts_LSK_StL>80

final_feats=left_join(final_feats,vst_selected_wide,
                      by=c(gene_name="lncRNA"))



rownames(final_feats)=final_feats$gene_name

final_feats=final_feats%>%arrange(-mean_corr_signature)
final_feats_annot=final_feats%>% dplyr::select(exonic_type,
                                               tr_level_support,
                                               counts_perkb,
                                               mean_corr_signature,
                                               best_class,
                                               intramod_connect_PCG)
final_feats_expression=final_feats[,colnames(final_feats)%in%sample_types]
final_feats_annot$tr_level_support=ifelse(final_feats_annot$tr_level_support,"yes","no")
#final_feats_annot$mean.counts_LSKStL_avobe80=ifelse(final_feats_annot$mean.counts_LSKStL_avobe80,"yes","no")

library(pheatmap)

tiff("heatmap_multiple_features.tiff",res = 300,width = 10,height = 8,units = "in")
pheatmap(final_feats_expression,
         fontsize_col=6,
         annotation_row = final_feats_annot,
         fontsize_row=8,cluster_cols = F,
         cluster_rows = F,
         scale = "none",legend = F)
dev.off()

