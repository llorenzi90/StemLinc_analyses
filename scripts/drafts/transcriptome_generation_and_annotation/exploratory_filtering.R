# filtering criteria: fast, untidy need solution now check prev scripts
# select genes and classify in: potNovel, annotated, mono, multi, in Klimmeck
# and plot ATAC-seq
#https://rockefelleruniversity.github.io/Genomic_HeatmapsAndProfiles/presentations/slides/GenomicHeatmapsAndProfiles.html#44

source("scripts/compare_two_transcriptomes.R")

# exploratory analysis on expression filters:
# options: expression, number of samples in which the transcript or gene appears
# expression above 0.5 TPM

N_samps_per_gene1=get_n_samples_per_gene_from_tracking(tracking = tr1,cols = 6:8)

get_n_samples_per_gene_from_tracking_gene_name <- function(tracking,gene_id="gene_name",cols=6:8){
  tracking <- as.data.frame(tracking)
  tracking$V2=tracking[,gene_id]
  get_n_samples_per_gene_from_tracking(tracking,cols = cols)
}

tr1 <- tr1 %>% mutate(gene_name = ifelse(V4%in%overlapping_class_codes,
                                         ref_gene_name,V2))

N_samps_per_gene1=get_n_samples_per_gene_from_tracking_gene_name(tracking = tr1)
max_Nsamps_per_gene1=tr1 %>% group_by(gene_name) %>% summarise(maxNsamps=max(Nsamps))
table(max_Nsamps_per_gene1$gene_name==names(N_samps_per_gene1))
table(max_Nsamps_per_gene1$gene_name%in%names(N_samps_per_gene1))
max_Nsamps_per_gene1$Nsamps=N_samps_per_gene1[match(max_Nsamps_per_gene1$gene_name,
                                                    names(N_samps_per_gene1))]

table(max_Nsamps_per_gene1$maxNsamps,max_Nsamps_per_gene1$Nsamps)
max_Nsamps_per_gene1$biotype=ref_annot$simplified_gene_biotype[match(max_Nsamps_per_gene1$gene_name,
                                                                     ref_annot$gene_name)]
# assign cc to potNovel genes
ordered_classcodes=c("=",
                     "j",
                     "k",
                     "c",
                     "m",
                     "n",
                     "e",
                     "o",
                     "p",
                     "s",
                     "x",
                     "i",
                     "y",
                     "u",
                     "r",
                     ".")
best_cc_per_gene=tr1 %>% group_by(gene_name) %>% arrange(match(V4,ordered_classcodes)) %>%
  summarise(best_cc=V4[1])

max_Nsamps_per_gene1 <- left_join(max_Nsamps_per_gene1,best_cc_per_gene)
max_Nsamps_per_gene1 <- max_Nsamps_per_gene1 %>% mutate(gene_class=ifelse(is.na(biotype),best_cc,biotype))
max_Nsamps_per_gene1_biot <- max_Nsamps_per_gene1 %>%filter(is.na(biotype)|biotype%in%c("lncRNA",
                                                                                   "protein_coding",
                                                                                   "pseudogene",
                                                                                   "TEC"))
table(max_Nsamps_per_gene1_biot$gene_class,max_Nsamps_per_gene1_biot$Nsamps)
table(max_Nsamps_per_gene1_biot$gene_class,max_Nsamps_per_gene1_biot$maxNsamps)
table(max_Nsamps_per_gene1_biot$gene_class,max_Nsamps_per_gene1_biot$maxNsamps>1)

table(max_Nsamps_per_gene1_biot$gene_class,max_Nsamps_per_gene1_biot$Nsamps>2)

gene_level_info=max_Nsamps_per_gene1_biot

## add max TPM
max_TPM1=get_max_per_group_by_arrange(tr1,group_var = gene_name,arr_var = mean_tpm,outname = "max_mean_tpm")
gene_level_info <- left_join(gene_level_info,max_TPM1)

table(gene_level_info$gene_class,
      gene_level_info$Nsamps==3&gene_level_info$max_mean_tpm>=0.2)

table(gene_level_info$gene_class,
      gene_level_info$Nsamps==3|(gene_level_info$Nsamps==2&gene_level_info$max_mean_tpm>=1))

gene_level_info %>% group_by(gene_class) %>%summarise(median(max_mean_tpm))

ggplot(gene_level_info,aes(x=as.factor(Nsamps),y=max_mean_tpm)) +
  geom_boxplot() + facet_wrap(~gene_class) + scale_y_log10() +
  geom_hline(yintercept = 0.2)

table(gene_level_info$gene_class,
      gene_level_info$Nsamps==3&gene_level_info$max_mean_tpm>=0.2)

# for now I will use this criteria:
# retained genes were required to be assembled in 3 replicates
# and to have a an expression value >= 0.2 TPM


filter_transcriptome <- function(annotated_tracking,Nsamps_per_gene=3,maxTPM=0.2){
  annotated_tracking <- annotated_tracking %>% mutate(gene_name = ifelse(V4%in%overlapping_class_codes,
                                           ref_gene_name,V2))
  N_samps_per_gene=get_n_samples_per_gene_from_tracking_gene_name(annotated_tracking)
  N_samps_per_gene <- data.frame(gene_name=names(N_samps_per_gene),Nsamps=N_samps_per_gene)
  max_TPM <- get_max_per_group_by_arrange(annotated_tracking,
                                           group_var = gene_name,
                                           arr_var = mean_tpm,
                                           outname = "max_mean_tpm")
  gene_level_info <- left_join(N_samps_per_gene,max_TPM)
  best_cc_per_gene=annotated_tracking %>% group_by(gene_name) %>% arrange(match(V4,ordered_classcodes)) %>%
    summarise(best_cc=V4[1])

  gene_level_info <- left_join(gene_level_info,best_cc_per_gene)
  gene_level_info$biotype=ref_annot$simplified_gene_biotype[match(gene_level_info$gene_name,
                                                                       ref_annot$gene_name)]
  gene_level_info <- gene_level_info %>% mutate(gene_class=ifelse(is.na(biotype),best_cc,biotype))

  # filter
  genes2keep <- gene_level_info %>% filter(Nsamps==Nsamps_per_gene&max_mean_tpm>=maxTPM) %>% pull(gene_name)

  filtered_tracking <- annotated_tracking %>% filter(gene_name %in% genes2keep)
  return(filtered_tracking)
}

filt_tr1 <- filter_transcriptome(tr1)
filt_tr2 <- filter_transcriptome(tr2)

#annotate_tracking_with_2_datasetstr2annotate_tracking_with_2_datasets <- function(tracking, tr1, tr2, cols=c(5:ncol(tracking))){
  # things to add
  #     - ref transcript biotype
  #     - Number of samples in which each transcript is assembled in each dataset (from each sample's tracking)
  #     - Number of exons of each transcript (from tracking)
  #     - expression of each transcript in each dataset (from each sample's tracking)

trs=list(filt_tr1,filt_tr2)
trs=lapply(trs,function(t){
    t=t%>%dplyr::select(V1,ref_gene_name,ref_biotype,mean_tpm,Nexons,Nsamps)
    return(t)
  })
trs=lapply(seq_along(trs),function(i){
    t=trs[[i]]
    colnames(t)[2:ncol(t)]=paste0(colnames(t)[2:ncol(t)],"_",i)
    return(t)
  })
tracking_path="outputs/gffcompare/SLvsKli_LSK.tracking"
tracking=read.table(tracking_path)

tracking <- trastools::split_samples_info(tracking,  remove = F)
tracking <- left_join(tracking,trs[[1]],by=c(q1_transcript="V1"))
tracking <- left_join(tracking,trs[[2]],by=c(q2_transcript="V1"))

  #remove rows that are not in any of both transcriptomes
tracking <- tracking %>% filter(!(is.na(ref_gene_name_1)&is.na(ref_gene_name_2)))
  #return(tracking)
#}

tracking_GL <- summarize_GL(tracking)
tracking_GL$biotype <- ref_annot$simplified_gene_biotype[match(tracking_GL$gene_name,
                                                                 ref_annot$gene_name)]


# Compare both transcriptomes ----
# copied from "compare_two_transcriptomes.R"
# Transcript level comparison ----

plot_pair_venn(tracking,V5,V6,V1,

               title = "Transcript level exact match")

# per classcode
for (cc in c("=","c","e","i","j","k","m","n","o","p","u","x")) {
  plot_pair_venn(tracking,V5,V6,V1,
                 filter_cond = tracking$V4==cc,sampnames = sample_names,
                 title = paste0("Transcript level overlap of ",cc," transcripts"))
}

plot_pair_venn(tracking,V5,V6,V1,
               filter_cond = tracking$V4%in%overlapping_class_codes,
               title = "Transcripts that overlap Ref exact match")
plot_pair_venn(tracking,V5,V6,V1,
               filter_cond = !tracking$V4%in%overlapping_class_codes,
               title = "Transcripts that do not overlap Ref exact match")

# Gene level comparison ----


## Gene level plots ----

### Heatmaps classcodes ----
# myplots <- htmltools::tagList()
# plot_heatmap <- function(cc_dat,sampnames=sample_names,
#                          title=""){
#   samp1=sampnames[1]
#   samp2=sampnames[2]
#   cc_dat <- cc_dat %>%
#     mutate(text = paste0(samp1,": ", Var1, "\n", samp2,": ", Var2,
#                          "\n", "Value: ",Freq, "\n", "Frac:", round(Freq/sum(Freq),2)))
#   p=ggplot(cc_dat, aes(Var1, Var2, fill= Freq,text=text)) +
#     geom_tile() + ggtitle(title) + xlab(samp1) +
#     ylab(samp2)
#
#   # p=plotly::as_widget(plotly::ggplotly(p, tooltip="text"))
#   return(p)
# }



# cc_dat=as.data.frame(table(tracking_GL$StemLinc_cc,
#                            tracking_GL$Klimmeck_cc))
#
# p=plot_heatmap(cc_dat = cc_dat, title = "All assembled genes")
# myplots[[1]]=plotly::as_widget(plotly::ggplotly(p, tooltip="text"))
#
# i=2
# for (biot in unique(tracking_GL$biotype)) {
#   tra=tracking_GL%>%filter(biotype%in%biot)
#   cc_dat=as.data.frame(table(tra$StemLinc_cc,
#                              tra$Klimmeck_cc))
#   p=plot_heatmap(cc_dat,title = paste0("Assembled genes in ",biot,
#                                        " biotype"))
#   myplots[[i]] <- plotly::as_widget(plotly::ggplotly(p, tooltip="text"))
#   i=i+1
# }
# myplots
### Pairwise Venns ----
plot_pair_venn(tracking_GL,StemLinc,Klimmeck,gene_name,title = "Gene level overlap")
plot_pair_venn(tracking_GL,StemLinc,Klimmeck,gene_name,filter_cond = !is.na(tracking_GL$biotype),
               title = "Gene level overlap - annotated genes")

plot_pair_venn(tracking_GL,StemLinc,Klimmeck,gene_name,
               filter_cond = is.na(tracking_GL$biotype),
               title = "Gene level overlap - potential novel genes")

for (biot in c("protein_coding","lncRNA","pseudogene")) {
  plot_pair_venn(tracking_GL,
                 StemLinc,
                 Klimmeck,
                 gene_name,
                 filter_cond = !is.na(tracking_GL$biotype)&tracking_GL$biotype==biot,
                 title = paste0("Gene level overlap - ",biot))

}

### Venns including ref ----
tracking_GL_with_ref <- full_join(ref_annot %>% mutate(biotype=simplified_gene_biotype,Ref=TRUE)%>%dplyr::select(gene_name,biotype,Ref),
                                  tracking_GL %>% dplyr::select(gene_name,StemLinc,Klimmeck,biotype)
)

tracking_GL_with_ref[is.na(tracking_GL_with_ref)] <- FALSE

plot_three_venn(tracking_GL_with_ref,StemLinc,Klimmeck,Ref,gene_name,title = "Gene level overlap")

for (biot in unique(tracking_GL_with_ref$biotype)) {
  p=plot_three_venn(tracking_GL_with_ref,StemLinc,Klimmeck,Ref,gene_name,
                    title = paste0("Gene level overlap - ", biot," and potential novel"),
                    filter_cond = tracking_GL_with_ref$biotype%in%c(biot,"FALSE"))
  print(p)
}

### Barplot of number of samples in each dataset ----
# per the top X most abundant groups in the heatmap
cc_dat=as.data.frame(table(tracking_GL$StemLinc_cc,
                           tracking_GL$Klimmeck_cc))
cc_dat$combined_cc=paste(cc_dat$Var1,cc_dat$Var2,sep = ";")
X=15
top_classes=cc_dat %>% arrange(-Freq) %>%slice(1:X)
top_classes
tracking_GL$combined_cc=paste(tracking_GL$StemLinc_cc,
                              tracking_GL$Klimmeck_cc,sep = ";")

tra <- tracking_GL %>% filter(combined_cc%in%top_classes$combined_cc)
tra <- tra%>%dplyr::select(combined_cc, max_Nsamps_1,max_Nsamps_2) %>%
  pivot_longer(cols = 2:3,names_to = "sample",values_to = "maxNsamps") %>%
  mutate(sample = recode(sample, "max_Nsamps_1" = "StemLinc", "max_Nsamps_2" = "Klimmeck"))

tra <- tra[!is.na(tra$maxNsamps),]

ggplot(tra,aes(x=sample,fill=as.factor(maxNsamps))) + geom_bar() + facet_wrap(~combined_cc)
ggplot(tra,aes(x=sample,fill=as.factor(maxNsamps))) + geom_bar(position = "fill") + facet_wrap(~combined_cc)


### UpSet plots ----
#### Intronic genes and number of samples ----
# intronic genes vs number of samples,

intronic_genes <- tracking_GL %>% filter(combined_cc%in%c("i;i","i;NA","NA;i"))
StemLinc_intronic=lapply(1:3,function(i)intronic_genes%>%filter(max_Nsamps_1==i) %>% pull(gene_name))
Klimmeck_intronic=lapply(1:3,function(i)intronic_genes%>%filter(max_Nsamps_2==i) %>% pull(gene_name))
intronic_genes_list <- c(StemLinc_intronic,Klimmeck_intronic)
names(intronic_genes_list)=c(paste0("SL_intronic_maxSamps=",1:3),
                             paste0("Kli_intronic_maxSamps=",1:3))

upset(fromList(intronic_genes_list), order.by = "freq",nsets = 6)
grid.text("Intronic genes and max number of assembled samples",x = 0.65, y=0.95, gp=gpar(fontsize=14))

#### Exonic type ----
exonic_type_list <- list(StemLinc_mono=tracking_GL%>%filter(max_Nexons_1==1) %>% pull(gene_name),
                         StemLinc_multi=tracking_GL%>%filter(max_Nexons_1>1) %>% pull(gene_name),
                         Klimmeck_mono=tracking_GL%>%filter(max_Nexons_2==1) %>% pull(gene_name),
                         Klimmeck_multi=tracking_GL%>%filter(max_Nexons_2>1) %>% pull(gene_name))
sets <- fromList(exonic_type_list)
sets$gene_name=unique(unlist(exonic_type_list,use.names = F))
tracking_GL$overlapRef=!is.na(tracking_GL$biotype)
sets$overlapRef=tracking_GL$overlapRef[match(sets$gene_name,
                                             tracking_GL$gene_name)]
upset(sets,
      query.legend = "bottom", nsets = 6, number.angles = 30, point.size = 3.5, line.size = 2,
      queries = list(

        list(
          query = elements,
          params = list("overlapRef",T),
          color =  "#20854EFF",
          active = T,
          query.name = "Overlaps Reference"
        )
      )
)
grid.text("Exonic type and overlap Ref",x = 0.65, y=0.95, gp=gpar(fontsize=14))

upset(sets %>% filter(!overlapRef),
      query.legend = "bottom", nsets = 6, number.angles = 30, point.size = 3.5, line.size = 2,

)
grid.text("Exonic type of potential novel genes",x = 0.65, y=0.95, gp=gpar(fontsize=14))

common_new_genes=tracking_GL$gene_name[tracking_GL$StemLinc&tracking_GL$Klimmeck&!tracking_GL$overlapRef]
common_new_genes_in_set=sets$gene_name[!sets$overlapRef&rowSums(sets[,1:4])>1]
missing_genes_in_set=common_new_genes[!common_new_genes%in%common_new_genes_in_set]

## Expression plots ----
tracking_GL <- tracking_GL %>% mutate(gene_class=ifelse(overlapRef,biotype,"potNovel"),
                                      in_sample=ifelse(StemLinc&Klimmeck,"both",ifelse(StemLinc,"StemLinc","Klimmeck")))
expression_dat <- tracking_GL%>%filter(gene_class%in%c("potNovel","protein_coding",
                                                       "lncRNA","TEC","pseudogene"))

expression_dat <- pivot_longer(expression_dat,cols = c("max_mean_tpm_1","max_mean_tpm_2"),
                               names_to = "sample", values_to = "TPM")

expression_dat <- expression_dat%>% mutate(sample = recode(sample,max_mean_tpm_1="StemLinc",
                                                           max_mean_tpm_2="Klimmeck"))

expression_dat <- expression_dat%>%filter(!is.na(TPM))
ggplot(expression_dat, aes(x=in_sample,fill=sample,y=TPM)) +geom_boxplot()+
  scale_y_log10() + facet_wrap(~gene_class) + theme_bw()

ggplot(expression_dat%>%filter(!duplicated(gene_name)), aes(x=in_sample)) +geom_bar()+
  facet_wrap(~gene_class) + theme_bw() + ylab("# genes")

# add max number of samples
expression_dat <- expression_dat %>% mutate(maxNsamps=ifelse(sample=="StemLinc",max_Nsamps_1,
                                                             max_Nsamps_2))

ggplot(expression_dat, aes(x=interaction(in_sample,maxNsamps))) +geom_bar()+
  facet_wrap(~gene_class) + theme_bw() + ylab("# genes") + theme(axis.text.x = element_text(angle = 45,hjust = 1))


ggplot(expression_dat, aes(x=interaction(in_sample,maxNsamps),fill=sample,y=TPM)) +geom_boxplot()+
  scale_y_log10() + facet_wrap(~gene_class) + theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust = 1))


