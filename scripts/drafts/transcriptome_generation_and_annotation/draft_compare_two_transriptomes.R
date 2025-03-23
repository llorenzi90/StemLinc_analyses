# to dos: change color palette

# read inputs
#library(VennDiagram)
library(tidyverse)
tracking_path="outputs/gffcompare/SLvsKli_LSK.tracking"
tracking=read.table(tracking_path)

table(tracking$V5!="-",tracking$V6!="-")
x= list(
  tracking %>% filter(V5!="-") %>% pull(V1) ,
  tracking %>% filter(V6!="-") %>% pull(V1)
)

names(x)=c("StemLinc","Klimmeck")
library(ggvenn)
ggvenn(
  x,
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
) + ggtitle("test")

plot_pair_venn <- function(tra, filter_cond=rep(T,nrow(tracking)),sampnames=c("StemLinc","Klimmeck"),title=""){
  tra <- tra %>% filter(filter_cond)
  vennlist <- list(tra %>% filter(V5!="-") %>% pull(V1) ,
                   tra %>% filter(V6!="-") %>% pull(V1))
  names(vennlist)=sampnames
  g=ggvenn(
    vennlist,
    fill_color = c("#0073C2FF", "#EFC000FF"),
    stroke_size = 0.5, set_name_size = 4
  ) +ggtitle(title)
  print(g)
  return(list(vennlist,g))
}


plot_pair_venn(tracking,
               filter_cond = tracking$V4=="=",
               title = "Overlap of '=' transcripts")


# gene level summarization ----
# summarise at gene level using some trastools functions and new ones

# First do it with given genes
# Then do it with ref genes for overlapping transcripts

# things to add previously:
#     - ref transcript biotype
#     - Number of samples in which each transcript is assembled in each dataset (from each sample's tracking)
#     - Number of exons of each transcript (from tracking)
#     - expression of each transcript in each dataset (from each sample's tracking)

# I want to have: occurrence of each gene in each dataset and each sample.
# Note that this is a more "fair" comparison of assembled
# transcripts because we are assessing in what of the same regions we find transcription, no matter if it is not
# exactly matching in structure (is this because coverage is hardly enough to retrieve complete isoforms?)
# Anyways...

# What I need to do:

# A) DATA TRANSFORMATION

# First, add the needed info to the tracking file (I can modify the transcriptome char
# Rmd to retrieve annotated tracking, or even better make a specific script for this? although more time,
# keep it simple
# summarise at gene level using gene_id, but later do it using the reference gene id for the
# overlapping classcodes, the given by gffcompare for the rest

# output of that function:
#   dataframe with:
# gene_id, best classcode, biotype, gene class, exonic type1, exonic type2, max Nsamples1, max Nsamples2,

# B) PLOTS
# 1. The per-gene classcodes heatmap is nice to compare genes between both samples,
# make it nicer, add
# number/percentage of loci
order_cc=c("=","j","k","c","m","n","e","o","p","s","x","i","y","u","r",".")
ordered_classcodes=order_cc
get_cc_per_gene_per_sample <- function(tracking,ordered_classcodes=order_cc){
  c2 <- tracking %>% mutate(V5=V5!="-",V6=V6!="-") %>%
    pivot_longer(cols = c(V5,V6),names_to = "sample",values_to = "is") %>%
    filter(is) %>% group_by(V2,sample) %>% arrange(match(V4,ordered_classcodes)) %>%
    summarise(best_cc=V4[1])
  c2$sample=factor(c2$sample,levels = c("V5","V6"))
  c3 <- c2 %>% pivot_wider(names_from = sample,values_from = best_cc)
  c3[is.na(c3)] = "NA" # try to preserve the sample order
  c3 <- c3 %>% select(V2,V5,V6)
  return(c3)
  # make nicer heatmap
}

get_cc_per_gene_name_per_sample <- function(tracking,ordered_classcodes=order_cc){
  c2 <- tracking %>% mutate(StemLinc=V5!="-",Klimmeck=V6!="-") %>%
    pivot_longer(cols = c(StemLinc,Klimmeck),names_to = "sample",values_to = "is") %>%
    filter(is) %>% group_by(gene_name,sample) %>% arrange(match(V4,ordered_classcodes)) %>%
    summarise(best_cc=V4[1])
  c3 <- c2 %>% pivot_wider(names_from = sample,values_from = best_cc)
  c3[is.na(c3)] = "NA" # try to preserve the sample order
  c3 <- c3 %>% select(gene_name,StemLinc,Klimmeck)
  return(c3)
  # make nicer heatmap
}

cc_dat=get_cc_per_gene_name_per_sample(tracking)
cc_dat=get_cc_per_gene_per_sample(tracking )

library(ggplot2)

cc_dat=as.data.frame(table(cc_dat$V5,cc_dat$V6))

p=ggplot(cc_dat, aes(Var1, Var2, fill= Freq)) +
  geom_tile()
p
cc_dat <- cc_dat %>%
  mutate(text = paste0("StemLinc: ", Var1, "\n", "Kli: ", Var2, "\n", "Value: ",Freq, "\n", "Frac:", round(Freq/sum(Freq),2)))

library(plotly)
p=ggplot(cc_dat, aes(Var1, Var2, fill= Freq,text=text)) +
  geom_tile()
ggplotly(p, tooltip="text")

library(ggplot2)
library(plotly)
DF <- data.frame(A = c(1:10), B = c(1:10), C = (1:10))
myplots <- htmltools::tagList()

for (i in 1: dim(DF)[2]) {
  p = ggplot(DF, aes(x = DF[, i], y = C)) +
    geom_point()
  myplots[[i]] = plotly::as_widget(plotly::ggplotly(p))
}

myplots

# 2. still can try venn, including even the reference genes, this is nice to have an idea
# of what fraction of known genes are covered. Venn per biotype and novel classcodes.

#Summarise at gene level
# first assign gene name from ref, if present, XLOC if not
tracking <- tracking %>% mutate(gene_name = ifelse(V4%in%overlapping_class_codes,
                                       ifelse(!is.na(ref_gene_name_1),
                                              ref_gene_name_1,
                                              ref_gene_name_2),
                                       V2))

# trb = tracking %>% filter(!is.na(ref_gene_name_1)&!is.na(ref_gene_name_2))
# table(trb$ref_gene_name_1!=trb$ref_gene_name_2&trb$V4%in%overlapping_class_codes)
# table(trb$ref_biotype_1!=trb$ref_biotype_2&trb$V4%in%overlapping_class_codes)
# View(trb[trb$ref_biotype_1!=trb$ref_biotype_2&trb$V4%in%overlapping_class_codes,])

# group by gene_name summarise each sample occurrence
tracking_GL <- tracking %>% group_by(gene_name) %>% summarise(StemLinc=any(q1_gene!="-"),
                                                              Klimmeck=any(q2_gene!="-"))

# The following are split by feature bc doing max(feat,na.rm=T) took a lot of time,
# I realized arrange method is better, need to generalize
# mean tpm
tragl1 <- tracking %>% group_by(gene_name) %>% arrange(-mean_tpm_1) %>% summarise(maxTPM1=mean_tpm_1[1])
tragl2 <- tracking %>% group_by(gene_name) %>% arrange(-mean_tpm_2) %>% summarise(maxTPM2=mean_tpm_2[1])
GL1=tragl1
GL2=tragl2

# max Nsamps
tragl1 <- tracking %>% group_by(gene_name) %>% arrange(-Nsamps_1) %>% summarise(maxNsamps1=Nsamps_1[1])
tragl2 <- tracking %>% group_by(gene_name) %>% arrange(-Nsamps_2) %>% summarise(maxNsamps2=Nsamps_2[1])
GL1 <- left_join(GL1,tragl1)
GL2 <- left_join(GL2,tragl2)

# max Nexons
tragl1 <- tracking %>% group_by(gene_name) %>% arrange(-Nexons_1) %>% summarise(maxNexons1=Nexons_1[1])
tragl2 <- tracking %>% group_by(gene_name) %>% arrange(-Nexons_2) %>% summarise(maxNexons2=Nexons_2[1])
GL1 <- left_join(GL1,tragl1)
GL2 <- left_join(GL2,tragl2)

# pull everything into a single table
tracking_GL <- left_join(tracking_GL,GL1)
tracking_GL <- left_join(tracking_GL,GL2)

# best classcode per sample
c2 <- tracking %>% mutate(StemLinc=V5!="-",Klimmeck=V6!="-") %>%
  pivot_longer(cols = c(StemLinc,Klimmeck),names_to = "sample",values_to = "is") %>%
  filter(is) %>% group_by(gene_name,sample) %>% arrange(match(V4,ordered_classcodes)) %>%
  summarise(best_cc=V4[1])
c3 <- c2 %>% pivot_wider(names_from = sample,values_from = best_cc)
#c3[is.na(c3)] = "NA" # try to preserve the sample order
c3 <- c3 %>% select(gene_name,StemLinc,Klimmeck)
colnames(c3)[2:3]=paste0(colnames(c3)[2:3],"_cc")

tracking_GL <- left_join(tracking_GL,c3)
# gene_id, best classcode, biotype, gene class, exonic type1, exonic type2, max Nsamples1, max Nsamples2,

# finally add biotype of reference
tracking_GL$biotype <- ref_annot$simplified_gene_biotype[match(tracking_GL$gene_name,
                                                               ref_annot$gene_name)]
source("scripts/nejm_palette.R")

plot_pair_venn(tracking_GL,StemLinc,Klimmeck,gene_name,title = "Gene level overlap")
plot_pair_venn(tracking_GL,StemLinc,Klimmeck,gene_name,filter_cond = !is.na(tracking_GL$biotype),
               title = "Gene level overlap - annotated genes")

plot_pair_venn(tracking_GL,StemLinc,Klimmeck,gene_name,
               filter_cond = is.na(tracking_GL$biotype),
               title = "Gene level overlap - potential novel genes")

for (biot in c("protein_coding","lncRNA","pseudogene")) {
  p=plot_pair_venn(tracking_GL,
                   StemLinc,
                   Klimmeck,
                   gene_name,
                 filter_cond = !is.na(tracking_GL$biotype)&tracking_GL$biotype==biot,
                 title = paste0("Gene level overlap - ",biot))
  print(p)
}

plot_three_venn <- function(traGL,
                            samp1_var,
                            samp2_var,
                            samp3_var,
                            feat_var,
                            filter_cond=rep(T,nrow(traGL)), # defaults the entire data
                            sampnames=c("StemLinc","Klimmeck","Ref"),
                            title=""){
  traGL <- traGL %>% filter(filter_cond)
  if(is.logical(traGL%>%pull({{samp1_var}}))){
    vennlist <- list(traGL %>% filter({{samp1_var}}) %>% pull({{feat_var}}) ,
                     traGL %>% filter({{samp2_var}}) %>% pull({{feat_var}}),
                     traGL %>% filter({{samp3_var}}) %>% pull({{feat_var}}))
  }else{
    vennlist <- list(traGL %>% filter({{samp1_var}}!="-") %>% pull({{feat_var}}) ,
                     traGL %>% filter({{samp2_var}}!="-") %>% pull({{feat_var}}),
                     traGL %>% filter({{samp3_var}}!="-") %>% pull({{feat_var}}))
  }

  names(vennlist)=sampnames

  g=ggvenn(
    vennlist,
    fill_color = c("#0072B5FF","#FFDC91FF","#20854EFF"),
    stroke_size = 0.5, set_name_size = 4
  ) +ggtitle(title)
  print(g)
  #return(list(vennlist,g))
} # makes venn diagram of two sets of transcripts


#"#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF" "#FFDC91FF" "#EE4C97FF"
# 3. barplot of number of samples in each dataset
# per the top X most abundant groups in the heatmap

# 4. UpSetR plots for visualizing some things e.g.
# intronic genes vs number of samples,
# exonic type

# 5. Expression levels of genes in common and unique for each sample
tracking_GL <- tracking_GL %>% mutate(gene_class=ifelse(overlapRef,biotype,"potNovel"),
                                      in_sample=ifelse(StemLinc&Klimmeck,"both",ifelse(StemLinc,"StemLinc","Klimmeck")))
expression_dat <- tracking_GL%>%filter(gene_class%in%c("potNovel","protein_coding",
                                                       "lncRNA","TEC","pseudogene"))

expression_dat <- pivot_longer(expression_dat,cols = c("maxTPM1","maxTPM2"),
                               names_to = "sample", values_to = "TPM")

expression_dat <- expression_dat%>% mutate(sample = recode(sample,maxTPM1="StemLinc",
                                           maxTPM2="Klimmeck"))

expression_dat <- expression_dat%>%filter(!is.na(TPM))
ggplot(expression_dat, aes(x=in_sample,fill=sample,y=TPM)) +geom_boxplot()+
  scale_y_log10() + facet_wrap(~gene_class) + theme_bw()

ggplot(expression_dat, aes(x=in_sample)) +geom_bar()+
   facet_wrap(~gene_class) + theme_bw() + ylab("# genes")


expression_dat <- expression_dat %>% mutate(maxNsamps=ifelse(sample=="StemLinc",maxNsamps1,
                                                             maxNsamps2))

ggplot(expression_dat, aes(x=maxNsamps,fill=sample)) +geom_bar()+
  facet_wrap(~in_sample) + theme_bw() + ylab("# genes")

ggplot(expression_dat, aes(x=interaction(in_sample,maxNsamps),fill=sample,y=TPM)) +geom_boxplot()+
  scale_y_log10() + facet_wrap(~gene_class) + theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust = 1))
# 6. Expression level for known lncRNAs, PCGs and pseudogenes : in both samples, in one or the other,


# Questions: how many genes are in common and unique for each sample?




# Idea only retain intronic genes that are very highly expressed and/or are expressed in all samples


# get sample occurence per gene
sopg <- get_sample_occurrance_per_gene_from_tracking(tracking)

# get overlap ref
oRefGL <- get_overlapRef_gene_level(tracking = tracking)

# define get best classcode per gene
# the rule is simple, make hierarchy of classcodes, the classcode of the
# gene is defined as the one that is higher in the hierarchy

order_cc=c("=","j","k","c","m","n","e","o","p","s","x","i","y","u","r",".")

get_cc_per_gene <- function(tracking,ordered_classcodes=order_cc){
  c1 <- tracking %>% group_by(V2) %>%
    arrange(match(V4,ordered_classcodes)) %>%
    summarise(best_cc=V4[1])
  c2 <- tracking %>% mutate(V5=V5!="-",V6=V6!="-") %>%
    pivot_longer(cols = c(V5,V6),names_to = "sample",values_to = "is") %>%
    filter(is) %>% group_by(V2,sample) %>% arrange(match(V4,ordered_classcodes)) %>%
    summarise(best_cc=V4[1])
  c2$sample=factor(c2$sample,levels = c("V5","V6"))
  c3 <- c2 %>% pivot_wider(names_from = sample,values_from = best_cc)
  c3[is.na(c3)] = "NA" # try to preserve the sample order
  c3 <- c3 %>% select(V2,V5,V6)

  # make nicer heatmap
}


# Gene level summarizarion using reference genes ----
# first create column with new gene ids:
# for those with overlap with reference

# Set hierarchy of classcodes
# Then, per gene, order classcodes according to hierarchy and take the first one (see transcriptome report Rmd script)
#
# Decide on order (think quickly, document if needed)
# =, j, c, etc, for non-overlapping: p, x, i, u, r
# Make quick check of what ccs can coexist in non-overlapping genes
#
# Should I first separate in overlap and non-olap and then rename ??
#
#   Do it also for each sample separately. In this case, what happens if not in sample?
#
#   I can do the same but with transcripts overlapping ref
#
# What if a gene is '=' for one sample and 'o' for the other sample, how do I compare them?
#   With the approach i am using the gene will get "=" classcode. i could do an overlap/non-overlap and then a kind of heatmap of classcodes vs classcodes and matches
#
