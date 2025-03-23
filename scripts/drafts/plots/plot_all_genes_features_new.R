## Head -------------------------------------
##
##
## Purpose of script:
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-05-10
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
library(rtracklayer)
library(ggpubr)
## Load data---------------------------
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/")
source("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/nejm_palette.R")
counts=read.table("data/expression_files/featureCounts/gene_expression_PCGs_lncRNAs_pseudo_unfiltered.240509.counts.txt")
TPMs=read.table("data/expression_files/featureCounts/gene_expression_PCGs_lncRNAs_pseudo_unfiltered.240509.TPM.txt")
counts_annot=read.table("data/expression_files/featureCounts/gene_expression_PCGs_lncRNAs_pseudo_unfiltered.240509.annot.txt",header = T)
gene_info=read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/data/gene_info_PCGs_lncRNAs_pseudo_unfiltered.240509.csv")
GTF=readGFF("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/data/various_gtf_files/merged_PCG_lncRNA_pseudo_transcriptome.gtf")
GTF%>%group_by(strand)%>%summarise(n=length(unique(gene_name)))

GTF$strand[GTF$strand=="*"] <- genes_strands[match(GTF$gene_name[GTF$strand=="*"],
                                                   names(genes_strands))]
any(GTF$strand=="*")
export(GTF,"test.gtf",format = "gtf")
# mv test.gtf merged_PCG_lncRNA_pseudo_transcriptome.gtf

genes_strands=c("XLOC_007631"="-",
                "XLOC_068623"="+",
                "XLOC_068641"="+",
                "XLOC_068650"="+",
                "XLOC_095109"="+",
                "XLOC_106083"="+",
                "XLOC_047358"="+",
                "XLOC_047370"="-",
                "XLOC_053951"="+",
                "XLOC_053960"="+"
)

gene_info$strand[match(names(genes_strands),gene_info$gene_name)] <- genes_strands


save_tiff_svg=function(plo,outdir=".",filename="plot",
                       h=10,w=10){
  tiff(res = 300, height = h ,width = w, units = "in",
      paste0(outdir,"/",filename,".tiff"))
  print(plo)
  dev.off()
  
  svg( height = h ,width = w,
       paste0(outdir,"/",filename,".svg"))
  print(plo)
  dev.off()
  
}
source("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/scripts/StemLinc_code/make_longer_make_summ_expression_data.R")

## change biotype by simplified biotype ----
gene_info$biotype_ext=gene_info$biotype
gene_info$biotype=gene_info$simpl_biotype
gene_info$simpl_biotype=gene_info$biotype_ext
colnames(gene_info)[colnames(gene_info)=="simpl_biotype"]="detailed_biotype"
gene_info=gene_info[,-ncol(gene_info)]

wgi <- function(gi) {
  write.csv(gene_info,"gene_info_PCGs_lncRNAs_pseudo_unfiltered.240511.csv")}

wgi(gene_info)

mainbiots=c("lncRNA","PotNovel","PCG")


# setwd ----
setwd("analyses/plots_May/plots_various_features_all_biotypes/")
# set levels of biotypes ----
gene_info$biotype=factor(gene_info$biotype,
                         levels = c("lncRNA","PotNovel","PCG","TEC","pseudogene"))

mycolors=nejm_pal[c(1,3,4,5,2,6,7,8)]
# expression plots ----
StemLinc_samps=paste0("LSK_StL_",1:3)

# LSK_StL_counts=counts[,colnames(counts)%in%StemLinc_samps]
# 
# cutoff=c()
# table(rowSums(LSK_StL_counts>=10)==3)
# table(rowSums(LSK_StL_counts>=5)==3)
# 
# mean_LSK_StL_counts=apply(LSK_StL_counts,1,mean)
# table(mean_LSK_StL_counts>=10)
# table(mean_LSK_StL_counts>=5)
# 
LSK_StL_TPMs=TPMs[,colnames(TPMs)%in%StemLinc_samps]

cutoff=c()
table(rowSums(LSK_StL_TPMs>=0.3)==3)
table(rowSums(LSK_StL_TPMs>=5)==3)

mean_LSK_StL_TPMs=apply(LSK_StL_TPMs,1,mean)
table(mean_LSK_StL_TPMs>=0.3)
table(mean_LSK_StL_TPMs>=0.1)

# # make longer
# LSK_StL_TPMs_long=make_longer(LSK_StL_TPMs)
# 
# colnames(LSK_StL_TPMs_long)[1]="gene_name"
# LSK_StL_TPMs_long=left_join(LSK_StL_TPMs_long,counts_annot%>%
#                               dplyr::select(gene_name,simpl_biotype))


setwd("analyses/plots_May/plots_various_features_all_biotypes/")

### boxplot ----

g=ggplot(gene_info,aes(x=biotype,y=mean_LSK_TPM)) + 
  geom_boxplot()+
  theme_classic() +scale_y_log10() +ggtitle("LSK StemLinc mean TPMs")
g


save_tiff_svg(g,filename = "Boxplot_mean_TPM_LSK_StL")

#### passfilter ----
g=ggplot(gene_info%>%filter(passfilter),aes(x=simpl_biotype,y=mean_LSK_TPM)) + 
  geom_boxplot()+
  theme_classic() +scale_y_log10() +ggtitle("LSK StemLinc mean TPMs - passfilter genes")
g


save_tiff_svg(g,filename = "Boxplot_mean_TPM_LSK_StL_filtered")

g=ggplot(gene_info,aes(x=simpl_biotype,y=mean_LSK_TPM,col=passfilter)) + 
  geom_boxplot() +
  theme_classic() +scale_y_log10() +ggtitle("LSK StemLinc mean TPMs")

g

save_tiff_svg(g,filename = "Boxplot_mean_TPM_LSK_StL_pass-nopass")


### violin plot ----
g=ggplot(gene_info,aes(x=biotype,y=mean_LSK_TPM,
                       fill=biotype)) + 
  geom_violin()+
  scale_fill_manual(values = mycolors)+
  geom_hline(yintercept = 0.5,linetype="dotted", 
             color = "grey", size=1.5)+
  theme_classic() +
  scale_y_log10() +ggtitle("")+ylab("Mean TPM")+
  theme(text = element_text(size = 20),
        legend.position = "none")+xlab("")
g

save_tiff_svg(g,h = 8,w = 10,
             
              filename = "Violin_mean_TPM_LSK_StL_with_cutoff0.5")


g=ggplot(gene_info%>%filter(mean_LSK_TPM>=0.5),aes(x=biotype,y=mean_LSK_TPM,
                       fill=biotype)) + 
  geom_violin()+
  scale_fill_manual(values = mycolors)+
  theme_classic() +scale_y_log10() +ggtitle("LSK StemLinc mean TPMs")
g

save_tiff_svg(g,filename = "Violin_mean_TPM_LSK_StL")


#### passfilter ----
g=ggplot(gene_info%>%filter(passfilter),aes(x=simpl_biotype,y=mean_LSK_TPM)) + 
  geom_violin()+
  theme_classic() +scale_y_log10() +ggtitle("LSK StemLinc mean TPMs")
g

save_tiff_svg(g,filename = "Violin_mean_TPM_LSK_StL_filtered")

### density plot ----

g=ggplot(gene_info,aes(x=mean_LSK_TPM,col=biotype)) + 
  geom_density()+
  theme_classic() +ggtitle("LSK StemLinc mean TPMs") +scale_x_log10()

g

save_tiff_svg(g,filename = "density_mean_TPM_LSK_StL")
#### passfilter ----
g=ggplot(gene_info%>%filter(passfilter),
         aes(x=mean_LSK_TPM,col=biotype)) + 
  geom_density()+
  theme_classic() +ggtitle("LSK StemLinc mean TPMs") +scale_x_log10()
g
save_tiff_svg(g,filename = "density_mean_TPM_LSK_StL_passfilter")

g=ggplot(gene_info%>%filter(passfilter),
         aes(x=mean_LSK_TPM,col=biotype)) + 
  stat_ecdf()+
  theme_classic() +ggtitle("LSK StemLinc mean TPMs") +scale_x_log10()
g

g=ggplot(gene_info%>%filter(passfilter),
         aes(x=biotype,
             y=mean_LSK_TPM)) + 
  geom_boxplot()+
  theme_classic() +
  ggtitle("LSK StemLinc mean TPMs") +scale_y_log10()
g

## filter by expression ----
gene_info_filt=gene_info%>%filter(mean_LSK_TPM>=0.5)

## N exons N genes and exonic type ----

N_exons=GTF%>% filter(type=="exon") %>% group_by(transcript_id) %>% 
  summarise(N_exons=n())
N_exons$gene_name=GTF$gene_name[match(N_exons$transcript_id,
                                  GTF$transcript_id)]
gene_level_info=N_exons%>% group_by(gene_name) %>% 
  summarise(N_exons_tr=paste(N_exons,collapse = ","),
            max_exons=max(N_exons),
            exonic_type=ifelse(any(N_exons!="1"),"multiexonic","monoexonic"),
            N_transcripts=n())
source("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/scripts/StemLinc_code/compute_gene_non_redundant_exonic_length.function.R")

exonic_length <- compute_gene_non_redundant_exonic_length(in_gtf = GTF)

gene_level_info$exonic_length <- exonic_length$exonic_length[match(gene_level_info$gene_name,
                                                           exonic_length$gene_id)]


gene_info <- left_join(gene_info,gene_level_info)

wgi(gene_info)

### fractions plot ----
gene_info_filt=gene_info%>%filter(mean_LSK_TPM>=0.5)

tmp_ginfo=gene_info_filt
tmp_ginfo$biot=tmp_ginfo$biotype
tmp_ginfo$biotype=as.character(tmp_ginfo$biotype)
tmp_ginfo_biots=unique(tmp_ginfo$biotype)
tmp_ginfo_biots

genes_per_biot=tmp_ginfo%>%group_by(biotype)%>%
  summarise(N_genes=n())
genes_per_biot$N_genes=as.character(genes_per_biot$N_genes)
genes_per_biot$N_genes=c("13.303",
                         "1.184",
                         "1.221",
                         "3.675",
                         "208")
tmp_ginfo_biots=paste0(genes_per_biot$biotype," (",
                       genes_per_biot$N_genes,")")
tmp_ginfo$biotype=tmp_ginfo_biots[match(tmp_ginfo$biotype,
                                        genes_per_biot$biotype)]

tmp_ginfo$biotype=factor(tmp_ginfo$biotype,
                         levels = unique(tmp_ginfo$biotype)[c(5,3,2,1,4)])
g=ggplot(tmp_ginfo,
         aes(x=biotype,fill=exonic_type)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = mycolors[c(6,5)])+

  theme_classic() +
  ggtitle("")+ylab("Fraction")+
  theme(text = element_text(size = 20))+xlab("") +coord_flip()
g
save_tiff_svg(g,
              filename = "Ngenes_exonic_type.all_biotypes")


g=ggplot(tmp_ginfo%>%filter(biot%in%mainbiots),
         aes(x=biotype,fill=exonic_type)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = mycolors[c(6,5)])+
  
  theme_classic() +
  ggtitle("")+ylab("Fraction")+
  theme(text = element_text(size = 20))+xlab("") +coord_flip()
g
save_tiff_svg(g,
              filename = "Ngenes_exonic_type.main_biotypes")

## length distributioin ----
gene_info_filt$biotype=factor(gene_info_filt$biotype,
                              levels = c("lncRNA","TEC",
                                         "PotNovel","PCG"))

median_Exonic_length=gene_info_filt%>%
  filter(biotype%in%mainbiots)%>% group_by(biotype)%>%
  summarise(median_Exonic_length=median(exonic_length/1000),
            ypos=median(exonic_length/1000))
median_Exonic_length$lab=c("1.4 Kb","2.5 Kb","4.2 Kb")
median_Exonic_length=gene_info_filt%>%group_by(biotype)%>%
  summarise(median_Exonic_length=median(exonic_length/1000),
            ypos=median(exonic_length/1000))
median_Exonic_length$lab=c("1.4 Kb","2.0 Kb","2.6 Kb","4.2 Kb")

table(gene_info_filt$biotype)

g=ggplot(gene_info_filt,aes(x=exonic_length,col=biotype))+ 
 stat_ecdf() +
  theme_classic() + 
  ggtitle("Gene length")+scale_x_log10()

g
g=ggplot(gene_info_filt,aes(x=exonic_length,col=biotype))+ 
  geom_density() +
  theme_classic() + 
  ggtitle("Gene length")+scale_x_log10()

g

g=ggplot(gene_info_filt,
         aes(x=biotype,fill=biotype,y=exonic_length/1000))+ 
  geom_violin() +
  theme_classic() + geom_text(data = median_Exonic_length,
                              aes(label=lab,y=ypos))+
  ggtitle("Gene exonic length")+scale_y_log10()+
  ylab("Exonic length (Kb)")+
  xlab("")+
scale_fill_manual(values = mycolors[c(1,4,2,3)])+
  theme(text=element_text(size=20))
g

save_tiff_svg(g,outdir = "tmp_plots/",
              filename = "exonic_length_violin_nopseudo")
g=ggplot(gene_info_filt%>%filter(biotype%in%mainbiots),
         aes(x=biotype,fill=biotype,y=exonic_length))+ 
  geom_boxplot() +
  theme_classic() + 
  ggtitle("Gene length")+scale_y_log10()

g

geom_text(data = labeldat, aes(label = labels, y = ypos), 
          position = position_dodge(width = .75), 
          show.legend = FALSE )
## coding potential ----
bestORF_matched_strand=read.table("analyses/CPAT/bestORF_gene_strand_all_genes.gene_level.may24.txt",header = T)
bestORF_matched_strand=bestORF_matched_strand[!duplicated(bestORF_matched_strand$gene_name),]
gene_info <- left_join(gene_info , bestORF_matched_strand%>% 
                         select(gene_name,Coding_prob))



g=ggplot(gene_info,aes(x=Coding_prob,col=biotype))+ 
  geom_histogram() +
  theme_classic() + 
  geom_vline(xintercept = 0.44, linetype="dotted", 
             color = "grey", size=1.5) + 
  ggtitle("Coding probability distribution - best ORF - gene level")

g

save_tiff_svg(plo = g,
              
              filename = "Coding_prob_histogram_all_biotypes")

g=ggplot(gene_info,aes(x=Coding_prob,col=biotype))+ 
  geom_density(size=2) +
  theme_classic() + 
  scale_color_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
  geom_vline(xintercept = 0.44, linetype="dotted", 
             color = "grey", size=1.5) + 
  ggtitle("Coding probability distribution - best ORF - gene level")

g

save_tiff_svg(plo = g,
              
              filename = "Coding_prob_density_all_biotypes")

### 3 main biots only ----
mainbiots=c("lncRNA","PotNovel","PCG")

g=ggplot(gene_info%>%filter(biotype%in%mainbiots),aes(x=Coding_prob,col=biotype))+ 
  geom_density(size=2) +
  theme_classic() + 
  scale_color_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
  geom_vline(xintercept = 0.44, linetype="dotted", 
             color = "grey", size=1.5) + 
  ggtitle("Coding probability distribution - best ORF - gene level")

g

save_tiff_svg(plo = g,
              
              filename = "Coding_prob_density_Main_biotypes_unfiltered")

### passfilter ----
g=ggplot(gene_info%>%filter(biotype%in%mainbiots,
                            mean_LSK_TPM>=0.5),
         aes(x=Coding_prob,col=biotype))+ 
  geom_density(size=2) +
  theme_classic() + 
  scale_color_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
  geom_vline(xintercept = 0.44, linetype="dotted", 
             color = "grey", size=1.5) + 
  ggtitle("Coding probability distribution - best ORF - gene level")

g

save_tiff_svg(plo = g,
              
              filename = "Coding_prob_density_Main_biotypes_filtered")

## conservation ----
gene_info=read.csv("gene_info_PCGs_lncRNAs_pseudo_unfiltered.240511.csv")
gene_info_unfiltered=gene_info
gene_info=gene_info%>%filter(passfilter)

# phastCons35wayPotNovelGencode <- read.table("Experimento4_LSK_190224/analyses/conservation/phastCons35way/potNovel_3reps.gencodelncRNAs.phastCons35way.txt",header = T)
# phastCons35wayRefSeq <- read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/APRIL/phastCons35way/RefSeq_finalgenes2add.phastCons35way.txt",header = T)
# phastCons35wayPCG <- read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/data/PCGs.phastCons35way.240510.txt",header = T)
# 
# View(gene_info[!gene_info$gene_name%in%c(phastCons35wayPCG$group_name,phastCons35wayPotNovelGencode$gene_name,phastCons35wayRefSeq$gene_name),])
# phastCons35wayPotNovelGencode$gene_name=
#   gene_info$gene_name[match(phastCons35wayPotNovelGencode$gene_id,
#                             gene_info$id_countTable)]
# 
# View(gene_info[!gene_info$gene_name%in%c(phastCons35wayPCG$group_name,phastCons35wayPotNovelGencode$gene_name,phastCons35wayRefSeq$gene_name),])
# 
# colnames(phastCons35wayPCG)[1]="gene_name"
# colnames(phastCons35wayRefSeq)[2]="mean_phastCons35way"
# colnames(phastCons35wayPotNovelGencode)[2]="mean_phastCons35way"
# phastCons_all=rbind(phastCons35wayPCG%>%dplyr::select(gene_name,mean_phastCons35way),
#                     phastCons35wayPotNovelGencode%>%dplyr::select(gene_name,mean_phastCons35way),
#                     phastCons35wayRefSeq%>%dplyr::select(gene_name,mean_phastCons35way))
# 
# 
# table(gene_info$gene_name%in%phastCons_all$gene_name)
# gene_names_to_match=gene_info$gene_name
# gene_names_to_match[!gene_names_to_match%in%phastCons_all$gene_name] <- 
#   gene_info$id_countTable[!gene_names_to_match%in%phastCons_all$gene_name]
# 
# table(gene_names_to_match%in%phastCons_all$gene_name)
# 
# gene_info$mean_phastCons35way=phastCons_all$mean_phastCons35way[
#   match(gene_names_to_match,phastCons_all$gene_name)
# ]
# 
# wgi(gene_info)

## add conservation of promoter regions ----
prom_phastCons=read.table("data/phastCons35way.Promoters.240512.txt",header = T)
colnames(prom_phastCons)[2]="promoter_phastCons35way"
gene_info=left_join(gene_info,prom_phastCons)

### plot conservation ----

g=ggplot(gene_info%>%filter(biotype%in%mainbiots,
                            mean_LSK_TPM>=0.5),
         aes(x=mean_phastCons35way,col=biotype))+ 
  geom_density(size=2) +
  theme_classic() + 
  scale_color_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
 
  ggtitle("")

g

save_tiff_svg(g,
              filename = "Conservation_main_biotypes.density")

g=ggplot(gene_info%>%filter(biotype%in%mainbiots,
                            mean_LSK_TPM>=0.5),
         aes(x=promoter_phastCons35way,col=biotype))+ 
  geom_density(size=2) +
  theme_classic() + 
  scale_color_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
  
  ggtitle("")


g

cons=gene_info%>%filter(passfilter)%>%dplyr::select(gene_name,biotype,mean_phastCons35way,promoter_phastCons35way)
colnames(cons)[3:4]=c("gene_body","promoter")
cons=pivot_longer(cons,cols = c(3,4),names_to = "conservation_at",
                  values_to = "mean.phastCons35way")
cons

cons$biotype=factor(cons$biotype,levels=c("lncRNA",
                                          "PotNovel",
                                          "PCG"))
# ecdf
g=ggplot(gene_info%>%filter(biotype%in%mainbiots,
                            mean_LSK_TPM>=0.5),
         aes(x=mean_phastCons35way,col=biotype))+ 
  stat_ecdf(size=2) +
  theme_classic() + 
  scale_color_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
  ylab("Cumulative fraction")+
  ggtitle("")

g

save_tiff_svg(g,
              filename = "Conservation_main_biotypes.ecdf")

g=ggplot(cons%>%filter(biotype%in%mainbiots),
         aes(x=mean.phastCons35way,col=biotype,linetype=conservation_at))+ 
  stat_ecdf(size=2) +
  theme_classic() + 
  scale_color_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
  ylab("Cumulative fraction")+
  ggtitle("")

g

setwd("analyses/plots_May/plots_various_features_all_biotypes/")
save_tiff_svg(g,
              filename = "Conservation_main_biotypes.promoter_and_gene_body.ecdf")

g=ggplot(cons%>%filter(biotype%in%mainbiots),
         aes(y=mean.phastCons35way,x=biotype,fill=conservation_at))+ 
  geom_boxplot() +
  theme_classic() + 
  scale_color_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
  ylab("PhastCons35way")+
  ggtitle("")

g

save_tiff_svg(g,
              filename = "Conservation_main_biotypes.promoter_and_gene_body.boxplot")

# density
g=ggplot(cons%>%filter(biotype%in%mainbiots),
         aes(x=mean.phastCons35way,col=biotype,linetype=conservation_at))+ 
  geom_density(size=2) +
  theme_classic() + 
  scale_color_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
  ylab("Cumulative fraction")+
  ggtitle("")

g

setwd("analyses/plots_May/plots_various_features_all_biotypes/")
save_tiff_svg(g,
              filename = "Conservation_main_biotypes.promoter_and_gene_body.density")

## expression ecdf ----

g=ggplot(gene_info%>%filter(biotype%in%mainbiots,
                            mean_LSK_TPM>=0.5),
         aes(x=mean_LSK_TPM,col=biotype))+ 
  stat_ecdf(size=2) +
  theme_classic() + 
  scale_color_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
  ylab("Cumulative fraction")+ scale_x_log10()+
  ggtitle("") + xlab("mean TPM")

g

save_tiff_svg(g,
              filename = "Expression_main_biotypes.ecdf")

gene_info$exonic_type=factor(gene_info$exonic_type,
                             levels = c("multiexonic",
                                        "monoexonic"))
g=ggplot(gene_info%>%filter(biotype%in%mainbiots,
                            mean_LSK_TPM>=0.5),
         aes(x=mean_LSK_TPM,col=biotype,linetype=exonic_type))+ 
  stat_ecdf(size=2) +
  theme_classic() + 
  scale_color_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
  ylab("Cumulative fraction")+ scale_x_log10()+
  ggtitle("") + xlab("mean TPM")

g

save_tiff_svg(g,
              filename = "Expression_main_biotypes.ecdf.exonic_type")

g=ggplot(gene_info%>%filter(biotype%in%mainbiots,
                            mean_LSK_TPM>=0.5),
         aes(x=interaction(exonic_type,biotype),y=mean_LSK_TPM,
             fill=biotype))+ 
 geom_boxplot() +
  theme_classic() + 
  scale_fill_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
  ylab("")+ scale_y_log10()+
  ggtitle("") + ylab("mean TPM")

g

save_tiff_svg(g,
              filename = "Expression_main_biotypes.boxplot.exonic_type")

# class distribution ----
lncRNA_class=read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/data/gene_classif/ncRNAs_classification.240512.txt",header = T)


gene_info=gene_info%>%filter(passfilter)
gene_info_ncRNAs=gene_info%>%filter(biotype!="PCG")
gene_info_ncRNAs_unfilt=gene_info_unfiltered%>%filter(biotype!="PCG")
colnames(lncRNA_class)[c(1:3,5)]=c("class","closest_gene","all_class","gene_name")
gene_info_ncRNAs <- left_join(gene_info_ncRNAs,
                              lncRNA_class%>%dplyr::select(gene_name,class,closest_gene))

gene_info_ncRNAs_unfilt <- left_join(gene_info_ncRNAs_unfilt,
                              lncRNA_class%>%dplyr::select(gene_name,class,closest_gene))

table(gene_info_ncRNAs_unfilt$passfilter,gene_info_ncRNAs_unfilt$clas)

table(gene_info$gene_name%in%old_gene_info$gene_name)

ggplot(gene_info_ncRNAs,aes(x=exonic_type,fill=class)) +
  geom_bar(position = "fill") +theme_classic() +facet_wrap(~ biotype) 

gene_info_ncRNAs$simpl_class=gene_info_ncRNAs$class
classes=sort(unique(lncRNA_class$class))

simpl_classes=c("antisense","convergent","convergent","divergent","divergent",
  "intergenic","intronic","intronic","sense","sense","sense","sense",
  "sense")

gene_info_ncRNAs$simpl_class=simpl_classes[match(gene_info_ncRNAs$class,
                                                 classes)]

g=ggplot(gene_info_ncRNAs,aes(x=biotype,fill=simpl_class)) +
  geom_bar(position = "fill") +theme_classic() +facet_wrap(~ exonic_type)+
  scale_fill_manual(values = mycolors) +
  theme(text = element_text(size=20))
g
save_tiff_svg(g,
              filename = "Class_distribution.barplot.facet_exonic_type")

g=ggplot(gene_info_ncRNAs,aes(x=exonic_type ,fill=simpl_class)) +
  geom_bar(position = "fill") +theme_classic() +facet_wrap(~ biotype)+
  scale_fill_manual(values = mycolors) +
  theme(text = element_text(size=20))

g
save_tiff_svg(g,
              filename = "Class_distribution.barplot.facet_biotype")

g=ggplot(gene_info_ncRNAs%>%filter(biotype=="lncRNA"),aes(x=exonic_type ,fill=simpl_class)) +
  geom_bar(position = "dodge") +theme_classic() +facet_wrap(~ biotype)+
  scale_fill_manual(values = mycolors) +
  theme(text = element_text(size=20))

g
save_tiff_svg(g,filename = "Barplot_classes_lncRNAsonly")

g=ggplot(gene_info_ncRNAs%>%filter(biotype=="PotNovel"),aes(x=exonic_type ,fill=simpl_class)) +
  geom_bar(position = "dodge") +theme_classic() +facet_wrap(~ biotype)+
  scale_fill_manual(values = mycolors) +
  theme(text = element_text(size=20))

g
save_tiff_svg(g,filename = "Barplot_classes_PotNovelonly")


g=ggplot(gene_info_ncRNAs%>%filter(biotype=="TEC"),aes(x=exonic_type ,fill=simpl_class)) +
  geom_bar(position = "dodge") +theme_classic() +facet_wrap(~ biotype)+
  scale_fill_manual(values = mycolors) +
  theme(text = element_text(size=20))

g

save_tiff_svg(g,filename = "Barplot_classes_TEConly")
g=ggplot(gene_info_ncRNAs%>%filter(biotype=="pseudogene"),aes(x=exonic_type ,fill=simpl_class)) +
  geom_bar(position = "dodge") +theme_classic() +facet_wrap(~ biotype)+
  scale_fill_manual(values = mycolors) +
  theme(text = element_text(size=20))

g
save_tiff_svg(g,filename = "Barplot_classes_pseudogeneonly")

# expression by class ----
g=ggplot(gene_info_ncRNAs,
         aes(x=simpl_class,y=mean_LSK_TPM,fill=exonic_type))+
  geom_boxplot() +theme_classic()+scale_fill_manual(values = mycolors[c(6,5)])+
  stat_compare_means() +scale_y_log10()

g
save_tiff_svg(g,"tmp_plots/","boxplot_expression_by_classandexonic_type_ncRNAs")

g=ggplot(gene_info_ncRNAs,
         aes(col=simpl_class,x=mean_LSK_TPM))+
  stat_ecdf() +theme_classic()+scale_color_manual(values = mycolors)+
  scale_x_log10()


g
save_tiff_svg(g,"tmp_plots/","ecdf_expression_by_classandexonic_type_ncRNAs")

### and exonic type ----

g=ggplot(gene_info_ncRNAs,
       aes(x=exonic_type,y=mean_LSK_TPM,fill=exonic_type))+
  geom_boxplot() +theme_classic() + facet_wrap(~ simpl_class)+
  scale_y_log10() +scale_fill_manual(values = mycolors[c(6,5)])+
  stat_compare_means()
  
g

dir.create("tmp_plots")
save_tiff_svg(g,"tmp_plots/","Boxplot_expression_by_class_ncRNAs")

g=ggplot(gene_info_ncRNAs,
       aes(x=mean_LSK_TPM,col=exonic_type))+
  stat_ecdf() +theme_classic() + facet_wrap(~ simpl_class)+
  scale_x_log10() 
g
save_tiff_svg(g,"tmp_plots/","ecdf_expression_by_class_ncRNAs")


g=ggplot(gene_info_ncRNAs%>%filter(biotype=="PotNovel"),
       aes(x=simpl_class,y=mean_LSK_TPM,fill=simpl_class)) +
  scale_fill_manual(values = mycolors)+
  geom_boxplot() +theme_classic()+facet_wrap(~ exonic_type) +
  scale_y_log10() +ggtitle("Potential novel loci")+
  theme(text=element_text(size=20))
g
save_tiff_svg(g,"tmp_plots/","boxplot_expression_by_class_PotNovel")

g=ggplot(gene_info_ncRNAs%>%filter(biotype=="lncRNA"),
         aes(x=simpl_class,y=mean_LSK_TPM,fill=simpl_class)) +
  scale_fill_manual(values = mycolors)+
  geom_boxplot() +theme_classic()+facet_wrap(~ exonic_type) +
  scale_y_log10() +ggtitle("annotated lncRNAs")+
  theme(text=element_text(size=20))
g
save_tiff_svg(g,"tmp_plots/","boxplot_expression_by_class_PotNovel")

#try vst
# characterize the highly expressed potnovel:
# what marks do they overlap? conservation, etc

## conservation by class ----
g=ggplot(gene_info_ncRNAs,
         aes(x=simpl_class,y=mean_phastCons35way,fill=exonic_type))+
  geom_boxplot() +theme_classic()+scale_fill_manual(values = mycolors[c(6,5)])+
  stat_compare_means()

g
save_tiff_svg(g,"tmp_plots/","boxplot_conservation_by_classandexonic_type_ncRNAs")

g=ggplot(gene_info_ncRNAs,
         aes(x=simpl_class,y=mean_phastCons35way,fill=exonic_type))+
  geom_boxplot() +theme_classic()+scale_fill_manual(values = mycolors[c(6,5)])+
  stat_compare_means()

g
save_tiff_svg(g,"tmp_plots/","boxplot_conservation_by_classandexonic_type_ncRNAs")

### and exonic type ----
g=ggplot(gene_info_ncRNAs,
         aes(x=exonic_type,y=mean_phastCons35way,fill=exonic_type))+
  geom_boxplot() +theme_classic() + facet_wrap(~ simpl_class)+
  scale_y_log10() +scale_fill_manual(values = mycolors[c(6,5)])+
  stat_compare_means()

g

g=ggplot(gene_info_ncRNAs,
         aes(x=mean_phastCons35way,col=exonic_type))+
  stat_ecdf() +theme_classic() + facet_wrap(~ simpl_class)+
  scale_x_log10() 
g
save_tiff_svg(g,"tmp_plots/","ecdf_genebody_phastcons_by_class_ncRNAs")


# overlap with marks ----
# how many lncRNAs overlap CAGE by class and include PCG
# how many lncRNA overlap enhancers 
# how many are likely enhancers?

# overlap with repeats ----
repeat_ol=
