## Head -------------------------------------
##
##
## Purpose of script: make venn diagram of overlap of genes with reference
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-05-13
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
library(rtracklayer)
## Load data---------------------------
source("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/nejm_palette.R")

setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/")
gtf_gencode_lncRNAs=readGFF("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/mouse/GENCODE_Release_M31_GRCm39_genomeAND_annot_for_marie_curie/gencode.vM31.long_noncoding_RNAs.gtf")
gene_info=read.csv("gene_info_PCGs_lncRNAs_pseudo_unfiltered.240511.csv")
# gencode overlap with StemLinc ----
mycolors=nejm_pal[c(1,3,4,5,2,6,7,8)]

ol_gencode <- read.table('/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/overlap_with_GENCODE_lncRNAs/gencodelncRNAs_overlap_StemLinc.gene_level.txt',header = T)
table(ol_gencode$N_exact_matchs==3)

colnames(ol_gencode)[4:8] <- c("NStemLinc_exact_match",
                               "NStemLinc_other_match",
                               "NStemLinc_any_match",
                               "NStemLinc_antisense",
                               "max_Nexons_StemLinc")

table(ol_gencode$gene_id%in%gene_info$gene_id)

# gene_info=left_join(gene_info,ol_gencode %>% 
#                       dplyr::select(gene_id,
#                                     NStemLinc_exact_match,
#                                     NStemLinc_other_match,
#                                     NStemLinc_any_match,
#                                     NStemLinc_antisense,
#                                     max_Nexons_StemLinc))




ol_RefSeq <- read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/APRIL/RefSeq_ol_StemLinc/RefSeq2add_overlap_StemLinc.gene_level.txt",header = T)
colnames(ol_RefSeq)[4:8] <- c("NStemLinc_exact_match",
                               "NStemLinc_other_match",
                               "NStemLinc_any_match",
                               "NStemLinc_antisense",
                               "max_Nexons_StemLinc")
table(ol_RefSeq$gene_name%in%gene_info$gene_name)
refseqgtf=readGFF("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/data/various_gtf_files/RefSeq_finaltranscripts2add.gtf")
gi=ol_RefSeq$gene_name[!ol_RefSeq$gene_name%in%gene_info$gene_name]
unique(refseqgtf$seqid[refseqgtf$gene_name%in%gi])
View(refseqgtf%>%filter(gene_name%in%gi))
ol_RefSeq$gene_id=gene_info$gene_id[match(ol_RefSeq$gene_name,
                                          gene_info$gene_name)]
ol_RefSeq$gene_id[is.na(ol_RefSeq$gene_id)] <- ol_RefSeq$gene_name[is.na(ol_RefSeq$gene_id)]

ol_Refs=rbind(ol_gencode%>%dplyr::select(gene_id,
                                         NStemLinc_exact_match,
                                         NStemLinc_other_match,
                                         NStemLinc_any_match,
                                         NStemLinc_antisense,
                                         max_Nexons_StemLinc),
              ol_RefSeq%>%dplyr::select(gene_id,
                                        NStemLinc_exact_match,
                                        NStemLinc_other_match,
                                        NStemLinc_any_match,
                                        NStemLinc_antisense,
                                        max_Nexons_StemLinc))
  
gene_info=left_join(gene_info,ol_Refs)

table(gene_info$NStemLinc_any_match==3,gene_info$passfilter)
gene_info_april=read.delim("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/APRIL/gene_level_info_APRIL.txt")
table(gene_info_april$gene_name%in%gene_info$gene_name)
gene_info=left_join(gene_info, gene_info_april%>%
                      dplyr::select(gene_name,
                                    Delas_overlap,
                                    Luo_gene))

gene_info_filt=gene_info%>%filter(passfilter)

table(gene_info_filt$NStemLinc_any_match,
      !is.na(gene_info_filt$Delas_overlap),
      !is.na(gene_info_filt$Luo_gene))

gene_info_filt$InRef=!is.na(gene_info_filt$NStemLinc_any_match)&gene_info_filt$NStemLinc_any_match==3
gene_info_filt_lncRNA_annotated=gene_info_filt%>%filter(biotype%in%c("lncRNA"))

table(gene_info_filt_lncRNA_annotated$NStemLinc_any_match==3)

gene_info_lncRNAs=gene_info%>%filter(biotype%in%c("lncRNA","PotNovel","TEC"))

table(gene_info_lncRNAs$NStemLinc_any_match==3)
gene_info_lncRNAs$NStemLinc_any_match[gene_info_lncRNAs$biotype=="PotNovel"]=3
table(gene_info_lncRNAs$NStemLinc_any_match==3,gene_info_lncRNAs$passfilter)



ol_RefSeq$biotype=gene_info$biotype[match(ol_RefSeq$gene_id,gene_info$gene_id)]
table(ol_RefSeq$biotype)


gene_info_lncRNAs_inSL=gene_info_lncRNAs[gene_info_lncRNAs$NStemLinc_any_match==3&
                                           !is.na(gene_info_lncRNAs$NStemLinc_any_match),]
table(gene_info_lncRNAs_inSL$biotype,gene_info_lncRNAs_inSL$passfilter)
gene_info_lncRNAs_inSL

gcRefseq_len=length(unique(gtf_gencode_lncRNAs$gene_id))+sum(ol_RefSeq$biotype=="lncRNA",na.rm = T)
Delas_Luo_length=3369+159

StemLinc_genes=gene_info_lncRNAs_inSL$gene_id

GENCODE_RefSeq_genes=gene_info_lncRNAs_inSL$gene_id[!gene_info_lncRNAs_inSL$biotype=="PotNovel"]
genes2add=gcRefseq_len - length(GENCODE_RefSeq_genes)
GENCODE_RefSeq_genes=c(GENCODE_RefSeq_genes,paste0("GR_",seq(1:genes2add)))


Delas_Luo_genes=gene_info_lncRNAs_inSL$gene_id[!is.na(gene_info_lncRNAs_inSL$Delas_overlap)|
                                                 !is.na(gene_info_lncRNAs_inSL$Luo_gene)]
genes2add=Delas_Luo_length - length(Delas_Luo_genes)
Delas_Luo_genes=c(Delas_Luo_genes,paste0("DL_",seq(1:genes2add)))

# library(VennDiagram)
# 
# # Generate 3 sets of 200 words
# set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
# set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
# set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
# 
# # Chart
# venn.diagram(
#   x = list(StemLinc_genes,GENCODE_RefSeq_genes,Delas_Luo_genes),
#   category.names = c("StemLinc assembly" , "GENCODE/Refseq lncRNAs" , "Delas/Luo et al. lncRNAs"),
#   filename = '#14_venn_diagramm.png',
#   output=TRUE
# )
# 
# 
# # Helper function to display Venn diagram
# display_venn <- function(x, ...){
#   library(VennDiagram)
#   grid.newpage()
#   venn_object <- venn.diagram(x, filename = NULL, ...)
#   grid.draw(venn_object)
# }
# 
# display_venn(list(StemLinc_genes,GENCODE_RefSeq_genes,Delas_Luo_genes))
# 
# library(ggvenn)
# ggvenn(
#   data = list(StemLinc_genes,GENCODE_RefSeq_genes,Delas_Luo_genes), 
#   fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
#   stroke_size = 0.5, set_name_size = 4
# )

library(gplots)

v.table <- venn(list(StemLinc_genes,GENCODE_RefSeq_genes,Delas_Luo_genes))
# the venn is not the best idea

# pie
prev_sel=read.csv("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/data/lncRNA_candidate_lists/final_summary_candidates_02May24.csv")
gene_info_prev_sel=gene_info%>%filter(gene_name%in%prev_sel$gene_name)
table(prev_sel$gene_name%in%gene_info$gene_name)
table(prev_sel$gene_name%in%gene_info_lncRNAs_inSL$gene_name)
table(prev_sel$gene_name%in%gene_info_lncRNAs_inSL$gene_name[gene_info_lncRNAs_inSL$passfilter])
View(prev_sel[prev_sel$gene_name%in%gene_info_lncRNAs_inSL$gene_name[gene_info_lncRNAs_inSL$passfilter],])

gene_info_lncRNAs_inSL$inDelas_Luo=ifelse(!is.na(gene_info_lncRNAs_inSL$Delas_overlap),
                                          T,ifelse(!is.na(gene_info_lncRNAs_inSL$Luo_gene),T,F))
table(gene_info_lncRNAs_inSL$inDelas_Luo)
table(gene_info_lncRNAs_inSL$inDelas_Luo,gene_info_lncRNAs_inSL$biotype)

gene_info_lncRNAs_inSL$annot=gene_info_lncRNAs_inSL$biotype
gene_info_lncRNAs_inSL$annot[gene_info_lncRNAs_inSL$biotype=="PotNovel"&
                               gene_info_lncRNAs_inSL$inDelas_Luo]="in_Delas/Luo"

tb=table(gene_info_lncRNAs_inSL$annot)

names(tb)=c("Delas/Luo et al.","GENCODE/RefSeq lncRNAs","StemLinc PotNovel","GENCODE TEC")
tb=tb[c(1,3,4,2)]
pie(tb,col=mycolors[c(5,2,6,1)],labels = paste0(names(tb),"\n(",
                                               tb,")"))
p=pie(tb,col=mycolors[c(5,2,6,1)],labels = paste0(names(tb),"\n(",
                                                  tb,")"))



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

save_tiff_svg(plo = pie(tb,col=mycolors[c(5,2,6,1)],labels = paste0(names(tb),"\n(",
                                                                    tb,")")),
              filename = "test.pie",h=5,w=8)

filename="test.pie"
svg( height = 5 ,width = 5,
     paste0("/",filename,".svg"))
print(plo)


dev.off()



#### expression filter ----

tb=table(gene_info_lncRNAs_inSL$annot[gene_info_lncRNAs_inSL$passfilter])

names(tb)=c("Delas/Luo et al.","GENCODE/RefSeq lncRNAs","StemLinc PotNovel","GENCODE TEC")
tb=tb[c(1,3,4,2)]
pie(tb,col=mycolors[c(5,2,6,1)],labels = paste0(names(tb),"\n(",
                                                tb,")"))
p=pie(tb,col=mycolors[c(5,2,6,1)],labels = paste0(names(tb),"\n(",
                                                  tb,")"))


save_tiff_svg(plo = pie(tb,col=mycolors[c(5,2,6,1)],labels = paste0(names(tb),"\n(",
                                                                    tb,")")),
              filename = "test.pie.filtered",h=5,w=8)


table(gene_info$biotype)
PCGs_passfilter=gene_info%>%filter(biotype=="PCG",passfilter)

PCGs_passfilter$inDelas_Luo=F
PCGs_passfilter$annot="PCG"
gene_info_filt=rbind(gene_info_lncRNAs_inSL%>%filter(passfilter),
                     PCGs_passfilter)
gene_info_filt=gene_info_filt[,-1]
write.csv(gene_info_filt,"data/gene_info_PCGs_lncRNAs_inStemLinc.expression_filtered.240513.csv",row.names = F)

g=ggplot(gene_info_filt,
         aes(x=mean_LSK_TPM,col=annot,linetype=exonic_type))+ 
  stat_ecdf(size=2) +
  theme_classic() + 
  scale_color_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
  ylab("Cumulative fraction")+ scale_x_log10()+
  ggtitle("") + xlab("mean TPM")

g

g=ggplot(gene_info_filt,
         aes(y=mean_LSK_TPM,x=annot,fill=annot))+ 
  geom_violin() +
  theme_classic() + 
  scale_fill_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
  ylab("Cumulative fraction")+ scale_y_log10()+
  ggtitle("") + xlab("mean TPM")

g

g=ggplot(gene_info_filt,
         aes(x=mean_LSK_TPM,col=annot))+ 
  stat_ecdf(size=2) +
  theme_classic() + 
  scale_color_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
  ylab("Cumulative fraction")+ scale_x_log10()+
  ggtitle("") + xlab("mean TPM")

g

g=ggplot(gene_info_filt,
         aes(y=mean_LSK_TPM,x=annot,fill=annot))+ 
  geom_boxplot() +
  theme_classic() + 
  scale_fill_manual(values = mycolors)+
  theme(text = element_text(size = 20))+
  ylab("Cumulative fraction")+ scale_y_log10()+
  ggtitle("") + xlab("mean TPM")

g

# conservation
ggplot(gene_info_filt,
       aes(x=mean_phastCons35way,col=annot))+ 
         stat_ecdf(size=2) +
         theme_classic() + 
         scale_color_manual(values = mycolors)+
         theme(text = element_text(size = 20))+
         ylab("Cumulative fraction")+ scale_x_log10()+
         ggtitle("") + xlab("mean phastCons35way")
       
g

prev_sel_remain=prev_sel%>%filter(gene_name%in%gene_info_filt$gene_name)
prev_sel_remain$annot=gene_info_filt$annot[match(prev_sel_remain$gene_name,
                                                 gene_info_filt$gene_name)]
table(prev_sel_remain$annot)
table(prev_sel_remain$annot,prev_sel_remain$exonic_type)
