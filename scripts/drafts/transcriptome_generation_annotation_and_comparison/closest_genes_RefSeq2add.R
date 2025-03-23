## ---------------------------
##
##
## Purpose of script: find closest genes including RefSeq genes to add
## 
## Author: Lucia Lorenzi
##
## Date Created: 2024-04-14 
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes: I want to find close genes to lncRNAs
## PCGs, other lncRNAs (known and novel)
##  create a bed file for each class of gene 
## for each lncRNA
## 1- find the 4 closest PCGs, known and novel lncRNAs
## 2- retrieve distance for each of them
## 3- find the closest gene overall and its distance
##
## I'll use bedtools closest with parameters:
##  . -N (Require that the query and the closest hit have different names. For BED, the 4th column is compared.)
##  . -k (Report the k closest hits. Default is 1. If tieMode = “all”, all ties will still be reported.)
##  . 
## ---------------------------

options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
require(rtracklayer)
require(bedtoolsr)
## ---------------------------

# setwd and load input data -----------------------------------------------


setwd("~/work_local/APRIL/")

#setwd("public_datasets/analyses/gffcompare_all/LT-HSC_LSK_MPPs/")
dir.create("closest_genes")

# define primary chromosomes ----
# (used to filter and sort)
primary_chrs=paste0("chr",
       c(1:19,"X","Y","M"))
#read GTF of potential novel genes
#potnovel_gtf=readGFF("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/analyses/gffcompare_all/LT-HSC_LSK_MPPs/LSK_genes/combined_gtfs/LT-HSC_LSK_MPPs.combined.non_ref_GENES_LSK.gtf")
potnovel_gtf=readGFF("/home/llorenzi/work_local/LSK_StemLinc.combined.potNovel_3reps.gtf")

#read GENCODE GTF for known lncRNAs
source("~/work_local/scripts/get_genomic_regions_function.R")
get_genomic_regions
gr_potnovel=get_genomic_regions(potnovel_gtf)

table(gr_potnovel$seqnames)

#GENCODE gtf
#Ref_gtf <- readGFF("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/mouse/GENCODE_Release_M31_GRCm39_genomeAND_annot_for_marie_curie/gencode.vM31.primary_assembly.annotation.gtf")
Ref_bed <- read.table('/home/llorenzi/references/gencode.vM31.primary_assembly.PCG_lncRNA.simplified.tsv',header = T)

# generate bed tables diff biotypes ---------------------------------------

#generate bed format tables for PCGs and lncRNAs from gencode
Ref_bed <- Ref_bed %>% filter(seqnames %in% primary_chrs)


Ref_bed_PCG=Ref_bed %>% filter(!biotype%in%c("TEC","lncRNA"))
Ref_bed_lncRNA <- Ref_bed %>% filter(biotype%in%c("TEC","lncRNA"))

#Add RefSeq genes to Refbed
RefSeqGTF=readGFF("/home/llorenzi/work_local/APRIL/RefSeq_finaltranscripts2add.gtf")

RefSeqbed=get_genomic_regions(RefSeqGTF,by = "gene_name")
RefSeqbed$biotype=RefSeqGTF$biotype[match(RefSeqbed$gene_name,
                                          RefSeqGTF$gene_name)]

#check if some genes are in Ref_bed
table(RefSeqbed$gene_name%in%Ref_bed$gene_name)
dupgenes=RefSeqbed$gene_name[RefSeqbed$gene_name%in%Ref_bed$gene_name]
# these are genes with biotypes in conflict!

# I'll check their expression across samples 
# setwd("/home/llorenzi/work_local/")
# 
# vst=read.table("expression_files/featureCounts/PCGslncRNAsPseudogenes.110424.vst.txt")
# vst$gene_name=rownames(vst)
# 
# plot_gene_Exp=function(ge){
#   vst_gene=vst[rownames(vst)==ge,]
#   vst_gene=pivot_longer(vst_gene,cols = 1:(ncol(vst_gene)-1),names_to = "samples",values_to = "vst")
#   vst_gene$samples=factor(vst_gene$samples,levels = unique(vst_gene$samples))
#   
#   vst_gene$sample_type=sapply(strsplit(as.character(vst_gene$samples),split = "_"),function(x)x[[1]])
#   vst_gene$sample_type=factor(vst_gene$sample_type,levels = unique(vst_gene$sample_type))
#   
#   g=ggplot(vst_gene,aes(x=samples,y=vst,fill=sample_type))+
#     geom_point()+theme_classic() +
#     theme(axis.text.x = element_text(angle = 90))+ggtitle(ge)
#   print(g)
#   
# }
# 
# sapply(dupgenes,plot_gene_Exp)

# "Gm3893" which is now an unprocessed_pseudogene in updated gencode(protein_coding in the one I use)
# and a pseudogene in RefSeq 
# has the most interesting expression pattern
# however when I check its intramodular connectivity it does not look interesting
# intramodular_connectivity <- read.delim("~/work_local/APRIL/WGCNA/intramodular_connectivity.txt")
# Therefore, not any of these genes are really interesting, I can just discard them

RefSeqbed=RefSeqbed%>%filter(!gene_name%in%dupgenes)

#sort bed tables
gr_potnovel <- gr_potnovel[order(gr_potnovel$seqnames,gr_potnovel$start),]
Ref_bed_PCG <- Ref_bed_PCG[order(match(Ref_bed_PCG$seqnames,primary_chrs),Ref_bed_PCG$start),]

Ref_bed_pseudogene=RefSeqbed%>%filter(biotype!="lncRNA")
Ref_bed_pseudogene=Ref_bed_pseudogene[order(match(Ref_bed_pseudogene$seqnames,
                                                  primary_chrs),
                                            Ref_bed_pseudogene$start),]

RefSeqbed_lncRNA=RefSeqbed%>%filter(biotype=="lncRNA")
RefSeqbed_lncRNA$gene_id=RefSeqbed_lncRNA$gene_name
RefSeqbed_lncRNA=RefSeqbed_lncRNA[,colnames(Ref_bed_lncRNA)]
Ref_bed_lncRNA=rbind(Ref_bed_lncRNA,RefSeqbed_lncRNA)
Ref_bed_lncRNA <- Ref_bed_lncRNA[order(match(Ref_bed_lncRNA$seqnames,primary_chrs),Ref_bed_lncRNA$start),]

## add biotype and gene_name to potnovel
gr_potnovel$gene_name=gr_potnovel$gene_id
gr_potnovel$biotype="potnew_lncRNA"

## add missing columns to pseudogenes
Ref_bed_pseudogene$gene_id=Ref_bed_pseudogene$gene_name
Ref_bed_pseudogene=Ref_bed_pseudogene[,colnames(Ref_bed_lncRNA)]

# potNovel bedtools closest genes -----------------------------------------

#compute closest gene in each biotype to each potential novel lncRNAs 
potNovel_closest_PCG=bt.closest(a=gr_potnovel,b=Ref_bed_PCG,d = T,k=4)
potNovel_closest_lncRNA=bt.closest(a=gr_potnovel,b=Ref_bed_lncRNA,d = T,k=4)
potNovel_closest_potNovel=bt.closest(a=gr_potnovel,b=gr_potnovel,d = T,N=T,k=4)
potNovel_closest_pseudo=bt.closest(a=gr_potnovel,b=Ref_bed_pseudogene,d = T,k=4)


#merge target biotypes into a single table
potNovel_closest_genes=rbind(potNovel_closest_PCG,
                             potNovel_closest_lncRNA,
                             potNovel_closest_potNovel,
                             potNovel_closest_pseudo)


#add/change column names 
cns=c("chrom","query_start","qery_end", 
      "q_gene_id","q_width","q_strand","q_gene_name","q_biotype","t_chrom",
      "t_start","t_end","t_gene_id","t_width","t_strand",
      "t_gene_name","t_biotype","distance")


colnames(potNovel_closest_genes)=cns


# lncRNAs bedtools closest ------------------------------------------------

lncRNA_closest_PCG=bt.closest(a=Ref_bed_lncRNA,b=Ref_bed_PCG,d = T,N=T,k=4)
lncRNA_closest_lncRNA=bt.closest(a=Ref_bed_lncRNA,b=Ref_bed_lncRNA,d = T,N=T,k=4)
lncRNA_closest_potNovel=bt.closest(a=Ref_bed_lncRNA,b=gr_potnovel,d = T,N=T,k=4)
lncRNA_closest_pseudo=bt.closest(a=Ref_bed_lncRNA,b=Ref_bed_pseudogene,d = T,N=T,k=4)

#merge target biotypes into a single table
lncRNA_closest_genes=rbind(lncRNA_closest_PCG,
                           lncRNA_closest_lncRNA,
                           lncRNA_closest_potNovel,
                           lncRNA_closest_pseudo)


#add/change column names 
colnames(lncRNA_closest_genes)=cns


# pseudogenes bedtools closest ------------------------------------------------

pseudo_closest_PCG=bt.closest(a=Ref_bed_pseudogene,b=Ref_bed_PCG,d = T,N=T,k=4)
pseudo_closest_lncRNA=bt.closest(a=Ref_bed_pseudogene,b=Ref_bed_lncRNA,d = T,N=T,k=4)
pseudo_closest_potNovel=bt.closest(a=Ref_bed_pseudogene,b=gr_potnovel,d = T,N=T,k=4)
pseudo_closest_pseudo=bt.closest(a=Ref_bed_pseudogene,b=Ref_bed_pseudogene,d = T,N=T,k=4)

#merge target biotypes into a single table
pseudo_closest_genes=rbind(pseudo_closest_PCG,
                           pseudo_closest_lncRNA,
                           pseudo_closest_potNovel,
                           pseudo_closest_pseudo)


#add/change column names 
colnames(pseudo_closest_genes)=cns


# merge bedtool results ---------------------------------------------------

knownANDpotNovel_closest_genes <- rbind(potNovel_closest_genes,
                                        lncRNA_closest_genes,
                                        pseudo_closest_genes)


# write table -------------------------------------------------------------
write.table(knownANDpotNovel_closest_genes,
            "~/work_local/APRIL/closest_genes/knownANDpotNovel_lncRNAs_pseudo_closest.txt",
            sep = "\t", quote = F,row.names = F)



# # genes within 10kb potNovel ------------------------------------------------------
# #get all the genes within 10 kb of each lncRNA (if any)
# 
# # for this, first create a bed file with flanking regions for each query gene 
# 
# #mm39_genome=read.table('/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/mouse/GRCm39_primary_chrs.bed')
# potnovel_10k_flank=gr_potnovel
# potnovel_10k_flank$start=potnovel_10k_flank$start - 10000
# potnovel_10k_flank$end=potnovel_10k_flank$end + 10000
# 
# #then run bedtools intersect with the extended regions against each 
# # biotype's genes
# PCG_within_10k_potNovel=bt.intersect(a = potnovel_10k_flank,b = Ref_bed_PCG,wo = T)
# lncRNA_within_10k_potNovel=bt.intersect(a = potnovel_10k_flank,b = Ref_bed_lncRNA,wo = T)
# potNovel_within_10k_potNovel=bt.intersect(a = potnovel_10k_flank,b = gr_potnovel,wo = T)
# 
# potNovel_within_10k_potNovel <- potNovel_within_10k_potNovel %>%filter(V4!=V12)
# 
# 
# genes_within_10k_potNovel <- rbind(PCG_within_10k_potNovel,
#                                    lncRNA_within_10k_potNovel,
#                                    potNovel_within_10k_potNovel)
# 
# 
# # genes within 10kb lncRNA ------------------------------------------------------
# #get all the genes within 10 kb of each lncRNA (if any)
# 
# # for this, first create a bed file with flanking regions for each query gene 
# 
# lncRNA_10k_flank=Ref_bed_lncRNA
# lncRNA_10k_flank$start=lncRNA_10k_flank$start - 10000
# lncRNA_10k_flank$end=lncRNA_10k_flank$end + 10000
# 
# #then run bedtools intersect with the extended regions against each 
# # biotype's genes
# PCG_within_10k_lncRNA=bt.intersect(a = lncRNA_10k_flank,b = Ref_bed_PCG,wo = T)
# lncRNA_within_10k_lncRNA=bt.intersect(a = lncRNA_10k_flank,b = Ref_bed_lncRNA,wo = T)
# potNovel_within_10k_lncRNA=bt.intersect(a = lncRNA_10k_flank,b = gr_potnovel,wo = T)
# 
# lncRNA_within_10k_lncRNA <- lncRNA_within_10k_lncRNA %>% filter(V4!=V12)
# 
# 
# genes_within_10k_lncRNA <- rbind(PCG_within_10k_lncRNA,
#                                    lncRNA_within_10k_lncRNA,
#                                    potNovel_within_10k_lncRNA)
# 
# 
# 
# # merge genes within 10k tables -------------------------------------------
# 
# genes_within_10k_knownANDpotNovel_lncRNAs <- rbind(genes_within_10k_potNovel,
#                                                    genes_within_10k_lncRNA)
# 
# 
# 
# # write table genes 10k ---------------------------------------------------
# 
# write.table(genes_within_10k_knownANDpotNovel_lncRNAs,
#             "closest_genes/genes_within_10k_knownANDpotNovel_lncRNAs.txt",
#             sep = "\t", quote = F,row.names = F)
# 
# 
# 
