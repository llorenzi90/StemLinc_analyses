## Head -------------------------------------
##
##
## Purpose of script: Calculation of average conservation scores for known and novel genes
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-04-11
##
## Email: lucialorenzi90@gmail.com
##
## Notes ---------------------------
## Which measure to use?
## The PhastCons score is a probability that each nucleotide 
## belongs to a conserved element, whereas abs(phyloP) is the -log(p-value) 
## under a null hypothesis of neutral evolution, and a negative sign indicates 
## faster-than expected evolution, while positive values imply conservation.May 27, 2010
## I can try and compare both   
##
## Setup ---------------------------

options(scipen = 999) 
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
library(GenomicScores)
library(rtracklayer)
## Load data---------------------------
setwd("~/work_local/")
avgs=availableGScores()
library(AnnotationHub)
rownames(avgs)[avgs$AnnotationHub]
phastcons <- getGScores("phastCons35way.UCSC.mm39")
source("/home/llorenzi/work_local/scripts/compute_gene_non_redundant_exonic_length.function.R")

# load my gene coords ----
#gtf=readGFF("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/analyses/gffcompare_all/LT-HSC_LSK_MPPs/gtfs/LT-HSC_LSK_MPPs.potentialNovel.NoCodingProb.gtf")
RefSeq_finalgenes2addGTF=readGFF("APRIL/RefSeq_finaltranscripts2add.gtf")

gtf_exons=RefSeq_finalgenes2addGTF %>% filter(type=="exon")

txdb <-makeGRangesFromDataFrame(as.data.frame(gtf_exons),keep.extra.columns = T)
per_gene_merged_exons <- reduce(split(txdb,txdb$gene_id))


# get phastCons scores ----
pcsco <- score(phastcons, unlist(per_gene_merged_exons))
summary(pcsco) 

per_gene_merged_exons=as.data.frame(per_gene_merged_exons)
per_gene_merged_exons$phastCons35way <- pcsco

per_gene_score <- per_gene_merged_exons %>% group_by(group_name) %>% 
  summarise(mean_phastCons35way=mean(phastCons35way,na.rm=T))



# A lot of lncRNAs have NA values
# rtracklayer approach ---- 
per_gene_merged_exons <- reduce(split(txdb,txdb$gene_id))
per_gene_merged_exons <- unlist(per_gene_merged_exons)

# import phastcons bigwig using the lncRNA exonic regions

phast_scores <- import.bw(con = "~/Downloads/mm39.phastCons35way.bw",
                          which=per_gene_merged_exons)
phast_scores$score

# Note that with this approach I get a value per position (bigWig) or per
# run with consecutive equal values (bedGraph), this gives the possibility
# of summarizing the scores in a different way if desired

# for now I just want to perform the average (removing NAs)
# and compare the ones that I obtain through this method
# with the ones that are not NA with the GenomicScores package

anyNA(phast_scores$score)
# there are not NA values, 
# but this is because import.bw() just don't retrieve
# anything if there is no value there
# I have to match the positions with the genes


# get the individual positions id for all nt covered by exons of interest
per_gene_merged_exons <- as.data.frame(reduce(split(txdb,txdb$gene_id)))

gene_posid=apply(per_gene_merged_exons[,c("seqnames","start","end")],1,function(x){
  return(paste(x[1],x[2]:x[3],sep = "_"))
})

gene_posid=unlist(gene_posid)

# get the same info, i.e individual position ids from phast_scores
phast_posid=paste(seqnames(phast_scores),start(phast_scores),sep = "_")

# Use these ids to match scores with original positions
table(gene_posid%in%phast_posid)

length(phast_scores$score)

matched_scores <- phast_scores$score[match(gene_posid,phast_posid,nomatch = NA)]

gene_ids=rep(per_gene_merged_exons$group_name,per_gene_merged_exons$width)
length(gene_ids)

scores_per_gene=data.frame(gene_id=gene_ids,
                           score=matched_scores) %>% group_by(gene_id) %>%
  summarise(mean_score=mean(score,na.rm=T))

table(is.na(scores_per_gene$mean_score))
# dir.create("public_datasets/analyses/gffcompare_all/LT-HSC_LSK_MPPs/conservation")
# write.table(scores_per_gene,"public_datasets/analyses/gffcompare_all/LT-HSC_LSK_MPPs/conservation/LT-HSC_LSK_MPPs.potentialNovel.NoCodingProb.phastCons35way.txt",
#             quote = F,row.names = F,sep = "\t")

dir.create("APRIL/phastCons35way")
scores_per_gene=left_join(scores_per_gene,RefSeq_finalgenes2addGTF%>%select(gene_id,gene_name,biotype))
scores_per_gene=scores_per_gene[!is.na(scores_per_gene$gene_name),]
write.table(scores_per_gene,"APRIL/phastCons35way/RefSeq_finalgenes2add.phastCons35way.txt",
            quote = F,row.names = F,sep = "\t")
