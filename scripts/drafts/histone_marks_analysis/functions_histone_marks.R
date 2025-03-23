## Setup ---------------------------

options(scipen = 999)
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
library(rtracklayer)
library(bedtoolsr)
## Load data---------------------------
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/")
source("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/scripts/StemLinc_code/new_scripts/compute_gene_non_redundant_exons_and_exonic_length.functions.R")
gene_info=read.csv("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/gene_info_PCGs_lncRNAs_pseudo_unfiltered.240511.csv")

primary_chrs=paste0("chr",c(rep(1:19),
                            "X","Y","M"))
# get TSS and promoteres for all genes ----

GTF=readGFF("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/data/various_gtf_files/merged_PCG_lncRNA_pseudo_transcriptome.gtf")


#

#GTF$start=GTF$start -1

# View(gene_info%>%filter(strand=="*"))
# genes_strands=c("XLOC_007631"="-",
#   "XLOC_068623"="+",
#   "XLOC_068641"="+",
#   "XLOC_068650"="+",
#   "XLOC_095109"="+",
#   "XLOC_106083"="+",
#   "XLOC_047358"="+",
#   "XLOC_047370"="-",
#   "XLOC_053951"="+",
#   "XLOC_053960"="+"
#   )
#
# gene_info$strand[match(names(genes_strands),gene_info$gene_name)] <- genes_strands
#

# define TSS and promoter regions ----

#first_exons=GTF%>%filter(exon_number==1) I realised this is not true for
# StringTie transcritps!!

get_first_exons <- function(GTF){
  exons=GTF%>%filter(type=="exon")
  exons$rowID=1:nrow(exons)
  exons_plus=exons%>%filter(strand=="+")  %>%group_by(transcript_id)%>%
    filter(start==min(start))
  exons_minus=exons%>%filter(strand=="-")  %>%group_by(transcript_id)%>%
    filter(end==max(end))
  first_exons=exons%>%filter(rowID%in%c(exons_plus$rowID,exons_minus$rowID))
  return(first_exons)
}

first_exons=get_first_exons(GTF)
first_exons%>%group_by(strand)%>%summarise(n=length(unique(gene_name)))

first_exons$start=first_exons$start - 1
# proximal promoter regions -500bp and +250 bp ----
first_exons_plus=first_exons%>%filter(strand=="+")
first_exons_plus=first_exons_plus%>%group_by(gene_name)%>% filter(!duplicated(start))
TSS_plus=first_exons_plus
TSS_plus$end=TSS_plus$start+1
first_exons_plus$end=first_exons_plus$start + 250
first_exons_plus$start=first_exons_plus$start -500

first_exons_minus=first_exons%>%filter(strand=="-")%>%group_by(gene_name)%>%filter(!duplicated(end))
TSS_minus=first_exons_minus
TSS_minus$start=TSS_minus$end -1

first_exons_minus$start=first_exons_minus$end - 250
first_exons_minus$end=first_exons_minus$end + 500

Promoters=rbind(first_exons_plus,first_exons_minus)

Promoters$type="promoter"
Promoterstw=Promoters
Promoterstw$start=Promoterstw$start+1
export(Promoterstw,"data/various_gtf_files/merged_PCG_lncRNA_pseudo_transcriptome.PROMOTERS.gtf",format = "GTF")

TSSs=rbind(TSS_plus,TSS_minus)

# functions ----

gtf2bed <- function(gtf,gi="gene_name",move_start=1){
  cols=c("seqid","start","end",gi,"score","strand")
  cols=cols[cols%in%colnames(gtf)]
  bed=gtf[,cols]
  bed$start=bed$start - move_start
  return(bed)
}

sort_bed <- function(bed,primary_chrs=paste0("chr",c(rep(1:19),
                                                     "X","Y","M"))){

  bed=as.data.frame(bed)

  bed=bed[order(match(bed[,1],
                      primary_chrs),
                bed[,2]),]
  return(bed)
}

merge_beds <- function(bed,
                       stranded=F,
                       primary_chrs=paste0("chr",c(rep(1:19),
                                                   "X","Y","M"))){

  bed=sort_bed(bed,primary_chrs)

  merged=bt.merge(bed,s = stranded)
  return(merged)
}

intersect_wo <- function(qbed,markbed){

  ol=bt.intersect(a = sort_bed(qbed),
                  b= sort_bed(markbed) , wo = T)
  return(ol)

}

intersect_wao <- function(qbed,markbed){

  ol=bt.intersect(a = sort_bed(qbed),
                  b= sort_bed(markbed) , wao = T)
  return(ol)

}

geneBody_intersect_fromGTF <- function(GTF,
                                       bed2merge,
                                       gene_level=T){
  ## Generate gene body coverage of marks:

  # Make overlap at exon level,
  # then aggregate the total covered length
  # for each transcript and divide by the transcript length

  # calculate :
  # transcript length (sum of it exons lengths)

  require(GenomicFeatures)
  merged_beds=merge_beds(bed2merge)
  mbed_cols=ncol(merged_beds)

  GTF=GTF%>%filter(type=="exon")

  trgtf= gtf2bed(GTF,
                 gi= "transcript_id")

  qbed_cols=ncol(trgtf)

  tr_bd_intersect=intersect_wao(trgtf,
                                merged_beds )

  tr_bd_intersect$exon_len=tr_bd_intersect$V3 - tr_bd_intersect$V2


  colnames(tr_bd_intersect)[qbed_cols+mbed_cols+1]="cov"
  tr_bd_intersect=
    tr_bd_intersect%>%
    group_by(V4)%>%
    summarise(total_cov=sum(cov),
              total_len=sum(exon_len),
              fraction_covered=total_cov/total_len)

  tr_bd_intersect$gene_name=GTF$gene_name[match(tr_bd_intersect$V4,
                                                GTF$transcript_id)]
  if(!gene_level){
    return(tr_bd_intersect)
  }else{
    max_fraction_covered=
      tr_bd_intersect%>%
      group_by(gene_name)%>%
      arrange(-fraction_covered)%>%
      summarise(max_fraction_covered=max(fraction_covered),
                max_tr_covered=V4[1],
                max_tr_len=total_len[1])

    gene_merged_exons=get_per_gene_merged_exons(GTF)
    colnames(gene_merged_exons)[6]="score"
    colnames(gene_merged_exons)[2]="gene_name"
    colnames(gene_merged_exons)[3]="seqid"

    gene_bd_intersect=intersect_wao(gtf2bed(gene_merged_exons),
                                    merged_beds)

    colnames(gene_bd_intersect)[qbed_cols+mbed_cols+1]="cov"
    gene_bd_intersect=
      gene_bd_intersect%>%
      group_by(V4)%>%
      summarise(total_cov=sum(cov),
                total_len=sum(V5),
                total_fraction_covered=total_cov/total_len)

    colnames(gene_bd_intersect)[1]="gene_name"
    gene_bd_intersect=left_join(gene_bd_intersect,
                                max_fraction_covered)

    return(gene_bd_intersect)
  }

  # returns
  # 1. max_fraction_covered:
  # for the transcript that is most highly covered
  # 2. total_gene_fraction_covered: covered length of merged
  # exonic sequence over total length
  # 3. total length covered
}


get_near_peaks_unstranded <- function(bed,mark,w=1000,cols="2,3,8,9"){
  win=bt.window(a = bed,mark,w = w)
  win.ov=bt.overlap(win,cols=cols)
  return(win.ov)
}

get_near_peaks_stranded <- function(bed,mark,cols="2,3,8,9",lwin=1000,rwin=0,sm=T){
  win=bt.window(a = bed,mark,l = lwin,r = rwin,sw = T,sm=sm)
  win.ov=bt.overlap(win,cols=cols)
  return(win.ov)
}
