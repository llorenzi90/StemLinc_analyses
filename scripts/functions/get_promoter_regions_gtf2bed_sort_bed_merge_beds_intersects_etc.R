# functions extracted from ol_multiple_public_histone_peaks.R
#source('/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/scripts/StemLinc_code/compute_gene_non_redundant_exons_and_exonic_length.functions.R')
library(trastools)
source("scripts/functions/functions.R")
primary_chrs=paste0("chr",c(rep(1:22),
                            "X","Y","M"))

get_first_exons <- function(GTF){
  exons=GTF%>%filter(type=="exon")
  exons$rowID=1:nrow(exons)
  exons_plus=exons%>%filter(strand=="+")  %>%group_by(transcript_id)%>%
    filter(start==min(start))
  exons_minus=exons%>%filter(strand=="-")  %>%group_by(transcript_id)%>%
    filter(end==max(end))
  first_exons=exons%>%filter(rowID%in%c(exons_plus$rowID,exons_minus$rowID))
  first_exons=first_exons[,-ncol(first_exons)]
  return(first_exons)
}

get_promoter_regions <- function(GTF,preTSS=500,postTSS=250){
  first_exons=get_first_exons(GTF)

  # promoter regions, default -500bp and +250 bp
  first_exons_plus=first_exons%>%filter(strand=="+")

  first_exons_plus$end=first_exons_plus$start + postTSS
  first_exons_plus$start=first_exons_plus$start - preTSS

  first_exons_minus=first_exons%>%filter(strand=="-")

  first_exons_minus$start=first_exons_minus$end - postTSS
  first_exons_minus$end=first_exons_minus$end + preTSS

  Promoters=rbind(first_exons_plus,first_exons_minus)

  Promoters$type="promoter"
  Promoters=Promoters[match(GTF$transcript_id,Promoters$transcript_id,nomatch = 0),]
  Promoters=Promoters[!duplicated(Promoters$transcript_id),]
  return(Promoters)
}

gtf2bed <- function(gtf,gi="gene_name",move_start=1){
  cols=c("seqid","start","end",gi,"score","strand")
  cols=cols[cols%in%colnames(gtf)]
  bed=gtf[,cols]
  bed$start=bed$start - move_start
  return(bed)
}

sort_bed <- function(bed,primary_chrs=paste0("chr",c(rep(1:22),
                                                     "X","Y","M"))){

  bed=as.data.frame(bed)

  bed=bed[order(match(bed[,1],
                      primary_chrs),
                bed[,2]),]
  return(bed)
}

merge_beds <- function(bed,
                       stranded=F,
                       primary_chrs=paste0("chr",c(rep(1:22),
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
                                       gene_level=T,
                                       gene_col="gene_id"){
  ## Generate gene body coverage of marks:

  # Make overlap at exon level for each transcript
  # then aggregate the total covered length
  # for each transcript and divide by the transcript length

  # calculate :
  # transcript length (sum of it exons lengths)

  require(GenomicFeatures)
  merged_beds=merge_beds(bed2merge)
  mbed_cols=ncol(merged_beds)

  GTF=GTF%>%filter(type=="exon")

  trgtf= extract_bed(GTF,
                   name= "transcript_id")

  qbed_cols=ncol(trgtf)

  tr_bd_intersect=intersect_wao(trgtf,
                                merged_beds )

  tr_bd_intersect$exon_len=tr_bd_intersect$V3 - tr_bd_intersect$V2


  colnames(tr_bd_intersect)[qbed_cols+mbed_cols+1]="cov"
  tr_bd_intersect=
    tr_bd_intersect%>%
    dplyr::group_by(V4)%>%
    dplyr::summarise(total_cov=sum(cov),
              total_len=sum(exon_len),
              fraction_covered=total_cov/total_len)

  GTF$gene_name=GTF[,colnames(GTF)==gene_col]
  tr_bd_intersect$gene_name=GTF$gene_name[match(tr_bd_intersect$V4,
                                                GTF$transcript_id)]
  GTF$gene_id=GTF$gene_name
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

    gene_merged_exons=trastools::get_per_gene_merged_exons(GTF)
    colnames(gene_merged_exons)[6]="score"
    colnames(gene_merged_exons)[2]="gene_name"
    colnames(gene_merged_exons)[3]="seqid"

    gene_bd_intersect=intersect_wao(extract_bed(gene_merged_exons,
                                                name = "gene_name"),
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
