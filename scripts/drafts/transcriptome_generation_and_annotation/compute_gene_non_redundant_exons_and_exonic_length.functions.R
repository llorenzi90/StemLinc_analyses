compute_gene_non_redundant_exonic_length <- function(in_gtf){
  require(tidyverse)
  require(GenomicFeatures)
  in_gtf_exons=in_gtf %>% filter(type=="exon")
  
  txdb <-makeGRangesFromDataFrame(as.data.frame(in_gtf_exons),keep.extra.columns = T)
  
  per_gene_merged_exons <- as.data.frame(reduce(split(txdb,txdb$gene_id)))
  
  gene_exonic_length <- aggregate(per_gene_merged_exons$width,
                                  by=list(gene_id=per_gene_merged_exons$group_name),
                                  function(x)sum(x))
  colnames(gene_exonic_length)[2]="exonic_length"
  
  return(gene_exonic_length)
}

compute_per_transcript_exonic_length <- function(in_gtf){
  require(tidyverse)
  in_gtf_exons=in_gtf %>% filter(type=="exon")
  in_gtf_exons$lens=in_gtf_exons$end - in_gtf_exons$start+1
  tr_len=in_gtf_exons%>%group_by(transcript_id)%>%summarise(len=sum(lens))
  return(tr_len)
}
get_per_gene_merged_exons <- function(in_gtf){
  require(tidyverse)
  require(GenomicFeatures)
  in_gtf_exons=in_gtf %>% filter(type=="exon")
  
  txdb <-makeGRangesFromDataFrame(as.data.frame(in_gtf_exons),keep.extra.columns = T)
  
  per_gene_merged_exons <- as.data.frame(IRanges::reduce(split(txdb,txdb$gene_id)))
  
  return(per_gene_merged_exons)
}


get_per_gene_merged_transcripts <- function(in_gtf,gene_id_col="gene_id"){
  require(tidyverse)
  require(GenomicFeatures)
  in_gtf_exons=in_gtf %>% filter(type=="transcript")
  
  txdb <-makeGRangesFromDataFrame(as.data.frame(in_gtf_exons),keep.extra.columns = T)
  
  per_gene_merged_transcripts <- as.data.frame(IRanges::reduce(split(txdb,txdb$gene_id)))
  
  return(per_gene_merged_transcripts)
}
# function to compute overlap with repeats including total gene fraction
# covered by repeats:

get_overlap_with_repeats_from_GTF=function(gtf){
  require(bedtoolsr)
  require(tidyverse)
  pgme=get_per_gene_merged_exons(gtf)
  pgme=pgme[,c(3:5,2,6,7)]
  pgme=pgme[order(pgme$seqnames,pgme$start),]
  exonic_length=pgme%>%group_by(group_name)%>%
    summarise(exonic_length=sum(width))
  ol=bt.intersect(a = pgme,b=repeats,wo = T)
  ol=ol%>% rowwise()%>% mutate(st=max(c(V2,V8)),en=min(c(V3,V9)))
  olsm=ol%>%group_by(V4)%>%
    summarise(repeats=paste(unique(V10[order(-V13)]),
                            collapse = ","),
              repeats_classes=paste(unique(V11[order(-V13)]),
                                    collapse = ","),
              t.overlap_length=sum(IRanges::width(IRanges::reduce(IRanges(
                start = st,end = en)))))
  olsm=left_join(olsm,exonic_length,by=c(V4="group_name"))
  olsm$repeat.fraction=olsm$t.overlap_length/olsm$exonic_length
  return(olsm)
}


exons_per_transcript=function(gtf, anno_col="biotype"){
  ept <-  gtf%>%filter(type=="exon")%>%group_by(transcript_id)%>% 
    summarise(N_exons=n())
  
  if(anno_col%in%colnames(gtf)){
    gtf=as.data.frame(gtf)
    ept$biotype=gtf[,anno_col][match(ept$transcript_id,
                                     gtf$transcript_id)]
    
    colnames(ept)[ncol(ept)]=anno_col
  }
  return(ept)
  
}


exons_per_gene=function(gtf, gene_col="gene_id", anno_col="biotype"){
  gtf$exon_id=paste0(gtf$seqid,":",gtf$start,"-",gtf$end,":",gtf$strand)
  gtf=as.data.frame(gtf)
  colnames(gtf)[colnames(gtf)==gene_col]="gene_id"
  epg <-  gtf%>%filter(type=="exon")%>%group_by(gene_id)%>% 
    summarise(N_exons=length(unique(exon_id)))
  
  if(anno_col%in%colnames(gtf)){
    gtf=as.data.frame(gtf)
    epg$biotype=gtf[,anno_col][match(epg$gene_id,
                                     gtf[,gene_col])]
    
    colnames(epg)[ncol(epg)]=anno_col
  }
  colnames(epg)[colnames(epg)=="gene_id"]=gene_col
  
  return(epg)
  
}


transcripts_per_gene=function(gtf, gene_col="gene_id", anno_col="biotype"){
  gtf=as.data.frame(gtf)
  colnames(gtf)[colnames(gtf)==gene_col]="gene_id"
  tpg <-  gtf%>%group_by(gene_id)%>% 
    summarise(N_transcripts=length(unique(transcript_id)))
  
  if(anno_col%in%colnames(gtf)){
    gtf=as.data.frame(gtf)
    tpg$biotype=gtf[,anno_col][match(tpg$gene_id,
                                     gtf[,gene_col])]
    
    colnames(tpg)[ncol(tpg)]=anno_col
  }
  colnames(tpg)[colnames(tpg)=="gene_id"]=gene_col
  
  return(tpg)
  
}
