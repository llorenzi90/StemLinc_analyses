require(data.table)
require(GenomicRanges)
get_genomic_regions=function(gffdf,by="gene_id"){
  df=as.data.frame(gffdf)
  df$gene_id=df[,by]
  df=df[!is.na(df$gene_id),]
  df=as.data.table(df)
  df[,':='(start=min(start),end=max(end)),by=gene_id]
  df=df[!duplicated(df$gene_id),]
  genecoords=makeGRangesFromDataFrame(df[,c("seqid" , "start","end","strand","gene_id")],keep.extra.columns = T)
  genecoords=as.data.frame(genecoords)
  colnames(genecoords)[6]=by
  #change format to bed file
  genecoords=genecoords[,c(1:3,6,4,5)]
  return(genecoords)
}
