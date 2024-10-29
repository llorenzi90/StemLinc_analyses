#functions
extract_bed <- function(df,name="gene_name",score="width"){
  #cnames <- c("seqid","start","end","name","score")
  df$start <- df$start - 1

  # chromosome column
  chr_c=which(colnames(df)%in%c("seqnames","seqid","chr","chromosome","seqname","regname"))[1]
  colnames(df)[chr_c] <- "seqid"
  # name column
  if(name!="name"&"name"%in%colnames(df))df <- df[,!colnames(df)=="name"]
  colnames(df)[colnames(df)==name] <- "name"
  # score column
  # if score == "width" but there's no width col generate it
  if(score=="width"&!"width"%in%colnames(df)) df$width <- df$end - df$start
  if(score!="score"&"score"%in%colnames(df)) df <- df[,!colnames(df)=="score"]
  colnames(df)[colnames(df)==score] <- "score"


  df <- df %>% select(seqid,start,end,name,score,strand)
  chr_order=paste0("chr",c(seq(1:22),"M","X","Y"))
  df <- df%>% arrange(match(seqid, chr_order,nomatch=1000),start)
  return(df)
}

write_bed <- function(bed,file){
  write.table(bed,file,sep = "\t",row.names = F,col.names = F,quote = F)
}

get_promoters_from_bed <- function(bed, before=500, after=250 ){
  bed %>% dplyr::filter(strand %in% c("+","-"))
  bed <- bed %>% dplyr::mutate(start=ifelse(strand=="+",start - before, end - after),
                               end=start + before + after)

  # correct negative coordinates that may occur
  bed <- bed %>% mutate(start = ifelse(start<0,0,start),
                        end = ifelse(end<0,0,end))
  return(bed)
}

get_promoters_bed_from_annotated_tracking <- function(anno_track, before=500, after=250, name="V1", score="gene_name"){
  bed <- extract_bed(anno_track,name = name,score = score )
  bed <- get_promoters_from_bed(bed)
  return(bed)
}


get_unique_names <- function(names){
  rle_names=rle(names)
  unique_names <- paste0(names,"_",unlist(sapply(rle_names$lengths,seq)))
}
