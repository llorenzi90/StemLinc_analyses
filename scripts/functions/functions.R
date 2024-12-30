#functions
extract_bed <- function(df,name="gene_name",score="width",move_start=1){
  #cnames <- c("seqid","start","end","name","score")
  df$start <- df$start - move_start

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


  df <- df %>% dplyr::select(seqid,start,end,name,score,strand)
  chr_order=paste0("chr",c(seq(1:22),"X","Y","M"))
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


get_timestamp <- function(inchar){
  timestamp <- gsub(".*?(([0-9]{8}_[0-9]{6})(_([0-9]{3}))?).*?$", "\\2.\\4", inchar)
  return(timestamp)
}


write_info_table <- function(data_mod, data_path,modify_path=T, suffix="last"){
  fiex=tools::file_ext(data_path)
  if(!modify_path) out_path=data_path else out_path <- gsub(fiex,paste0(suffix,".",fiex),data_path)
  write.table(data_mod, out_path, quote = F, row.names = F, sep = "\t")
}


pivot_longer_from_matrix <- function(da,val2,nam2,rname="gene_name"){
  da=as.data.frame(da)
  da$gene_name=rownames(da)
  if(rname!="gene_name")colnames(da[ncol(da)])=rname
  da=pivot_longer(da,cols=1:(ncol(da)-1),
                  values_to = val2,names_to = nam2)
}

print_inclined_theme_axis_text_x <- function(){
  print(" theme(axis.text.x = element_text(angle = 45, hjust = 1)) ")}

print_inclined_theme_axis_title_x <- function(){
  print(" theme(axis.title.x = element_text(angle = 45, hjust = 1)) ")}


separate_data <- function(data,col,split=","){
  data <- as.data.frame(data)
  col2sep <- data[,col]
  seplist <- strsplit(col2sep,split = split)
  lens <- sapply(seplist,length)
  long_data <- data[rep(1:nrow(data),lens),]
  if(is.numeric(col)) col=colnames(data)[col]
  new_col=paste0(col,"_split")
  long_data[,new_col] <- unlist(seplist)
  return(long_data)
}
