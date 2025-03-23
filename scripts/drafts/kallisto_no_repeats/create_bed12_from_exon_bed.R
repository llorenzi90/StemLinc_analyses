exon_bed_b12 <- exon_bed %>% arrange(seqid,start) %>% mutate(size=end-start ) %>%
  group_by(transcript_id) %>% summarise(chr=unique(seqid),
                                        Tstart=min(start),
                                        Tend=max(end),
                                        strand=strand[1],
                                        blockcount=n(),
                                        score=100,
                                        itemRGB="255,0,0",
                                        blocksizes=paste0(size,collapse = ","),
                                        blockstarts=paste0(start - Tstart,collapse = ","))

exon_bed_b12 <- exon_bed_b12 %>% select(chr,Tstart,Tend,transcript_id,score,strand,Tstart,Tend,itemRGB,blockcount,blocksizes,blockstarts)
write.table(exon_bed_b12, "outputs/bed_files/test_bed12.bed",quote = F,col.names = F,row.names = F,sep = "\t")
sum(c(160,32,66,48,71,74,102,177,1777))

#exon_bed_b12=read.table("outputs/bed_files/test_bed12.bed")
