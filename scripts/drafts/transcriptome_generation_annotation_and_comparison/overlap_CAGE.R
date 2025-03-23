# Load an process data ----
source("scripts/load_functions_and_data_for_overlap_marks.R")

CAGE=read.table("data/public_data/CAGE/TSS_mouse.liftOver.mm39.bed")
chr_order=paste0("chr",c(seq(1:22),"M","X","Y"))
CAGE <- CAGE %>% arrange(match(V1,chr_order),V2)

# CAGE ----
# get CAGE peaks within -500 +250 of TSS
CAGE_500_250=get_near_peaks_stranded(TSS_bed,CAGE,lwin = 500,rwin = 250)
# this table has the closest CAGE (if present within -500 of TSS to +250) per transcript
# column V16 has the amount of overlap with TSS, 1 if overlap, negative value if offset

# Set distance to closest CAGE: 0 if overlap,
#                               -X if CAGE upstream (relative to transcript strand) limit: -500
#                               +X if CAGE downstream (relative to transcript strand) limit: +250
CAGE_500_250 <- CAGE_500_250 %>% mutate(V16=ifelse(V16==1,0,ifelse(V8>V2,abs(V16),V16))) # set 0 if overlap, positive value if CAGE after TSS, negative if before
CAGE_500_250 <- CAGE_500_250 %>% mutate(V16=ifelse(V6=="-",-V16,V16)) # set distance relative to transcript strand negative if before, positive if after
# closest whatsoever:
closest_CAGE2 <- get_closest_peaks_stranded(bed = TSS_bed,mark = CAGE)

# summarise at gene level
closest_CAGE <- CAGE_500_250 %>% group_by(V5) %>%
  summarise(closest_CAGE=ifelse(any(V16==0),0,V16[which.min(abs(V16))][1]))

closest_CAGE2 <- closest_CAGE2 %>% group_by(V5) %>%
  arrange(V16) %>% summarise(closest_CAGE=V16[which.min(abs(V16))][1])

colnames(closest_CAGE)[1] <- "gene_name"
colnames(closest_CAGE2)[1] <- "gene_name"

write.table(closest_CAGE,paste0(outdir,gl_out_pref,".CAGE_500_250.tsv"),row.names = F,quote = F,sep = "\t")
write.table(closest_CAGE2,paste0(outdir,gl_out_pref,".closest_CAGE.tsv"),row.names = F,quote = F,sep = "\t")

gene_level_info <- left_join(gene_level_info,closest_CAGE2)

