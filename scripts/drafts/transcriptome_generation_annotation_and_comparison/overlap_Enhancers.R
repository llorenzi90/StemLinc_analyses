source("scripts/load_functions_and_data_for_overlap_marks.R")
outdir="outputs/overlap_marks/Enhancers/"
#Enhancers ----
## Enhancer Atlas ----
Enhancer_atlas=read.table("data/public_data/Enhancer_data/Enhancer_atlas_LSK.mm39.bed")
# in this case I will use a symmetric window of 500 bp

# closest whatsoever:
closest_Enhancer_Atlas <- get_closest_peaks_stranded(bed = TSS_bed,
                                                     mark = Enhancer_atlas,s = F)

# summarise at gene level


distance_to_closest_enhancer <- closest_Enhancer_Atlas %>% group_by(V5) %>% arrange(V11) %>%
  summarise(closest_Enhancer_Atlas=V11[which.min(abs(V11))][1])

colnames(distance_to_closest_enhancer)[1] <- "gene_name"

write.table(distance_to_closest_enhancer,
            paste0(outdir,gl_out_pref,".closest_Enhancer_Atlas.tsv"),row.names = F,quote = F,sep = "\t")


## FANTOM enhancers ----
FANTOM_enhancers <- read.table("data/public_data/Enhancer_data/mouse_permissive_enhancers_phase_1_and_2.mm39.bed")

# closest whatsoever:
closest_FANTOM_enhancers <- get_closest_peaks_stranded(bed = TSS_bed,
                                                     mark = FANTOM_enhancers,s = F)


# summarise at gene level

distance_to_closest_enhancer <- closest_FANTOM_enhancers %>% group_by(V5) %>% arrange(V19) %>%
  summarise(closest_FANTOM_enhancer=V19[which.min(abs(V19))][1])

colnames(distance_to_closest_enhancer)[1] <- "gene_name"

write.table(distance_to_closest_enhancer,paste0(outdir,gl_out_pref,".closest_FANTOM_enhancers.tsv"),row.names = F,quote = F,sep = "\t")


# approach similar to the overlap with histones ----
gtf_path="data/raw/LSK_StemLinc.combined.gtf"
annot_path="outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.gene_classif.tsv"
outdir="outputs/overlap_marks/Enhancers/"
#dir.create(outdir)

# read raw sample gtf
GTF=readGFF(gtf_path)
# read filtered annot to select filtered transcripts:
annot=read.table(annot_path,header = T)

# filter GTF to retain only transcripts in the filtered set
# and assign gene_names to gene_id from the annot file
GTF <- GTF%>%filter(transcript_id%in%annot$V1) %>%
  mutate(gene_id=annot$gene_name[match(transcript_id,
                                       annot$V1)])

# generate a bed with transcript coordinates
tr_bed <- extract_bed(GTF%>%filter(type=="transcript"),
                      name = "transcript_id")

# extract promoter regions from this
proms_bed <- get_promoters_from_bed(tr_bed,
                                    before = 500,
                                    after = 250)

## Enhancer Atlas ----
mark="Enhancer_Atlas"
### at TSS----
# compute intersect with TSS
Enhancer_Atlas_TSS=bt.intersect(TSS_bed,
                            Enhancer_atlas,
                            wao = T)
# the last column contains 0 if no overlap, 1 if overlap
colnames(Enhancer_Atlas_TSS)[ncol(Enhancer_Atlas_TSS)] <- "ol"

# summarise if any overlap at transcript_level
ol_Enhancer_Atlas_TSS=Enhancer_Atlas_TSS%>%group_by(V4)%>%
  summarise(TSS_ol_Enhancer_Atlas=any(ol==1))

# add gene names
ol_Enhancer_Atlas_TSS$gene_name <-
  transcript_level_data$gene_name[match(ol_Enhancer_Atlas_TSS$V4,
                                        transcript_level_data$V1)]

# summarize any overlap at gene level
ol_Enhancer_Atlas_TSS_gene_level <- ol_Enhancer_Atlas_TSS %>%
  group_by(gene_name) %>%
  summarise(TSS_ol_Enhancer_Atlas = any(TSS_ol_Enhancer_Atlas))

# write TSS overlap
write.table(ol_Enhancer_Atlas_TSS_gene_level,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark,
                   "atTSS.tsv"),
            sep = "\t",row.names = F,quote = F)

### at promoters ----
# compute intersect with promoters
ol_Enhancer_Atlas_Prom=bt.intersect(sort_bed(proms_bed),
                             Enhancer_atlas,wao = T)

# the last column gives the amount of overlap:
colnames(ol_Enhancer_Atlas_Prom)[ncol(ol_Enhancer_Atlas_Prom)] <- "ol"

# summarize first at transcript level, in case a single
# promoter is overlapped by more
# than one peak:
ol_Enhancer_Atlas_Prom <- ol_Enhancer_Atlas_Prom %>%
  group_by(V4) %>%
  mutate(ol=sum(ol))

# remove duplicated transcript_ids (just in case):
ol_Enhancer_Atlas_Prom <- ol_Enhancer_Atlas_Prom %>%
  filter(!duplicated(V4))

# compute fraction of overlap with promoter region length (750 by default)
ol_Enhancer_Atlas_Prom$frac_ol <-
  ol_Enhancer_Atlas_Prom$ol / (ol_Enhancer_Atlas_Prom$V3[1] - ol_Enhancer_Atlas_Prom$V2[1])

# add gene names
ol_Enhancer_Atlas_Prom$gene_name <-
  transcript_level_data$gene_name[match(ol_Enhancer_Atlas_Prom$V4,
                                        transcript_level_data$V1)]


# summarise maximum fraction
# of overlap (i.e the coverage of the
# promoter of the transcript with maximum mark coverage)
ol_Enhancer_Atlas_Prom_gene_level <- ol_Enhancer_Atlas_Prom %>%
  group_by(gene_name) %>%
  summarise(Enhancer_Atlas_promoter_cov=max(frac_ol))

# write promoter data
write.table(ol_Enhancer_Atlas_Prom_gene_level,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark
                   ,"atPromoters.tsv"),
            sep = "\t",row.names = F,quote = F)


### at gene body ----
# compute peak coverage at gene body
# this function does all the work, check source function
Enhancer_Atlas_gene_body_intersect=
  geneBody_intersect_fromGTF(GTF,Enhancer_atlas)

# modify column names so they refer to the current mark:
colnames(Enhancer_Atlas_gene_body_intersect)[-1] <-
  paste0(mark,"_",colnames(Enhancer_Atlas_gene_body_intersect)[-1])

# write gene body data
write.table(Enhancer_Atlas_gene_body_intersect,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark,
                   "atGeneBody.tsv"),
            sep = "\t",row.names = F,quote = F)



## FANTOM enhancers ----
mark="FANTOM_enhancers"
### at TSS----
# compute intersect with TSS
FANTOM_enhancers_TSS=bt.intersect(TSS_bed,
                                FANTOM_enhancers,
                                wao = T)
# the last column contains 0 if no overlap, 1 if overlap
colnames(FANTOM_enhancers_TSS)[ncol(FANTOM_enhancers_TSS)] <- "ol"

# summarise if any overlap at transcript_level
ol_FANTOM_enhancers_TSS=FANTOM_enhancers_TSS%>%group_by(V4)%>%
  summarise(TSS_ol_FANTOM_enhancers=any(ol==1))

# add gene names
ol_FANTOM_enhancers_TSS$gene_name <-
  transcript_level_data$gene_name[match(ol_FANTOM_enhancers_TSS$V4,
                                        transcript_level_data$V1)]

# summarize any overlap at gene level
ol_FANTOM_enhancers_TSS_gene_level <- ol_FANTOM_enhancers_TSS %>%
  group_by(gene_name) %>%
  summarise(TSS_ol_FANTOM_enhancers = any(TSS_ol_FANTOM_enhancers))

# write TSS overlap
write.table(ol_FANTOM_enhancers_TSS_gene_level,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark,
                   "atTSS.tsv"),
            sep = "\t",row.names = F,quote = F)

### at promoters ----
# compute intersect with promoters
ol_FANTOM_enhancers_Prom=bt.intersect(sort_bed(proms_bed),
                                    FANTOM_enhancers,wao = T)

# the last column gives the amount of overlap:
colnames(ol_FANTOM_enhancers_Prom)[ncol(ol_FANTOM_enhancers_Prom)] <- "ol"

# summarize first at transcript level, in case a single
# promoter is overlapped by more
# than one peak:
ol_FANTOM_enhancers_Prom <- ol_FANTOM_enhancers_Prom %>%
  group_by(V4) %>%
  mutate(ol=sum(ol))

# remove duplicated transcript_ids (just in case):
ol_FANTOM_enhancers_Prom <- ol_FANTOM_enhancers_Prom %>%
  filter(!duplicated(V4))

# compute fraction of overlap with promoter region length (750 by default)
ol_FANTOM_enhancers_Prom$frac_ol <-
  ol_FANTOM_enhancers_Prom$ol / (ol_FANTOM_enhancers_Prom$V3[1] - ol_FANTOM_enhancers_Prom$V2[1])

# add gene names
ol_FANTOM_enhancers_Prom$gene_name <-
  transcript_level_data$gene_name[match(ol_FANTOM_enhancers_Prom$V4,
                                        transcript_level_data$V1)]


# summarise maximum fraction
# of overlap (i.e the coverage of the
# promoter of the transcript with maximum mark coverage)
ol_FANTOM_enhancers_Prom_gene_level <- ol_FANTOM_enhancers_Prom %>%
  group_by(gene_name) %>%
  summarise(FANTOM_enhancers_promoter_cov=max(frac_ol))

# write promoter data
write.table(ol_FANTOM_enhancers_Prom_gene_level,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark
                   ,"atPromoters.tsv"),
            sep = "\t",row.names = F,quote = F)


### at gene body ----
# compute peak coverage at gene body
# this function does all the work, check source function
FANTOM_enhancers_gene_body_intersect=
  geneBody_intersect_fromGTF(GTF,FANTOM_enhancers)

# modify column names so they refer to the current mark:
colnames(FANTOM_enhancers_gene_body_intersect)[-1] <-
  paste0(mark,"_",colnames(FANTOM_enhancers_gene_body_intersect)[-1])

# write gene body data
write.table(FANTOM_enhancers_gene_body_intersect,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark,
                   "atGeneBody.tsv"),
            sep = "\t",row.names = F,quote = F)
