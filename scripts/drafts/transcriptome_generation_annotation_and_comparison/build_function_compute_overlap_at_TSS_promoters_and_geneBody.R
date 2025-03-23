compute_overlap_at_TSS_promoters_and_geneBody <- function(TSS_bed,GTF,mark_bed)

### at TSS----
# compute intersect with TSS
ol_H3K27me3_TSS=bt.intersect(TSS_bed,
                             H3K27me3,
                             wao = T)
# the last column contains 0 if no overlap, 1 if overlap
colnames(ol_H3K27me3_TSS)[ncol(ol_H3K27me3_TSS)] <- "ol"

# summarise if any overlap at transcript_level
ol_H3K27me3_TSS <- ol_H3K27me3_TSS %>%
  group_by(V4) %>%
  summarise(TSS_ol_H3K27me3=any(ol==1))

# add gene names
ol_H3K27me3_TSS$gene_name <-
  transcript_level_data$gene_name[match(ol_H3K27me3_TSS$V4,
                                        transcript_level_data$V1)]

# summarize any overlap at gene level
ol_H3K27me3_TSS_gene_level <- ol_H3K27me3_TSS %>%
  group_by(gene_name) %>%
  summarise(TSS_ol_H3K27me3 = any(TSS_ol_H3K27me3))

# write TSS overlap
write.table(ol_H3K27me3_TSS_gene_level,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark,
                   "atTSS.tsv"),
            sep = "\t",row.names = F,quote = F)

### at promoters ----
# compute intersect with promoters
ol_H3K27me3_Prom=bt.intersect(sort_bed(proms_bed),
                              H3K27me3,wao = T)

# the last column gives the amount of overlap:
colnames(ol_H3K27me3_Prom)[ncol(ol_H3K27me3_Prom)] <- "ol"

# summarize first at transcript level, in case a single
# promoter is overlapped by more
# than one peak:
ol_H3K27me3_Prom <- ol_H3K27me3_Prom %>%
  group_by(V4) %>%
  mutate(ol=sum(ol))

# remove duplicated transcript_ids (just in case):
ol_H3K27me3_Prom <- ol_H3K27me3_Prom %>%
  filter(!duplicated(V4))

# compute fraction of overlap with promoter region length (750 by default)
ol_H3K27me3_Prom$frac_ol <-
  ol_H3K27me3_Prom$ol / (ol_H3K27me3_Prom$V3[1] - ol_H3K27me3_Prom$V2[1])

# add gene names
ol_H3K27me3_Prom$gene_name <-
  transcript_level_data$gene_name[match(ol_H3K27me3_Prom$V4,
                                        transcript_level_data$V1)]


# summarise maximum fraction
# of overlap (i.e the coverage of the
# promoter of the transcript with maximum mark coverage)
ol_H3K27me3_Prom_gene_level <- ol_H3K27me3_Prom %>%
  group_by(gene_name) %>%
  summarise(H3K27me3_promoter_cov=max(frac_ol))

# write promoter data
write.table(ol_H3K27me3_Prom_gene_level,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark
                   ,"atPromoters.tsv"),
            sep = "\t",row.names = F,quote = F)


### at gene body ----
# compute peak coverage at gene body
# this function does all the work, check source function
H3K27me3_gene_body_intersect=
  geneBody_intersect_fromGTF(GTF,H3K27me3)

# modify column names so they refer to the current mark:
colnames(H3K27me3_gene_body_intersect)[-1] <-
  paste0(mark,"_",colnames(H3K27me3_gene_body_intersect)[-1])

# write gene body data
write.table(H3K27me3_gene_body_intersect,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark,
                   "atGeneBody.tsv"),
            sep = "\t",row.names = F,quote = F)
