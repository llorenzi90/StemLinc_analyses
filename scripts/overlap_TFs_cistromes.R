# Load data ----
source("scripts/load_functions_and_data_for_overlap_marks.R")
lwin=500
rwin=250
## TFs ----

TFs <- read.table("data/public_data/cistromes/ALL_TFs_mm39_cistrome.bed")

# In this case I do not care about strand (is this correct??)
# I want to make a window of 250 bp upstream
# of the transcript start, 100 downstream

# keep only unique TSS so we do not overload next command
TSS_id <- paste0(TSS_bed$seqid,":",
                 TSS_bed$start,
                 "-",
                 TSS_bed$end,
                 ":",
                 TSS_bed$strand)

length(unique(TSS_id))
TSS_annot <- data.frame(TSS_id=TSS_id,
                        transcript_id=TSS_bed$name,
                        gene_name=TSS_bed$score)

TSS_bed <- TSS_bed %>% filter(!duplicated(TSS_id))
TFs_unstranded_ov <- get_near_peaks_stranded(TSS_bed,TFs,lwin = lwin,rwin = rwin,sm = F)

# save this file that takes long to create:
saveRDS(TFs_unstranded_ov, file = paste0("outputs/overlap_marks/",tl_out_pref,".TFs_overlap_",lwin,"_",rwin,".rds"))
TFs_unstranded_ov$TF_rel_level=sapply(strsplit(TFs_unstranded_ov$V11,split = "_"),function(x)x[2])
length(unique(TFs_unstranded_ov$V4))
# get TSS_ids from the overlap


TSS_ids_ov <- paste0(TFs_unstranded_ov$V1,":",
                     TFs_unstranded_ov$V2,
                               "-",
                     TFs_unstranded_ov$V3,
                               ":",
                     TFs_unstranded_ov$V6)

length(unique(TSS_ids_ov))
TFs_unstranded_ov$TSS_id=TSS_ids_ov

TFs_summary <- TFs_unstranded_ov %>% group_by(V4,TF_rel_level) %>%
  summarise(TFs=paste(unique(V10),collapse = ","))

TFs_summary <- pivot_wider(TFs_summary,names_from = TF_rel_level,values_from = TFs,
                           names_prefix="TF_rel.lev_")

colnames(TFs_summary)[1] <- "transcript_id"

TFs_summary$TSS_id <- TSS_annot$TSS_id[match(TFs_summary$transcript_id,
                                             TSS_annot$transcript_id)]

TFs_summary <- ungroup(TFs_summary)
TSS_annot <- left_join(TSS_annot,TFs_unstranded_ov %>% select(TSS_id,V10,TF_rel_level))

# gene level


TFs_summary_gene_level <- TSS_annot %>% group_by(gene_name,TF_rel_level) %>%
  summarise(TFs=paste(unique(V10),collapse = ","))

TFs_summary_gene_level <- pivot_wider(TFs_summary_gene_level,names_from = TF_rel_level,values_from = TFs,
                                      names_prefix="TF_rel.lev_")
TFs_summary_gene_level <- TFs_summary_gene_level[,-ncol(TFs_summary_gene_level)]
#colnames(TFs_summary_gene_level)[1] <- "gene_name"



# write data ----
write.table(TFs_summary,paste0(outdir,tl_out_pref,".TFs_rel_level_",lwin,"_",rwin,".tsv"),
            row.names = F,quote = F,sep = "\t")

write.table(TFs_summary_gene_level,paste0(outdir,gl_out_pref,".TFs_rel_level_",lwin,"_",rwin,".tsv"),
            row.names = F,quote = F,sep = "\t")


# ChIP-seq analysis and TF motif mapping
# ChIP-seq files were downloaded from the Cistrome Data Browser
# (www.cistrome.org) (Mei et al. 2017) for 771 human TFs
# (Supplemental Table S7)—218 of which overlapped with the set
# of 519 JASPAR motifs. BEDTools (Quinlan and Hall 2010) was
# used to merge peaks for a given TF and then intersect the merged
# ChIP peaks with our set of promoters. Since Cistrome peaks were
# in hg38 and our promoters were in hg19, we first used liftOver
# (Hinrichs et al. 2006) to convert our promoters to hg38 coordi-
#   nates. Motifs were mapped in sequences using FIMO (version
# 4.11.2) (Grant et al. 2011) with a P-value threshold of 1 × 10−5.
# Motifs were assigned to ChIP-seq peaks if there was a FIMO motif
# mapped within 250 bp of the ChIP-seq peak.
