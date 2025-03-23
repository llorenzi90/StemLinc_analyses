# Wikipedia: The cleavage is catalysed by the enzyme CPSF[13][18]
#and occurs 10–30 nucleotides downstream of its binding site.[19]
# This site often has the polyadenylation signal sequence AAUAAA on the RNA,
# but variants of it that bind more weakly to CPSF exist.[18][20]
# Two other proteins add specificity to the binding to an RNA: CstF and CFI.
# CstF binds to a GU-rich region further downstream of CPSF's site.[21]
# CFI recognises a third site on the RNA (a set of UGUAA sequences in mammals[22][23][24])
# and can recruit CPSF even if the AAUAAA sequence is missing.

# then polyAsite db marks the polyA sites based on 3' RNA seq, (this is the site of cleavage)
# it also incorporates polyA signals as part of the identification

# check histones,

# locate motifs for polyA signal and histone motifs (and CPSF, CstF and CFI)


# Load an process data ----

source("scripts/load_functions_and_data_for_overlap_marks.R")


polyA_path <- "data/public_data/polyAsite2.0/atlas.clusters.2.0.chr.liftOver.mm39.bed"
polyA <- read.table(polyA_path)
# colnames(polyA) <- c("chr","start","end","unique_cluster_ID","mean_tpm",
#                      "strand","fraction_samples","N_seq_protocol","mean_tpm2",
#                      "cluster_annotation",  # in order of decreasing priority: TE,
#                                             # terminal exon; EX, exonic; IN, intronic;
#                                             # DS, 1,000 nt
#                                             # downstream of an annotated terminal exon;
#                                             # AE, anti-sense to an exon; AI, anti-sense
#                                             # to an intron; AU, 1,000 nt upstream in
#                                             # anti-sense direction of a transcription
#                                             # start site; IG, intergenic)
#
#                      "polyA_signal")
polyA <- polyA %>% arrange(match(V1,chr_order),V2)

# polyA site distance----
# polyA sites should be found at the 3' (similar approach to
# what I did with CAGE)
# get polyA peaks within -100 +100 of TSS
polyA_100_100=get_near_peaks_stranded(TES_bed,polyA,lwin = 100,rwin = 100)
ncol(polyA_100_100)
# this table has the closest polyA (if present within -100 of TSS to +100) per transcript
# column V18 has the amount of overlap with TSS, 1 if overlap, negative value if offset

# Set distance to closest polyA: 0 if overlap,
#                               -X if polyA upstream (relative to transcript strand) limit: -100
#                               +X if polyA downstream (relative to transcript strand) limit: +100
polyA_100_100 <- polyA_100_100 %>%
  mutate(V18=ifelse(V18==1,0,ifelse(V8>V2,abs(V18),V18))) # set 0 if overlap, positive value if polyA after TSS, negative if before
polyA_100_100 <- polyA_100_100 %>%
  mutate(V18=ifelse(V6=="-",-V18,V18)) # set distance relative to transcript strand negative if before, positive if after
# closest whatsoever:
closest_polyA2 <- get_closest_peaks_stranded(bed = TES_bed,mark = polyA)

# summarise at gene level
closest_polyA <- polyA_100_100 %>% group_by(V5) %>%
  summarise(closest_polyA=ifelse(any(V18==0),0,V18[which.min(abs(V18))][1]))

closest_polyA2 <- closest_polyA2 %>% group_by(V5) %>%
  arrange(desc(V11)) %>% summarise(closest_polyA=V18[which.min(abs(V18))][1], #arrange by mean tpm of polyA site
                                   polyA_signal=V17[which.min(abs(V18))][1])

colnames(closest_polyA)[1] <- "gene_name"
colnames(closest_polyA2)[1] <- "gene_name"

write.table(closest_polyA,paste0(outdir,gl_out_pref,".polyA_100_100.tsv"),row.names = F,quote = F,sep = "\t")
write.table(closest_polyA2,paste0(outdir,gl_out_pref,".closest_polyA.tsv"),row.names = F,quote = F,sep = "\t")

gene_level_info <- left_join(gene_level_info,closest_polyA2)

