## Notes ---------------------------
## activation-associated histone marks, histone H3 Lys4
##   trimethylation (H3K4me3) and histone H3 Lys36 trimethylation (H3K36me3)
##   For LncHSCs, H3K4me3 was typically located at their predicted
##   transcriptional start site (TSS) and H3K36me3 along their gene bodies
##
## enhancer regions, marked by histone H3 Lys27 acetylation (H3K27ac) or
##     histone H3 Lys4 monomethylation (H3K4me1) but not H3K4me3 and H3K27me3
##   Also check Enhancer Atlas and FANTOM5
##
## 0. liftOver to mm39 (bash scripts in each directory)
## 1. Get TSS for all genes
## 2. Check peaks sizes
## 3. Overlap of TSS with either mark
## 4. Overlap of gene body with either mark
##
##

## NOTES on histone marks used ----
# Histone marks H3K4me3, H3K27me3 and H3K36me3
# for mouse hematopoietic stem cells - GSE47765

# 24 month H3K27ac Hematopoietic stem cells (HSC) GSM1544999

# 24 month H3K4me1 Hematopoietic stem cells (HSC) GSM1545000


# Setup ----
library(tidyverse)
source("scripts/source_all_functions.R")
source("scripts/load_functions_and_data_for_overlap_marks.R")
gtf_path="data/raw/LSK_StemLinc.combined.gtf"
annot_path="outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.gene_classif.tsv"
outdir="outputs/overlap_marks/histone_marks/"
dir.create(outdir)

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

# Overlap promoters with activation marks
# The below code is repetitive but carefully checked
# (should make this into a function)
# For each histone mark it generates and writes gene level tables with:
# 1- if there is overlap of any histone peak with any of the gene's TSS
# 2- the fraction of overlap of a gene's promoter with maximum overlap
# 3- the gene body coverage (fraction of exonic
#   sequence covered) by histone peaks

## H3K4me3 ----
### load data ----
mark="H3K4me3"
tfiles=list.files("data/public_data/ChIP_seq_data/",pattern = mark,full.names = T)

print(tfiles)
H3K4me3_beds=lapply(tfiles,read.table)

# check median peak length dist of each file
lapply(H3K4me3_beds,function(x)summary(x$V3 - x$V2))
H3K4me3_beds=do.call("rbind",H3K4me3_beds)
summary(H3K4me3_beds$V3 - H3K4me3_beds$V2)

# merge peaks into continous peaks
H3K4me3 <- merge_beds(H3K4me3_beds)

### at TSS----
# compute intersect with TSS
ol_H3K4me3_TSS=bt.intersect(TSS_bed,
                            H3K4me3,
                            wao = T)
# the last column contains 0 if no overlap, 1 if overlap
colnames(ol_H3K4me3_TSS)[ncol(ol_H3K4me3_TSS)] <- "ol"

# summarise if any overlap at transcript_level
ol_H3K4me3_TSS=ol_H3K4me3_TSS%>%group_by(V4)%>%
  summarise(TSS_ol_H3K4me3=any(ol==1))

# add gene names
ol_H3K4me3_TSS$gene_name <-
  transcript_level_data$gene_name[match(ol_H3K4me3_TSS$V4,
                                        transcript_level_data$V1)]

# summarize any overlap at gene level
ol_H3K4me3_TSS_gene_level <- ol_H3K4me3_TSS %>%
  group_by(gene_name) %>%
  summarise(TSS_ol_H3K4me3 = any(TSS_ol_H3K4me3))

# write TSS overlap
write.table(ol_H3K4me3_TSS_gene_level,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark,
                   "atTSS.tsv"),
            sep = "\t",row.names = F,quote = F)

### at promoters ----
# compute intersect with promoters
ol_H3K4me3_Prom=bt.intersect(sort_bed(proms_bed),
                             H3K4me3,wao = T)

# the last column gives the amount of overlap:
colnames(ol_H3K4me3_Prom)[ncol(ol_H3K4me3_Prom)] <- "ol"

# summarize first at transcript level, in case a single
# promoter is overlapped by more
# than one peak:
ol_H3K4me3_Prom <- ol_H3K4me3_Prom %>%
  group_by(V4) %>%
  mutate(ol=sum(ol))

# remove duplicated transcript_ids (just in case):
ol_H3K4me3_Prom <- ol_H3K4me3_Prom %>%
  filter(!duplicated(V4))

# compute fraction of overlap with promoter region length (750 by default)
ol_H3K4me3_Prom$frac_ol <-
  ol_H3K4me3_Prom$ol / (ol_H3K4me3_Prom$V3[1] - ol_H3K4me3_Prom$V2[1])

# add gene names
ol_H3K4me3_Prom$gene_name <-
  transcript_level_data$gene_name[match(ol_H3K4me3_Prom$V4,
                                        transcript_level_data$V1)]


# summarise maximum fraction
# of overlap (i.e the coverage of the
# promoter of the transcript with maximum mark coverage)
ol_H3K4me3_Prom_gene_level <- ol_H3K4me3_Prom %>%
  group_by(gene_name) %>%
  summarise(H3K4me3_promoter_cov=max(frac_ol))

# write promoter data
write.table(ol_H3K4me3_Prom_gene_level,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark
                   ,"atPromoters.tsv"),
            sep = "\t",row.names = F,quote = F)


### at gene body ----
# compute peak coverage at gene body
# this function does all the work, check source function
H3K4me3_gene_body_intersect=
  geneBody_intersect_fromGTF(GTF,H3K4me3)

# modify column names so they refer to the current mark:
colnames(H3K4me3_gene_body_intersect)[-1] <-
  paste0(mark,"_",colnames(H3K4me3_gene_body_intersect)[-1])

# write gene body data
write.table(H3K4me3_gene_body_intersect,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark,
                   "atGeneBody.tsv"),
            sep = "\t",row.names = F,quote = F)

# H3K27me3 ----
#   H3K27me3 is an epigenetic modification to the DNA
#    packaging protein Histone H3. It is a mark that
#    indicates the tri-methylation of lysine 27 on histone H3 protein.
#    This tri-methylation is associated with the downregulation of nearby genes via the formation of heterochromatic regions.

### load data ----
mark="H3K27me3"
tfiles=list.files("data/public_data/ChIP_seq_data/",pattern = mark,full.names = T)
print("processing files ...")
print(tfiles)

# read files
H3K27me3_beds=lapply(tfiles,read.table)

# check median peak length dist of each file
lapply(H3K27me3_beds,function(x)summary(x$V3 - x$V2))
# bind into a single data frame
H3K27me3_beds=do.call("rbind",H3K27me3_beds)
summary(H3K27me3_beds$V3 - H3K27me3_beds$V2)

# merge peaks into continous peaks
H3K27me3 <- merge_beds(H3K27me3_beds)

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

## H3K36me3 ----
### load data ----
mark="H3K36me3"
tfiles=list.files("data/public_data/ChIP_seq_data/",pattern = mark,full.names = T)
print(tfiles)
H3K36me3_beds=lapply(tfiles,read.table)

lapply(H3K36me3_beds, function(x)summary(x$V3 - x$V2))
H3K36me3_beds=do.call("rbind",H3K36me3_beds)
summary(H3K36me3_beds$V3 - H3K36me3_beds$V2)

# merge peaks into continous peaks
H3K36me3 <- merge_beds(H3K36me3_beds)

### at TSS----
# compute intersect with TSS
ol_H3K36me3_TSS=bt.intersect(TSS_bed,
                             H3K36me3,
                             wao = T)
# the last column contains 0 if no overlap, 1 if overlap
colnames(ol_H3K36me3_TSS)[ncol(ol_H3K36me3_TSS)] <- "ol"

# summarise if any overlap at transcript_level
ol_H3K36me3_TSS <- ol_H3K36me3_TSS %>%
  group_by(V4) %>%
  summarise(TSS_ol_H3K36me3=any(ol==1))

# add gene names
ol_H3K36me3_TSS$gene_name <-
  transcript_level_data$gene_name[match(ol_H3K36me3_TSS$V4,
                                        transcript_level_data$V1)]

# summarize any overlap at gene level
ol_H3K36me3_TSS_gene_level <- ol_H3K36me3_TSS %>%
  group_by(gene_name) %>%
  summarise(TSS_ol_H3K36me3 = any(TSS_ol_H3K36me3))

# write TSS overlap
write.table(ol_H3K36me3_TSS_gene_level,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark,
                   "atTSS.tsv"),
            sep = "\t",row.names = F,quote = F)

### at promoters ----
# compute intersect with promoters
ol_H3K36me3_Prom=bt.intersect(sort_bed(proms_bed),
                              H3K36me3,wao = T)

# the last column gives the amount of overlap:
colnames(ol_H3K36me3_Prom)[ncol(ol_H3K36me3_Prom)] <- "ol"

# summarize first at transcript level, in case a single
# promoter is overlapped by more
# than one peak:
ol_H3K36me3_Prom <- ol_H3K36me3_Prom %>%
  group_by(V4) %>%
  mutate(ol=sum(ol))

# remove duplicated transcript_ids (just in case):
ol_H3K36me3_Prom <- ol_H3K36me3_Prom %>%
  filter(!duplicated(V4))

# compute fraction of overlap with promoter region length (750 by default)
ol_H3K36me3_Prom$frac_ol <-
  ol_H3K36me3_Prom$ol / (ol_H3K36me3_Prom$V3[1] - ol_H3K36me3_Prom$V2[1])

# add gene names
ol_H3K36me3_Prom$gene_name <-
  transcript_level_data$gene_name[match(ol_H3K36me3_Prom$V4,
                                        transcript_level_data$V1)]


# summarise maximum fraction
# of overlap (i.e the coverage of the
# promoter of the transcript with maximum mark coverage)
ol_H3K36me3_Prom_gene_level <- ol_H3K36me3_Prom %>%
  group_by(gene_name) %>%
  summarise(H3K36me3_promoter_cov=max(frac_ol))

# write promoter data
write.table(ol_H3K36me3_Prom_gene_level,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark
                   ,"atPromoters.tsv"),
            sep = "\t",row.names = F,quote = F)


### at gene body ----
# compute peak coverage at gene body
# this function does all the work, check source function
H3K36me3_gene_body_intersect=
  geneBody_intersect_fromGTF(GTF,H3K36me3)

# modify column names so they refer to the current mark:
colnames(H3K36me3_gene_body_intersect)[-1] <-
  paste0(mark,"_",colnames(H3K36me3_gene_body_intersect)[-1])

# write gene body data
write.table(H3K36me3_gene_body_intersect,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark,
                   "atGeneBody.tsv"),
            sep = "\t",row.names = F,quote = F)


# H3K27ac ----
### load data ----
H3K27ac=read.table("data/public_data/ChIP_seq_data/GSM1544999_m24_H3K27ac.narrowPeak.mm39.bed")
mark="H3K27ac"

### at TSS----
# compute intersect with TSS
ol_H3K27ac_TSS=bt.intersect(TSS_bed,
                            H3K27ac,
                            wao = T)
# the last column contains 0 if no overlap, 1 if overlap
colnames(ol_H3K27ac_TSS)[ncol(ol_H3K27ac_TSS)] <- "ol"

# summarise if any overlap at transcript_level
ol_H3K27ac_TSS <- ol_H3K27ac_TSS %>%
  group_by(V4) %>%
  summarise(TSS_ol_H3K27ac=any(ol==1))

# add gene names
ol_H3K27ac_TSS$gene_name <-
  transcript_level_data$gene_name[match(ol_H3K27ac_TSS$V4,
                                        transcript_level_data$V1)]

# summarize any overlap at gene level
ol_H3K27ac_TSS_gene_level <- ol_H3K27ac_TSS %>%
  group_by(gene_name) %>%
  summarise(TSS_ol_H3K27ac = any(TSS_ol_H3K27ac))

# write TSS overlap
write.table(ol_H3K27ac_TSS_gene_level,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark,
                   "atTSS.tsv"),
            sep = "\t",row.names = F,quote = F)

### at promoters ----
# compute intersect with promoters
ol_H3K27ac_Prom=bt.intersect(sort_bed(proms_bed),
                             H3K27ac,wao = T)

# the last column gives the amount of overlap:
colnames(ol_H3K27ac_Prom)[ncol(ol_H3K27ac_Prom)] <- "ol"

# summarize first at transcript level, in case a single
# promoter is overlapped by more
# than one peak:
ol_H3K27ac_Prom <- ol_H3K27ac_Prom %>%
  group_by(V4) %>%
  mutate(ol=sum(ol))

# remove duplicated transcript_ids (just in case):
ol_H3K27ac_Prom <- ol_H3K27ac_Prom %>%
  filter(!duplicated(V4))

# compute fraction of overlap with promoter region length (750 by default)
ol_H3K27ac_Prom$frac_ol <-
  ol_H3K27ac_Prom$ol / (ol_H3K27ac_Prom$V3[1] - ol_H3K27ac_Prom$V2[1])

# add gene names
ol_H3K27ac_Prom$gene_name <-
  transcript_level_data$gene_name[match(ol_H3K27ac_Prom$V4,
                                        transcript_level_data$V1)]


# summarise maximum fraction
# of overlap (i.e the coverage of the
# promoter of the transcript with maximum mark coverage)
ol_H3K27ac_Prom_gene_level <- ol_H3K27ac_Prom %>%
  group_by(gene_name) %>%
  summarise(H3K27ac_promoter_cov=max(frac_ol))

# write promoter data
write.table(ol_H3K27ac_Prom_gene_level,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark
                   ,"atPromoters.tsv"),
            sep = "\t",row.names = F,quote = F)


### at gene body ----
# compute peak coverage at gene body
# this function does all the work, check source function
H3K27ac_gene_body_intersect=
  geneBody_intersect_fromGTF(GTF,H3K27ac)

# modify column names so they refer to the current mark:
colnames(H3K27ac_gene_body_intersect)[-1] <-
  paste0(mark,"_",colnames(H3K27ac_gene_body_intersect)[-1])

# write gene body data
write.table(H3K27ac_gene_body_intersect,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark,
                   "atGeneBody.tsv"),
            sep = "\t",row.names = F,quote = F)






# H3K4me1 ----
### load data ----
H3K4me1=read.table("data/public_data/ChIP_seq_data/GSM1545000_m24_H3K4me1.narrowPeak.mm39.bed")
mark="H3K4me1"

### at TSS----
# compute intersect with TSS
ol_H3K4me1_TSS=bt.intersect(TSS_bed,
                            H3K4me1,
                            wao = T)
# the last column contains 0 if no overlap, 1 if overlap
colnames(ol_H3K4me1_TSS)[ncol(ol_H3K4me1_TSS)] <- "ol"

# summarise if any overlap at transcript_level
ol_H3K4me1_TSS <- ol_H3K4me1_TSS %>%
  group_by(V4) %>%
  summarise(TSS_ol_H3K4me1=any(ol==1))

# add gene names
ol_H3K4me1_TSS$gene_name <-
  transcript_level_data$gene_name[match(ol_H3K4me1_TSS$V4,
                                        transcript_level_data$V1)]

# summarize any overlap at gene level
ol_H3K4me1_TSS_gene_level <- ol_H3K4me1_TSS %>%
  group_by(gene_name) %>%
  summarise(TSS_ol_H3K4me1 = any(TSS_ol_H3K4me1))

# write TSS overlap
write.table(ol_H3K4me1_TSS_gene_level,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark,
                   "atTSS.tsv"),
            sep = "\t",row.names = F,quote = F)

### at promoters ----
# compute intersect with promoters
ol_H3K4me1_Prom=bt.intersect(sort_bed(proms_bed),
                             H3K4me1,wao = T)

# the last column gives the amount of overlap:
colnames(ol_H3K4me1_Prom)[ncol(ol_H3K4me1_Prom)] <- "ol"

# summarize first at transcript level, in case a single
# promoter is overlapped by more
# than one peak:
ol_H3K4me1_Prom <- ol_H3K4me1_Prom %>%
  group_by(V4) %>%
  mutate(ol=sum(ol))

# remove duplicated transcript_ids (just in case):
ol_H3K4me1_Prom <- ol_H3K4me1_Prom %>%
  filter(!duplicated(V4))

# compute fraction of overlap with promoter region length (750 by default)
ol_H3K4me1_Prom$frac_ol <-
  ol_H3K4me1_Prom$ol / (ol_H3K4me1_Prom$V3[1] - ol_H3K4me1_Prom$V2[1])

# add gene names
ol_H3K4me1_Prom$gene_name <-
  transcript_level_data$gene_name[match(ol_H3K4me1_Prom$V4,
                                        transcript_level_data$V1)]


# summarise maximum fraction
# of overlap (i.e the coverage of the
# promoter of the transcript with maximum mark coverage)
ol_H3K4me1_Prom_gene_level <- ol_H3K4me1_Prom %>%
  group_by(gene_name) %>%
  summarise(H3K4me1_promoter_cov=max(frac_ol))

# write promoter data
write.table(ol_H3K4me1_Prom_gene_level,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark
                   ,"atPromoters.tsv"),
            sep = "\t",row.names = F,quote = F)


### at gene body ----
# compute peak coverage at gene body
# this function does all the work, check source function
H3K4me1_gene_body_intersect=
  geneBody_intersect_fromGTF(GTF,H3K4me1)

# modify column names so they refer to the current mark:
colnames(H3K4me1_gene_body_intersect)[-1] <-
  paste0(mark,"_",colnames(H3K4me1_gene_body_intersect)[-1])

# write gene body data
write.table(H3K4me1_gene_body_intersect,
            paste0(outdir,
                   gl_out_pref,
                   ".",mark,
                   "atGeneBody.tsv"),
            sep = "\t",row.names = F,quote = F)



