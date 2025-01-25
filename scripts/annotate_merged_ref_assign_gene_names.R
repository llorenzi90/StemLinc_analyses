## Head -------------------------------------
##
##
## Purpose of script: annotate merged refs transcriptome
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-07-18 / updated and extended 2025-01-21
##
## Email: lucialorenzi90@gmail.com
##
## Notes ---------------------------
##
## As it is the reference transcriptome I am using,
## I want to annotate it and clean it up a bit, in order to
## identify and keep track of possible conflicts in biotypes
## (e.g. merged genes of different biotypes or loci that become too long)
##
## I'll aggregate by gene name, GENCODE gene name first
##
## Input files: - GTF of merged GENCODE vM31 and RefSeq - ANNOTATION RELEASE DATE:	10-April-2023
##              merging was performed with gffcompare
##              - individual GTFs
##              - gffcompare tracking file
##

## Setup ---------------------------

options(scipen = 999)
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
require(GenomicRanges)
require(rtracklayer)
library(trastools)
## Load data---------------------------
#reference datasets
refGTF=readGFF('/home/llorenzi/Documentos/references/merged_refs.combined.gtf')
gencode=readGFF('/home/llorenzi/Documentos/references/gencode.vM31.primary_assembly.annotation.gtf')
refseq=readGFF('/home/llorenzi/Documentos/references/GCF_000001635.27_GRCm39_genomic.nocmm.chr.gff')
merged_refs_track=read.table('/home/llorenzi/Documentos/references/merged_refs.tracking')

# Annotate tracking file ----

# in terms of gene ids, gene name, transcript ids, gene biotype in each
# reference transcriptome

# each row in the tracking file correspond to a transcript present either in
# gencode, refseq or both.
# gffcompare assigns new transcript and gene ids to these transcripts,
# but from now on
# I don't want to use them because it can be misleading
# (e.g. mixing different overlapping genes together)
# I'll aggregate by gene name, which could also have some drawbacks (redundance),
# but I'll be aware of them and keep track of conflicts

# First check that all transcript ids in tracking file correspond to
# transcript ids in the merged gtf
table(merged_refs_track$V1==refGTF$transcript_id[refGTF$type=="transcript"])
table(merged_refs_track$V1%in%refGTF$transcript_id)
# they all correspond, although the order is different

# Split the info stored in q1 and q2 columns in merged refs tracking to annotate the
# transcripts in the merged transcriptome back to their source
# for this use custom function "split_samples_info"
merged_ref_annot=trastools::split_samples_info(merged_refs_track,qnames = c("gencode","refseq"))

# Check that gencode genes correspond all to gencode gene_ids in gencode gtf
table(unique(merged_ref_annot$gencode_gene)%in%gencode$gene_id)

### Parse RefSeq gene ids ----
# As RefSeq annotation is a bit more complex,
# RefSeq gene ids in the tracking file are a combination of RefSeq gene and Parent
table(unique(merged_ref_annot$refseq_gene)%in%c(refseq$gene,
                                               unlist(refseq$Parent)))

# check what are the ones that are "Parent"
rsg=unique(merged_ref_annot$refseq_gene)[!unique(merged_ref_annot$refseq_gene)%in%refseq$gene]
head(rsg,50)
table(grepl("gene-",rsg))
grep("gene-",rsg,invert = T,value = T)
table(grepl("rna-",grep("gene-",rsg,invert = T,value = T)))
# most of them are "gene-" and others are "rna-"

# All refseq transcript ids do match terms in the
# ID column in refseq transcriptome:
table(unique(merged_ref_annot$refseq_transcript)%in%refseq$ID)

# extract subset of refseq gtf with assigned gene biotype
refseq_biot=refseq[!is.na(refseq$gene_biotype),]
refseq_biot=as.data.frame(refseq_biot)
length(unique(refseq_biot$gene))
length(unique(refseq$gene))
table(merged_ref_annot$refseq_transcript%in%refseq_biot$ID)

# if I can match every RefSeq transcript with a RefSeq gene, then
# I have a biotype for each one
table(is.na(refseq$gene[match(merged_ref_annot$refseq_transcript,refseq$ID)]))

table(sapply(strsplit(rsg,split = "-"),function(sp)paste0(sp[-1],
                                                          collapse = "-"))%in%
        refseq$gene)
# If I remove the leading xxx- to the Parent ids, a subset of them
# now match refseq genes
# 1233 gene ids, still do not match:
tt=rsg[!sapply(strsplit(rsg,split = "-"),function(sp)paste0(sp[-1],collapse = "-"))%in%refseq$gene]

# I can get their corresponding gene ids by matching their refseq_transcript in
# tracking file with their refseq ID in refseq gtf
tt_IDs=merged_ref_annot$refseq_transcript[match(tt,merged_ref_annot$refseq_gene)]
anyNA(refseq$gene[refseq$ID%in%tt_IDs])
# all these IDs do have an assigned refseq gene

# remove the leading "gene-" for the ones that do not have a refseq gene
tcond=!merged_ref_annot$refseq_gene%in%
  refseq$gene&!merged_ref_annot$refseq_gene=="-"
merged_ref_annot$refseq_gene[tcond]  <-
  sapply(strsplit(merged_ref_annot$refseq_gene[tcond],split = "-"),function(sp)paste0(sp[-1],collapse = "-"))

table(unique(merged_ref_annot$refseq_gene)%in%refseq$gene)

# there are still the 1233 ones that do not match a refseq gene
# for these, use the refseq ID to get the gene:
tcond=!merged_ref_annot$refseq_gene%in%refseq$gene&!merged_ref_annot$refseq_gene=="-"
merged_ref_annot$refseq_gene[tcond] <-
  refseq$gene[match(merged_ref_annot$refseq_transcript[tcond],
                    refseq$ID)]

# Now all refseq_gene in merged_ref_annot should match a refseq gene
table(unique(merged_ref_annot$refseq_gene)%in%refseq$gene)
unique(merged_ref_annot$refseq_gene)[!unique(merged_ref_annot$refseq_gene)%in%refseq$gene] #as expected this is "-"
table(unique(merged_ref_annot$refseq_gene)%in%refseq_biot$gene) # now all refseq genes have a correspoinding  biotype

### Assign biotypes from both annotations ----
merged_ref_annot$gencode_trbiotype=gencode$transcript_type[match(merged_ref_annot$gencode_transcript,
                                                                gencode$transcript_id)]

merged_ref_annot$gencode_biotype=gencode$gene_type[match(merged_ref_annot$gencode_gene,
                                                        gencode$gene_id)]

merged_ref_annot$refseq_biotype=refseq_biot$gene_biotype[match(merged_ref_annot$refseq_gene,
                                                              refseq_biot$gene)]

table(merged_ref_annot$gencode_biotype,merged_ref_annot$refseq_biotype)

# Summarise at gene level ----
# General rule is to take the gencode gene name first (in case of differences between both annotations)

# Add gencode gene names
gencode_tr=gencode%>%filter(type=="transcript")

# add gencode gene names to tracking file (so far it has only gencode gene_id)
merged_ref_annot$gencode_gene_name=gencode_tr$gene_name[match(merged_ref_annot$gencode_transcript,
                                                             gencode_tr$transcript_id)]

table(merged_ref_annot$refseq_gene!="-"&merged_ref_annot$gencode_gene!="-")

table(merged_ref_annot$refseq_gene!="-",merged_ref_annot$gencode_gene!="-")
#        FALSE   TRUE
# FALSE      0 104309
# TRUE   96231  44046
# there are only 44046 transcripts (18%) that are in both ref transcriptomes!

table(merged_ref_annot$gencode_gene_name==merged_ref_annot$refseq_gene)
#FALSE   TRUE
#106093  42262
# from these 96% have matching gene names

tr_gene_no_match=merged_ref_annot$V1[!is.na(merged_ref_annot$gencode_biotype)&
                                      !is.na(merged_ref_annot$refseq_biotype)&
                                      merged_ref_annot$gencode_gene_name!=merged_ref_annot$refseq_gene]

a_nomatch=merged_ref_annot%>%filter(V1%in%tr_gene_no_match)
length(unique(a_nomatch$refseq_gene))

# I'll summarise at gene level

# I'll use the gencode gene names first:
merged_ref_annot$gene_name=merged_ref_annot$gencode_gene_name
# When not in gencode, use refseq gene name,

merged_ref_annot$gene_name[is.na(merged_ref_annot$gene_name)]=merged_ref_annot$refseq_gene[is.na(merged_ref_annot$gene_name)]
length(unique(merged_ref_annot$gene_name))
length(unique(gencode$gene_name))
length(unique(gencode$gene_id))
# check if any gene spans more than one chr, many biotypes
merged_ref_annot$chr=refGTF$seqid[match(merged_ref_annot$V1,
                                       refGTF$transcript_id)]


# here starts: /run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/scripts/StemLinc_code/inciso_parse_merged_Ref_annot.R

#merged_ref_annot$gene_name=merged_ref_annot$gencode_gene_name
table(is.na(merged_ref_annot$gene_name))

# how to assign gene names to gffcompare loci/transcripts
# cases:
# gffcmp.gene_id:
#  a) contains transcripts from gencode alone
#  b) contains transcripts from refseq alone (check if gene name is also in gencode)
#  c) contains transcripts from both and:
#     c.1) their gene_names are unique after removing NAs (= gene match between gencode and refseq, no matter if only diff isoforms,
#            biotypes could also differ, we later mark it as conflict )
#     c.2) their gene_names are different because they have different gene synonyms that correspond to the same gene (indirect
#         match, can we check this with refseq synonyms)
#     c.3) different overlapping genes (same or different biotype, biotypes also can differ)
#         c.3.1) refseq has only one biotype and gencode has only another, use gencode biotype
#         c.3.2) refseq has two or more biotypes and gencode has only one and matches one gene name (or viceversa), decision: refseq has extra genes in other biotypes
#     c.3) gene names between ref are not matching in a matching transcript, then keep one gencode gene name
#
# special cases:
#  A) genes that have more than one genomic location (there are few) --> need to modify their ids
#  B) genes that, in a given reference have overlapping genes from the same biotype --> split


# aggregate unique gencode gene name and unique refseq gene name
a <- character(0)
identical(a,refseq$gene_synonym[[1]])
ischar0=sapply(refseq$gene_synonym,function(s)identical(a,s))

synonyms=refseq[!ischar0,c("gene","Name","gene_synonym")]
table(synonyms$gene==synonyms$Name)
synonyms_ul=unlist(synonyms$gene_synonym)
synonymsdf=as.data.frame(synonyms)
synonymsdf=synonymsdf[rep(1:nrow(synonymsdf),sapply(synonyms$gene_synonym,length)),]
synonymsdf$gene_synonym=synonyms_ul

merged_ref_annot$gencode_gene_name[is.na(merged_ref_annot$gencode_gene_name)]="-"

# For each transfrag in the merged refs transcriptome
# get the combination of original gencode gene name and RefSeq gene name
merged_ref_annot$merged_gene_names=paste0(merged_ref_annot$gencode_gene_name,";",merged_ref_annot$refseq_gene)

# For each gffcompare XLOC id (merged loci)
# get all merged_gene_names
unique_gencode_refseq_names_per_XLOC=merged_ref_annot%>%group_by(V2)%>% reframe(gencode_refseq=unique(merged_gene_names))

# select XLOC that have more than one combination:
dup_XLOC=unique(unique_gencode_refseq_names_per_XLOC$V2[duplicated(unique_gencode_refseq_names_per_XLOC$V2)])
length(unique(dup_XLOC))
#[1] 25089
# out of 62449 unique XLOC loci (40%) have more than one combination

### simple gff loci ----
# Let's start by the simplest case,
# keep only genes with only one combination (easiest case)
# These are 37360 XLOCs

non_dup_XLOC=unique_gencode_refseq_names_per_XLOC%>%filter(!V2%in%dup_XLOC)
non_dup_XLOC=separate(non_dup_XLOC,col=2,sep = ";",into = c("gencode","refseq"))
non_dup_XLOC$type_gene=ifelse(non_dup_XLOC$gencode==non_dup_XLOC$refseq,"match",
                              ifelse(non_dup_XLOC$refseq=="-","gencode",
                                     ifelse(non_dup_XLOC$gencode=="-","refseq","no_match")))
non_dup_XLOC_nomatch=non_dup_XLOC%>%filter(type_gene=="no_match")

synonym_match=mapply(function(x,y) any(synonymsdf$gene_synonym[synonymsdf$gene==y]==x),
                     non_dup_XLOC_nomatch$gencode,non_dup_XLOC_nomatch$refseq)

non_dup_XLOC_nomatch$synonym_match=synonym_match

non_dup_XLOC$type_gene[non_dup_XLOC$V2%in%non_dup_XLOC_nomatch$V2[non_dup_XLOC_nomatch$synonym_match]]="synonym_match"

table(non_dup_XLOC$type_gene)

# For this set of "no_match" genes, I actually do not have to do anything, because
# my first rule applies here: if differing names, use gencode name

# I do have to check if the refseq gene names that do not match any gencode name
table(non_dup_XLOC$refseq[non_dup_XLOC$type_gene=="refseq"]%in%merged_ref_annot$gencode_gene_name)
# only 13
refseq_dup_13=non_dup_XLOC$refseq[non_dup_XLOC$type_gene=="refseq"][non_dup_XLOC$refseq[non_dup_XLOC$type_gene=="refseq"]%in%merged_ref_annot$gencode_gene_name]

View(refGTF%>%filter(gene_name%in%refseq_dup_13)%>% select(seqid,start,end,gene_id,gene_name)%>%arrange(gene_name))

duplicated_gencode_gene_in_nondup_XLOC=unique(non_dup_XLOC$gencode[duplicated(non_dup_XLOC$gencode)])
duplicated_gencode_gene_in_nondup_XLOC=duplicated_gencode_gene_in_nondup_XLOC[-1]
View(non_dup_XLOC%>%filter(gencode%in%duplicated_gencode_gene_in_nondup_XLOC))

# all of them are contiguous genes, therefore I can collapse them with the
# gencode transcripts, so I do not have to do anything extra so far

# add biotypes
non_dup_XLOC$gencode_biotype=merged_ref_annot$gencode_biotype[match(non_dup_XLOC$gencode,
                                                                    merged_ref_annot$gencode_gene_name)]

non_dup_XLOC$refseq_biotype=merged_ref_annot$refseq_biotype[match(non_dup_XLOC$refseq,
                                                                  merged_ref_annot$refseq_gene)]


### more complex gff loci ----
# Now, move to the XLOCs that have more than one gencode-RefSeq combination
complex_XLOC=unique_gencode_refseq_names_per_XLOC%>%filter(V2%in%dup_XLOC)
# from these, remove the ones that are still unique,
# i.e that have a unique gene name when parsing combinations
complex_XLOC=separate(complex_XLOC,col = 2,sep = ";",into = c("gencode","refseq"))
N_gene_names_complex <- complex_XLOC %>% group_by(V2) %>% summarise(N_gene_names=sum(unique(c(gencode,refseq))!="-"))
table(N_gene_names_complex$N_gene_names==1)
# the big majority still has only one gene_name!
# only 4466 XLOCs merge things with different gene names and need to be revised
still_complex_XLOC <- complex_XLOC%>%filter(V2%in%N_gene_names_complex$V2[N_gene_names_complex$N_gene_names!=1])
# reasons can be:
# A) different naming in GENCODE and RefSeq
# B) different overlapping genes (same biotype)
# C) different overlapping genes (different biotype)

# check refseq names in an XLOC that do not match any
# gencode name in the same XLOC
still_complex_XLOC_refseq_only=still_complex_XLOC%>%filter(gencode=="-")

# Check if the refseq gene matches a gencode gene name in the same loci,
# if so, it means that the refseq is an alternative isoform
# of the same gene, I just use this refseq gene, no conflict
matches_gencode_in_same_loci=mapply(function(x,y)y%in%still_complex_XLOC$gencode[still_complex_XLOC$V2==x],
                                    still_complex_XLOC_refseq_only$V2,still_complex_XLOC_refseq_only$refseq)
still_complex_XLOC_refseq_only$matches_gencode_name_in_same_XLOC=matches_gencode_in_same_loci
table(still_complex_XLOC_refseq_only$matches_gencode_name_in_same_XLOC)

# So I only have to continue with the refseq only ones that DO NOT match
# any gencode gene in the same XLOC:
still_complex_XLOC_refseq_only_nomatch=still_complex_XLOC_refseq_only%>%filter(!matches_gencode_name_in_same_XLOC)

# add coordinates to merged_ref_annot ----
merged_ref_annot=left_join(merged_ref_annot,refGTF%>%filter(type=="transcript")%>%dplyr::select(seqid,start,end,strand,transcript_id),by=c(V1="transcript_id"))

# add new synonyms from matching transcripts ----
# those transcripts that have different gene names in gencode and
# refseq provide new synonyms to unify the gene names
new_synonymdf=synonymsdf%>%select(Name,gene_synonym)
colnames(new_synonymdf)[1]="RefSeq_name"

synonyms2=merged_ref_annot[,c("refseq_gene","gencode_gene_name")]
synonyms2=synonyms2%>%filter(gencode_gene_name!="-"&refseq_gene!="-")
synonyms2=synonyms2%>%filter(gencode_gene_name!=refseq_gene)
synonyms2=synonyms2[!duplicated(synonyms2),]
colnames(synonyms2)=colnames(new_synonymdf)
new_synonymdf=rbind(new_synonymdf,synonyms2)
new_synonymdf=new_synonymdf[!duplicated(new_synonymdf),]

synonymsdf=new_synonymdf
write.table(synonymsdf,"/home/llorenzi/Documentos/references/RefSeq_gene_synonyms.txt",sep = "\t",quote = F,row.names = F)

# check if refseq gene names match with synonyms
matches_synonym_in_same_loci=mapply(function(x,y){
  synons=synonymsdf$gene_synonym[synonymsdf$RefSeq_name%in%y]
  return(any(synons%in%still_complex_XLOC$gencode[still_complex_XLOC$V2==x]))
},
still_complex_XLOC_refseq_only$V2,still_complex_XLOC_refseq_only$refseq)

table(matches_synonym_in_same_loci|matches_gencode_in_same_loci)
table(matches_gencode_in_same_loci)
View(still_complex_XLOC_refseq_only[matches_synonym_in_same_loci&!matches_gencode_in_same_loci,])

# select those complex loci that match only synonym
still_complex_XLOC_refseq_only_matches_synonym <-
  still_complex_XLOC_refseq_only[matches_synonym_in_same_loci&!matches_gencode_in_same_loci,]

# using the XLOC and the RefSeq gene that has a synonym,
# retrieve the corresponding gencode gene name
new_name=mapply(function(x,y){
  synons=synonymsdf$gene_synonym[synonymsdf$RefSeq_name%in%y]
  return(unique(synons[synons%in%still_complex_XLOC$gencode[still_complex_XLOC$V2==x]]))

}, still_complex_XLOC_refseq_only_matches_synonym$V2,
still_complex_XLOC_refseq_only_matches_synonym$refseq)

table(sapply(new_name,length))
which(sapply(new_name,length)>1)
#XLOC_018368
#81
new_name[[81]]
# [1] "Nupl1"   "Gm49336"

# for those cases (this case, actually, is only one), take the first synonym
new_name=sapply(new_name,
               function(x)x[1])

still_complex_XLOC_refseq_only_matches_synonym$new_name=new_name

still_complex_XLOC_refseq_only_nomatch <-
  still_complex_XLOC_refseq_only %>%
  filter(!matches_gencode_name_in_same_XLOC&!matches_synonym_in_same_loci)

View(merged_ref_annot %>% filter(refseq_gene %in% still_complex_XLOC_refseq_only_nomatch$refseq))
#add biotypes
# refseq
still_complex_XLOC_refseq_only_nomatch$refseq_biotype=
  merged_ref_annot$refseq_biotype[match(still_complex_XLOC_refseq_only_nomatch$refseq,
                                        merged_ref_annot$refseq_gene)]

# gencode
still_complex_XLOC_refseq_only_nomatch$gencode_biotype=sapply(still_complex_XLOC_refseq_only_nomatch$V2,
                                                              function(x){
                                                                gcbiots=merged_ref_annot$gencode_biotype[merged_ref_annot$V2==x]
                                                                gcbiots=gcbiots[!is.na(gcbiots)]
                                                                return(paste0(sort(unique(gcbiots)),collapse = ","))})


# correct names for matched synonyms ----
# so at the end there are only 304 cases in which I should change the gene name
# in refseq to match the desired synonym in gencode, these are the ones
# in table still_complex_XLOC_refseq_only_matches_synonym
# For the no match cases, I decided I will leave the ids from refseq as separate
# transcriptional units. It may happen that there are some redundancy
# but I prefer to keep the conflicts/discrepancies between both ref annotations

XLOC_refseq_combined_ids=paste(still_complex_XLOC_refseq_only_matches_synonym$V2,
                               still_complex_XLOC_refseq_only_matches_synonym$refseq,sep = ";")

merged_ref_annot_XLOC_refseq_combined_ids=paste(merged_ref_annot$V2,
                                                merged_ref_annot$refseq_gene,
                                                sep=";")
table(XLOC_refseq_combined_ids%in%merged_ref_annot_XLOC_refseq_combined_ids)

merged_ref_annot$old_gene_name=merged_ref_annot$gene_name

gn2c <- merged_ref_annot$gene_name[merged_ref_annot_XLOC_refseq_combined_ids%in%XLOC_refseq_combined_ids&merged_ref_annot$gencode_gene_name=="-"]

table(gn2c%in%merged_ref_annot$gencode_gene_name)
View(merged_ref_annot[merged_ref_annot_XLOC_refseq_combined_ids%in%XLOC_refseq_combined_ids,])

merged_ref_annot$gene_name[merged_ref_annot_XLOC_refseq_combined_ids%in%XLOC_refseq_combined_ids&
                             merged_ref_annot$gencode_gene_name=="-"] <- still_complex_XLOC_refseq_only_matches_synonym$new_name[match(
                               merged_ref_annot_XLOC_refseq_combined_ids[merged_ref_annot_XLOC_refseq_combined_ids%in%XLOC_refseq_combined_ids&
                                                                           merged_ref_annot$gencode_gene_name=="-"],
                               XLOC_refseq_combined_ids
                             )]
View(merged_ref_annot[merged_ref_annot_XLOC_refseq_combined_ids%in%XLOC_refseq_combined_ids&merged_ref_annot$gencode_gene_name=="-",])

# Check loci consistency ----
# Now I can check the transcripts belonging to each gene_name below to see if there
# are loci that get too long or span differnt locations, and change their names accordingly
# for each gene name check if there are transcripts that do not overlap
# that are further than 1000 bp, if so assign gene_name_1 gene_name_2 to them

# 1) write function to merge transcripts within each gene_name (to check if some of them should be separated)
# 2) modify gene names that are duplicated after doing this
# 3) aggregate by current gene_name and check gene lengths
# 3) plot length distribution and check outliers
# 4) solve conflicts between RefSeq and GENCODE
# 5) simplify biotypes

refGTF_new_gene_id=refGTF
refGTF_new_gene_id$gene_id <-
  merged_ref_annot$gene_name[match(refGTF_new_gene_id$transcript_id,
                                   merged_ref_annot$V1)]

# 1) write function to merge transcripts within each gene_name ----
# (to check if some of them should be separated)
bt_merge_transcripts <- function(bed,dist=1000,mycols="4,5", operations="distinct"){
  require(bedtoolsr)
  options(bedtools.path = "/home/llorenzi/software/bedtools2/bin")

  sorted_bed=bed[order(bed[,1],bed[,2]),]
  merged_tr=bt.merge(i=sorted_bed,s = T,c = mycols,o=operations,d = dist)
  return(merged_tr)
}



bt_merge_transcripts_per_gene <- function(in_bed,gene_id="gene_name"){
  colnames(in_bed)[colnames(in_bed)==gene_id]="gene_name"
  listres=lapply(unique(in_bed$gene_name),function(x){
    bt_merge_transcripts(in_bed[in_bed$gene_name==x,])

  })
  return(do.call("rbind",listres))
}
#bt_merged_transcripts_per_gene=bt_merge_transcripts_per_gene(in_bed)
#write.table(bt_merged_transcripts_per_gene,"bt_merged_transcripts_per_gene_merged_refs.txt",quote = F,sep = "\t",row.names = F,col.names=F)
bt_merged_transcripts_per_gene=read.table("bt_merged_transcripts_per_gene_merged_refs.txt")
sort(table(bt_merged_transcripts_per_gene$V1[grepl("chr",bt_merged_transcripts_per_gene$V1)]))

dupgenes=unique(bt_merged_transcripts_per_gene$V5[duplicated(bt_merged_transcripts_per_gene$V5)])
View(bt_merged_transcripts_per_gene%>%filter(V5%in%dupgenes))
duplicated_bt_merged_transcripts_per_gene <-
  bt_merged_transcripts_per_gene%>%filter(V5%in%dupgenes)

# 2) modify gene names that are duplicated after doing this ----

### deduplication of gene names:

#keep track of the order of original index when
#gene names are alphabeticaly sorted (this is because next
#function will sort names by alphabetic order)
ord <- order(duplicated_bt_merged_transcripts_per_gene$V5)
head(duplicated_bt_merged_transcripts_per_gene[ord,])

#count how many times each duplicated gene is duplicated
ta <- table(duplicated_bt_merged_transcripts_per_gene$V5)
table(ta)

#check that the gene names in original matrix when
#sorted with the ord index are the same than the names in ta
table(unique(duplicated_bt_merged_transcripts_per_gene$V5[ord])==names(ta))

#generate gene suffixes for each gene and the corresponding
#modified gene names
ta_vect <- sapply(ta, function(x)return(paste0("_",seq(1:x))))
modified_names <- paste0(duplicated_bt_merged_transcripts_per_gene$V5[ord],unlist(ta_vect))

#assign the modified gene names to the original matrix in the
#correct order
duplicated_bt_merged_transcripts_per_gene[ord,"modified_gene_names"] <- modified_names

# split the transcripts that correspond to each new gene_name
tr_ids=unlist(strsplit(duplicated_bt_merged_transcripts_per_gene$V4,split = ","))

duplicated_bt_merged_transcripts_per_gene <- duplicated_bt_merged_transcripts_per_gene[rep(1:nrow(duplicated_bt_merged_transcripts_per_gene),
                                                                                        sapply(strsplit(duplicated_bt_merged_transcripts_per_gene$V4,split = ","),length)),]

duplicated_bt_merged_transcripts_per_gene$transcript_id=tr_ids

# change name in annotation for the 118 target transcripts:
merged_ref_annot$gene_name[merged_ref_annot$V1%in%duplicated_bt_merged_transcripts_per_gene$transcript_id] <-
  duplicated_bt_merged_transcripts_per_gene$modified_gene_names[match(merged_ref_annot$V1[merged_ref_annot$V1%in%duplicated_bt_merged_transcripts_per_gene$transcript_id],
                                                                      duplicated_bt_merged_transcripts_per_gene$transcript_id)]





# 3) aggregate by current gene_name and check gene lengths ----
refGTF_new_gene_id$gene_id <- merged_ref_annot$gene_name[match(refGTF_new_gene_id$transcript_id,
                                                               merged_ref_annot$V1)]

refGTF_genomic_coords=trastools::get_genomic_range_by_gene(refGTF_new_gene_id)

gp=ggplot(refGTF_genomic_coords,aes(y=width)) + geom_boxplot() +scale_y_log10()
gp
ggplot(refGTF_genomic_coords, aes(y = width)) +
  geom_boxplot() + scale_y_log10()+
  stat_summary(
    aes(x=0,label = round(stat(y), 1)),
    geom = "text",
    fun.y = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
    hjust = -1
  )
summary(refGTF_genomic_coords$width)

# The longest gene is Fgfr2, but I checked in UCSC and the length is correct...

# 4) Conflicts within GENCODE ----
gencode_biotypes_per_Gene=merged_ref_annot%>%group_by(gene_name) %>%
  summarise(N_biots=length(unique(gencode_biotype[!is.na(gencode_biotype)])),
            biotypes=paste(sort(unique(gencode_biotype)),collapse = ","))

gencode_biotypes_per_Gene=merged_ref_annot%>%group_by(gene_name, gencode_gene) %>%
  summarise(N_biots=length(unique(gencode_biotype[!is.na(gencode_biotype)])),
            biotypes=paste(sort(unique(gencode_biotype)),collapse = ","))

gencode_biotypes_per_Gene <- gencode_biotypes_per_Gene %>% filter(N_biots>0)
gencode_biotypes_per_Gene <- gencode_biotypes_per_Gene[!duplicated(gencode_biotypes_per_Gene[,c("gene_name","biotypes")]),]
multi_biots_gene_names <- gencode_biotypes_per_Gene$gene_name[duplicated(gencode_biotypes_per_Gene$gene_name)]
gencode_biotypes_per_Gene <- gencode_biotypes_per_Gene %>% filter(gene_name%in%multi_biots_gene_names)
length(multi_biots_gene_names)
# There are 62 gencode gene names that have more than one biotype
# I will take the gencode gene id for the secondary biotypes in this cases
unique(gencode_biotypes_per_Gene$biotypes)
ordered_biots <- unique(gencode_biotypes_per_Gene$biotypes)[c(3,2,1,4,7,5,6)]
ordered_biots
# [1] "protein_coding"                     "lncRNA"
# [3] "transcribed_unprocessed_pseudogene" "transcribed_processed_pseudogene"
# [5] "processed_pseudogene"               "unprocessed_pseudogene"
# [7] "TEC"

gencode_biotypes_per_Gene <- gencode_biotypes_per_Gene %>%
  arrange(gene_name,match(biotypes,ordered_biots))

# remove the first biotype from each gene (this is the one that remains unchanged)
gencode_biotypes_per_Gene <- gencode_biotypes_per_Gene %>%
  group_by(gene_name) %>% dplyr::slice(-1)

gencode_biotypes_per_Gene$gene_name_gencode_id <- paste(gencode_biotypes_per_Gene$gene_name,
                                                        gencode_biotypes_per_Gene$gencode_gene)
# now with this table I can change the gene names by gene id in these 62 cases:
gene_name_gencode_id <- paste(merged_ref_annot$gene_name,merged_ref_annot$gencode_gene)

# assign genecode gene ids to these conflicted gene names:
merged_ref_annot$gene_name[gene_name_gencode_id%in%gencode_biotypes_per_Gene$gene_name_gencode_id] <-
  gencode_biotypes_per_Gene$gencode_gene[match(gene_name_gencode_id[gene_name_gencode_id%in%gencode_biotypes_per_Gene$gene_name_gencode_id],
                                               gencode_biotypes_per_Gene$gene_name_gencode_id)]

# 5) Conflicts between RefSeq and GENCODE ----
# Now assign gencode gene biotype and if not available, refseq biotype
merged_ref_annot$gene_biotype=ifelse(!is.na(merged_ref_annot$gencode_biotype),
                                     merged_ref_annot$gencode_biotype,
                                merged_ref_annot$refseq_biotype)

# Now, solve possible conflicts, similar to what I did for gencode biotypes:
biotypes_per_Gene=merged_ref_annot%>%group_by(gene_name, gencode_gene,refseq_gene,refseq_transcript) %>%
  summarise(N_biots=length(unique(gene_biotype[!is.na(gene_biotype)])),
            biotypes=paste(sort(unique(gene_biotype)),collapse = ","))

# put the gencode genes first
biotypes_per_Gene <- biotypes_per_Gene %>% arrange(gene_name,gencode_gene=="-")
# remove duplicated gene_name-biotype
biotypes_per_Gene <- biotypes_per_Gene[!duplicated(biotypes_per_Gene[,c("gene_name","biotypes")]),]
# Most genes now should be only once, we are interested in duplicated ones:
# get duplicated gene_names:
multi_biots_gene_names <- biotypes_per_Gene$gene_name[duplicated(biotypes_per_Gene$gene_name)]
length(multi_biots_gene_names)
# these are 2045, most of these are "Gm..." genes
biotypes_per_Gene <- biotypes_per_Gene %>% filter(gene_name%in%multi_biots_gene_names)
# I will take the GENCODE biotypes for these cases:
biotypes_per_Genef <- biotypes_per_Gene %>% filter(gencode_gene!="-")
length(unique(biotypes_per_Genef$gene_name))
# correct biotype in merged_ref_annot:
old_biots <- merged_ref_annot$gene_biotype[merged_ref_annot$gene_name%in%biotypes_per_Gene$gene_name]
new_biots <- biotypes_per_Gene$biotypes[match(merged_ref_annot$gene_name[merged_ref_annot$gene_name%in%biotypes_per_Gene$gene_name],
                                              biotypes_per_Gene$gene_name)]
View(cbind(old_biots,new_biots))

merged_ref_annot$gene_biotype[merged_ref_annot$gene_name%in%biotypes_per_Gene$gene_name] <-
  biotypes_per_Gene$biotypes[match(merged_ref_annot$gene_name[merged_ref_annot$gene_name%in%biotypes_per_Gene$gene_name],
                                   biotypes_per_Gene$gene_name)]

# 6) simplify biotypes ----
ref_biotypes=unique(merged_ref_annot$gene_biotype)
#simplify biotypes
simpl_biotypes <- list("protein_coding" = c("protein_coding",
                                            grep("_gene|_V_",
                                                 grep("pseudogene",
                                                      ref_biotypes,
                                                      invert = T,value = T),value = T)),
                       "lncRNA" = "lncRNA",
                       "pseudogene" = grep("pseudogene",
                                           ref_biotypes,
                                           value = T),
                       "TEC" = grep("TEC",
                                    ref_biotypes,
                                    value = T),
                       "miRNA" =grep("miRNA",
                                     ref_biotypes,
                                     value = T),
                       "sno-snRNA" = c("snoRNA"  ,
                                       "snRNA",
                                       "ncRNA;snoRNA"),
                       "rRNA" = grep("rRNA",
                                     ref_biotypes,
                                     value = T),
                       "tRNA" = grep("tRNA",
                                     ref_biotypes,
                                     value = T))
#"other" = )

table(ref_biotypes%in%unlist(simpl_biotypes))
ref_biotypes=merged_ref_annot$gene_biotype

ref_biotypes[!ref_biotypes%in%unlist(simpl_biotypes)]

simpl_biotypes <- c(simpl_biotypes,
                    list(other = unique(ref_biotypes)[!unique(ref_biotypes)%in%
                                                        unlist(simpl_biotypes)]))
simpl_biotypes_table <- data.frame(orig.biot=unlist(simpl_biotypes),
                                   final.biot=rep(names(simpl_biotypes),sapply(simpl_biotypes,length)))

merged_ref_annot$simplified_gene_biotype=simpl_biotypes_table$final.biot[match(merged_ref_annot$gene_biotype,
                                                                               simpl_biotypes_table$orig.biot)]

# double check that all genes have a unique biotype
biotypes_per_Gene=merged_ref_annot%>%group_by(gene_name) %>%
  summarise(N_biots=length(unique(gene_biotype[!is.na(gene_biotype)])),
            biotypes=paste(sort(unique(gene_biotype)),collapse = ","))
table(biotypes_per_Gene$N_biots)

colnames(merged_ref_annot)[c(1:2,4)] <- c("gffC_transcript_id",
                                        "gffC_gene_id","gffC_cc")
merged_ref_annot <- merged_ref_annot[,-3]
write.table(merged_ref_annot,"~/Documentos/references/annotated_tracking_file.updated_gene_names.20250121.txt",
            ,sep = "\t",row.names = F,quote = F)

refGTF$gene_name <- merged_ref_annot$gene_name[match(refGTF$transcript_id,
                                                     merged_ref_annot$gffC_transcript_id)]

export(refGTF,"~/Documentos/references/merged_refs.combined.updated_gene_names.20250121.gtf",format = "gtf")

# make symlinks to this data in corresponding script

# Make R object to store updated ref_annot and refGTF
ref_annot=merged_ref_annot
save(ref_annot,refGTF,file="data/RData/merged_ref_annot_GTF.RData")
