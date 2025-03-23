library(trastools)
library(optparse)
library(tidyverse)
library(rtracklayer)

#ref_annot_path="~/Rprojects/StemLinc_analyses/data/references/merged_refs_annotation/annotated_tracking_file.updated_gene_names.txt"
ref_annot_path="data/references/merged_refs_annotation/annotated_tracking_file.updated_gene_names.20250121.txt"


ordered_classcodes <- c("=", "j", "k", "c", "m", "n", "e", "o", "p", "s", "x", "i", "y", "u", "r", ".")
ref_annot=read.table(ref_annot_path,header = T)

get_n_samples_per_gene_from_tracking_gene_name <- function(tracking,gene_id="gene_name",cols=6:8){
  tracking <- as.data.frame(tracking)
  tracking$V2=tracking[,gene_id]
  get_n_samples_per_gene_from_tracking(tracking,cols = cols)
}

get_max_per_group_by_arrange <- function(tracking,
                                         group_var,
                                         arr_var,
                                         summ_var={{arr_var}},
                                         descending = T,
                                         outname="max_val"){

  tracking <- tracking %>% group_by({{group_var}})

  if(descending){
    tracking <- tracking %>% arrange(desc({{arr_var}}))
  } else tracking <- tracking %>% arrange({{arr_var}})

  tracking <- tracking %>% summarise(!!outname:= {{summ_var}}[1])

  return(tracking)
}


filter_transcriptome <- function(annotated_tracking,Nsamps_per_gene=3,maxTPM=0.2){
  annotated_tracking <- annotated_tracking %>% mutate(gene_name = ifelse(V4%in%trastools::overlapping_class_codes,
                                                                         ref_gene_name,V2))
  N_samps_per_gene=get_n_samples_per_gene_from_tracking_gene_name(annotated_tracking)
  N_samps_per_gene <- data.frame(gene_name=names(N_samps_per_gene),Nsamps=N_samps_per_gene)

  max_TPM <- get_max_per_group_by_arrange(annotated_tracking,
                                          group_var = gene_name,
                                          arr_var = mean_tpm,
                                          outname = "max_mean_tpm")
  gene_level_info <- left_join(N_samps_per_gene,max_TPM)
  best_cc_per_gene=annotated_tracking %>% group_by(gene_name) %>% arrange(match(V4,ordered_classcodes)) %>%
    summarise(best_cc=V4[1])

  gene_level_info <- left_join(gene_level_info,best_cc_per_gene)
  gene_level_info$biotype=ref_annot$simplified_gene_biotype[match(gene_level_info$gene_name,
                                                                  ref_annot$gene_name)]
  gene_level_info <- gene_level_info %>% mutate(gene_class=ifelse(is.na(biotype),best_cc,biotype))

  # filter
  gene_level_genomic_range <- trastools::get_genomic_range_by_gene(annotated_tracking,by = "gene_name")
  gene_level_info <- left_join(gene_level_info, gene_level_genomic_range)
  genes2keep <- gene_level_info %>% filter(Nsamps==Nsamps_per_gene&max_mean_tpm>=maxTPM) %>% pull(gene_name)
  gene_level_info$pass_filter <- gene_level_info$gene_name%in%genes2keep
  filtered_tracking <- annotated_tracking %>% filter(gene_name %in% genes2keep)

  return(list(filtered_tracking=filtered_tracking,
              gene_level_info=gene_level_info))
}

# code ----
#file_path=opt$file_path
file_path="outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.tsv"
N_samps=3
tpm_cutoff=0.2
suffix="test"
# Read data
annotated_tracking <- read.table(file_path,header = T)

annotated_tracking$ref_gene_name <- ref_annot$gene_name[match(annotated_tracking$ref_transcript_id,
                                                              ref_annot$gffC_transcript_id)]

annotated_tracking$ref_biotype <- ref_annot$simplified_gene_biotype[match(annotated_tracking$ref_transcript_id,
                                                              ref_annot$gffC_transcript_id)]
annotated_tracking <- annotated_tracking %>% mutate(gene_name = ifelse(V4%in%trastools::overlapping_class_codes,
                                                                       ref_gene_name,V2))

annotated_tracking$ref_strand <- ref_annot$strand[match(annotated_tracking$ref_transcript_id,
                                                    ref_annot$gffC_transcript_id)]

annotated_tracking$width <- annotated_tracking$end - annotated_tracking$start +1


annotated_tracking <- annotated_tracking %>% mutate(biotype=ifelse(overlap_Ref,ref_biotype,"potNovel"))

# remove monoexonic isoforms whose closest match is multiexonic ----
ref_gtf <- readGFF("data/references/merged_refs_annotation/merged_refs.combined.updated_gene_names.20250121.gtf")
# table(ref_gtf$gene_name == ref_annot$gene_name[match(ref_gtf$transcript_id,
#                                                      ref_annot$gffC_transcript_id)])

ref_annot_Nexons <- ref_gtf %>% filter(type=="exon") %>% group_by(transcript_id) %>%
  summarise(Nexons=n())

annotated_tracking$Nexons_ref <- ref_annot_Nexons$Nexons[match(annotated_tracking$ref_transcript_id,
                                                               ref_annot_Nexons$transcript_id)]

monoexons_overlap_multiexonRef <- annotated_tracking %>%filter(overlap_Ref &
                                                                 Nexons==1 &
                                                                 Nexons_ref!=1)
# 7307 transcripts
View(monoexons_overlap_multiexonRef)
table(monoexons_overlap_multiexonRef$V4)
# c    e    m    n    o
# 2604 2806  190  575 1132
# the majority are 'e' transfrags labeled as "possible pre-mRNA fragment", then contained,
# "other" and finally the rest are intron retains (m and n)
table(monoexons_overlap_multiexonRef$ref_biotype)
# lncRNA          other protein_coding     pseudogene            TEC           tRNA
# 3198              8           3760            274             64              3

monoexons_overlap_multiexonRef %>% group_by(gene_name) %>%
  summarise(biotype=paste(unique(ref_biotype),collapse = ",")) %>%
  group_by(biotype) %>% summarise(n())

# A tibble: 6 × 2
# biotype        `n()`
# <chr>          <int>
# 1 TEC               15
# 2 lncRNA          1930
# 3 other              6
# 4 protein_coding  2318
# 5 pseudogene       179
# 6 tRNA               2

# there are 1930 lncRNA genes that get monoexonic isoforms where they have multiexonic isofrms

# check how many of the recovered lncRNAs have only monoexonic assembled isoforms
# select lncRNA loci for which at least one of its multiexonic isoforms was recovered by monoexon assembled isoforms
monoexons_overlap_multiexonRef_lncRNA_genes <- monoexons_overlap_multiexonRef$ref_gene_name[monoexons_overlap_multiexonRef$ref_biotype=="lncRNA"]
# select all other transcripts
complement_monoexons_overlap_multiexonRef <- annotated_tracking %>% filter(!V1%in%monoexons_overlap_multiexonRef$V1)
# how many of the first are still found on the second, meaning they are recovered by a multiexonic isoform?
table(unique(monoexons_overlap_multiexonRef_lncRNA_genes)%in%complement_monoexons_overlap_multiexonRef$ref_gene_name)
# most of them are still recovered, only 441(23%) are lost

# therefore is a good idea to filter out these
# 1st filter ----
# transfrag level
# remove monoexon isoforms in annotated multiexons
filtered_annot <- complement_monoexons_overlap_multiexonRef

unique(monoexons_overlap_multiexonRef_lncRNA_genes)[!unique(monoexons_overlap_multiexonRef_lncRNA_genes)%in%complement_monoexons_overlap_multiexonRef$ref_gene_name]

# check classcodes distribution remaining after this first filter ----
table(filtered_annot$V4)
# I will remove spurious classcodes: e, p, o, and s, while s are minimal, all e have been removed in
#  the previous step, there are 5510 "p" and 994 "o". Most 'p' isoforms will be sense_downstream
table(filtered_annot$ref_biotype,filtered_annot$V4)
# Most of them are close to PCGs. I want to avoid these two close sense_downstream
# potential novel lncRNAs because it is difficult to determine if they are really independent transcriptional units
# OK, I could use CAGE, but should perform a separate analyses for those. For now, exclude.
# 2nd filter ----
# transfrag level
spurious_ccs <- c("e","p","o","s")
filtered_annot <- filtered_annot %>% filter(!V4%in%spurious_ccs)

# Now I also want to remove the intronic sense transcripts to protein_coding, for the same reasons above
# 3rd filter ----
# transfrag level
intronic_sense <- filtered_annot %>% filter(V4=="i"&ref_biotype=="protein_coding" &ref_strand==strand)
table(intronic_sense$Nexons)
length(unique(intronic_sense$gene_name))
#These are 17856 transfrags from 15841 loci, the great majority of which are monoeoxnic.
summary(intronic_sense$width)
# with a median length of 357 and a mean of 1105

filtered_annot <- filtered_annot %>% filter(!V1%in%intronic_sense$V1)


# Decide on how to remove monoexons from mixed new loci
# check how many mixed new loci
new_transfrags <- filtered_annot %>% filter(!overlap_Ref)

new_transfrags_exonic_types <- new_transfrags %>% group_by(gene_name) %>%
  summarise(Nmono=sum(Nexons==1),
            Nmulti=sum(Nexons>1),
            exonic_class=ifelse(Nmono>0,ifelse(Nmulti>0,"mixed","monoexonic"),"multiexonic"))

table(new_transfrags_exonic_types$exonic_class)
# View "mixed" loci
View(filtered_annot %>% filter(gene_name%in%new_transfrags_exonic_types$gene_name[new_transfrags_exonic_types$exonic_class=="mixed"]))

N_samps_per_gene=get_n_samples_per_gene_from_tracking_gene_name(filtered_annot)
new_transfrags_exonic_types$Nsamps <-
  N_samps_per_gene[match(new_transfrags_exonic_types$gene_name,names(N_samps_per_gene))]

table(new_transfrags_exonic_types$exonic_class,new_transfrags_exonic_types$Nsamps)
# The great majority of not annotated novel loci are found in only one sample
# This is not true for mixed class
# in which case there are more loci present in 2 samples

# test how many multiexonic are left if I remove monoexonic variants in new multiexonic loci ----
# split by loci again (to be done)
#old_filtered_annot_GL <- read.delim("~/Documentos/references/old_outputs_transcriptome_char_StemLinc/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.biotypes_of_interest.tsv")

mixed_gene_names <- new_transfrags_exonic_types$gene_name[new_transfrags_exonic_types$exonic_class=="mixed"]
temp_filtered_annot <- filtered_annot %>% filter(!(gene_name%in%mixed_gene_names&Nexons==1))
# recalculate number of samples after having removed the monoexons from mixed loci:
N_samps_per_gene=get_n_samples_per_gene_from_tracking_gene_name(temp_filtered_annot)
new_transfrags_exonic_types$Nsamps_removed_monoexons <-
  N_samps_per_gene[match(new_transfrags_exonic_types$gene_name,names(N_samps_per_gene))]

new_transfrags_exonic_types_mixed <- new_transfrags_exonic_types %>%
  filter(exonic_class=="mixed")
table(new_transfrags_exonic_types_mixed$Nsamps,
      new_transfrags_exonic_types_mixed$Nsamps_removed_monoexons)
# 1   2   3
# 1  17   0   0
# 2 392  31   0
# 3 173 133  78
# If I do this, I will remove 173 multiexonic new loci that are supported by 1 sample
# and 133 that are supported by 2 samples

# Although a bit complex, I could easily implement the rule:
# for mixed potential novel loci, require that after removing monoexonic isoforms,
# the gene is still assembled in at least two samples, i.e if in two out
# of three samples I find splicing, I am ok with a third sample only supporting
# with monoexonic isoform in the same locus

# So the steps would be: filter by number of samples and then
# remove monoexons from mixed, further filter those that are left with only one sample

# 4th filter ----
# loci level
# retain only genes assembled in 3 samples
length(unique(filtered_annot$gene_name))
#[1] 95759
N_samps_per_gene=get_n_samples_per_gene_from_tracking_gene_name(filtered_annot)

filtered_annot <- filtered_annot %>% filter(gene_name%in%names(N_samps_per_gene)[N_samps_per_gene==3])
length(unique(filtered_annot$gene_name))
#[1] 27754
# this leaves us with a nicely reduced number of genes
table(filtered_annot$V4)
table(filtered_annot$overlap_Ref)

table(filtered_annot$biotype)
table(filtered_annot$biotype,filtered_annot$V4)

# Now filter by minimum expression, as I did before
# 5th filter ----
# gene level
# Expression filter
max_TPM <- get_max_per_group_by_arrange(filtered_annot,
                                        group_var = gene_name,
                                        arr_var = mean_tpm,
                                        outname = "max_mean_tpm")

filtered_annot <- filtered_annot %>%
  filter(gene_name%in%max_TPM$gene_name[max_TPM$max_mean_tpm>=0.2])

filtered_annot %>% group_by(gene_name) %>% summarise(biotype=unique(biotype)) %>%
  group_by(biotype) %>% summarise(n())


# 6th filter ----
# transfrag level
# remove monoexonic isoforms from multiexonic potNovel
new_transfrags <- filtered_annot %>% filter(!overlap_Ref)

new_transfrags_exonic_types <- new_transfrags %>% group_by(gene_name) %>%
  summarise(Nmono=sum(Nexons==1),
            Nmulti=sum(Nexons>1),
            exonic_class=ifelse(Nmono>0,ifelse(Nmulti>0,"mixed","monoexonic"),"multiexonic"))

table(new_transfrags_exonic_types$exonic_class)
# View "mixed" loci
View(filtered_annot %>% filter(gene_name%in%new_transfrags_exonic_types$gene_name[new_transfrags_exonic_types$exonic_class=="mixed"]))

N_samps_per_gene=get_n_samples_per_gene_from_tracking_gene_name(filtered_annot)
new_transfrags_exonic_types$Nsamps <-
  N_samps_per_gene[match(new_transfrags_exonic_types$gene_name,names(N_samps_per_gene))]

table(new_transfrags_exonic_types$exonic_class,new_transfrags_exonic_types$Nsamps)
mixed_gene_names <- new_transfrags_exonic_types$gene_name[new_transfrags_exonic_types$exonic_class=="mixed"]

filtered_annot <- filtered_annot %>% filter(!(gene_name%in%mixed_gene_names&Nexons==1))

# recalculate number of samples after having removed the monoexons from mixed loci:
N_samps_per_gene=get_n_samples_per_gene_from_tracking_gene_name(filtered_annot)
new_transfrags_exonic_types$Nsamps_removed_monoexons <-
  N_samps_per_gene[match(new_transfrags_exonic_types$gene_name,names(N_samps_per_gene))]

new_transfrags_exonic_types_mixed <- new_transfrags_exonic_types %>%
  filter(exonic_class=="mixed")
table(new_transfrags_exonic_types_mixed$Nsamps,
      new_transfrags_exonic_types_mixed$Nsamps_removed_monoexons)

# 7th filter
# I will leave it like this for now...


# summarise at gene level ----
filtered_annot_gene_level <- filtered_annot %>% group_by(gene_name) %>%
  arrange(match(V4,ordered_classcodes))%>%
  summarise(exonic_type=ifelse(any(Nexons>1),"multiexonic","monoexonic"),
            biotype=unique(biotype),
            best_cc=V4[1],
            max_mean_tpm=max(mean_tpm),
            chr=seqid[1],
            start=min(start),
            end=max(end),
            strand=unique(strand),
            width=end - start +1
                             )

# write data ----
dir.create("outputs/transcriptome_characterization/jan25")
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_name <- gsub(".tsv",paste0(".filtered.",ts,".tsv"), basename(file_path))
write.table(filtered_annot,paste0("outputs/transcriptome_characterization/jan25/",out_name),
            quote = F,row.names = F,sep="\t")
out_name_gl <- gsub(".tsv",".gene_level.tsv",out_name)
write.table(filtered_annot_gene_level,
            paste0("outputs/transcriptome_characterization/jan25/",out_name_gl),
            quote = F,row.names = F,sep="\t")

# select only biotypes of interest
biots=c("protein_coding","lncRNA","potNovel","pseudogene","TEC")
filtered_annot_sel_biots <- filtered_annot %>% filter(biotype%in%biots)
filtered_annot_gene_level_sel_biots <- filtered_annot_gene_level %>% filter(biotype%in%biots)
out_name <- gsub(".tsv",paste0(".filtered.sel_biots.",ts,".tsv"), basename(file_path))
write.table(filtered_annot_sel_biots,
            paste0("outputs/transcriptome_characterization/jan25/",out_name),
            quote = F,row.names = F,sep="\t")
out_name_gl <- gsub(".tsv",".gene_level.tsv",out_name)
write.table(filtered_annot_gene_level_sel_biots,
            paste0("outputs/transcriptome_characterization/jan25/",out_name_gl),
            quote = F,row.names = F,sep="\t")


# write filtered gtf ----
gtf <- readGFF("data/raw/LSK_StemLinc.combined.gtf")
table(filtered_annot$V1%in%gtf$transcript_id)
gtf <- gtf%>% filter(transcript_id%in%filtered_annot$V1)
gtf$gene_name <- filtered_annot$gene_name[match(gtf$transcript_id,
                                                filtered_annot$V1)]

export(gtf,"outputs/gtf/LSK_StemLinc.combined.filtered.20250127_112131.gtf",format = "gtf")
gtf$gene_id=gtf$gene_name
export(gtf,"outputs/gtf/LSK_StemLinc.combined.filtered.20250127_112131.gene_name.gtf",format = "gtf")
