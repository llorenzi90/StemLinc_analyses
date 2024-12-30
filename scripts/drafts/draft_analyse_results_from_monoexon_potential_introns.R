############################################
##
## Purpose of script: explore potential splicing motifs in monoexons (miss monique in the background)
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-12-27
##
## Email: lucialorenzi90@gmail.com
##
#  Notes ---------------------------
##  Exploratory analysis of the test
##  results from "~/Documentos/draft_find_splicing_sites_in_monoexons.py"
##
##  The main idea is to see if there are differences between occurrence of splicing motifs
##  between monoexon sequences (mainly from lncRNAs) and random sequences of similar size
##  (with no overlap with  any transcribed region in the assembled genome),
##  if we see differences, e.g random regions have stronger potential introns,
##  what does it mean?
##  We may also see no differences, what does it mean?
##
##  Improvements if promising: 1. add other motifs, not only splicing, even de novo discovery, if found
##  any enrichment relative to random regions, then try to map this motif to known motifs
##  (TFs, RBPs, miRNAs, what esle?)
##  2. correct for base composition: select random regions with
##  similar GC content and also dinucleotide compoisition
##
##
#  Setup ---------------------------

options(scipen = 999)
# Load packages---------------------------
# load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
library(rtracklayer)
library(bedtoolsr)
library(ggplot2)
library(DBI)
# Load data---------------------------

db <- dbConnect(RSQLite::SQLite(), "/home/llorenzi/Rprojects/StemLinc_analyses/outputs/dbs/StemLinc.db")
setwd("~/Documentos/test_outputs_100/")

input_bed <- read.table("/home/llorenzi/Rprojects/StemLinc_analyses/outputs/bed_files/tr_level/LSK_StemLinc.combined.sorted.chr1_monoexons.100_random.bed")
donors <- read.delim("test_outputs/monoexon_analysis_test_donor_sites.tsv")
acceptors <- read.delim("test_outputs/monoexon_analysis_test_acceptor_sites.tsv")
potential_introns <- read.delim("test_outputs/monoexon_analysis_test_potential_introns.tsv")
# script python: revisar extended end y extended start in donors, add monoexon length,
# add bed file of random sequences, add GC content of each sequence, maybe add masking of NNN
# fasta of random regions, check 0-based...

table(grepl("random",donors$seq_ID)[!duplicated(donors$seq_ID)])
table(grepl("random",acceptors$seq_ID)[!duplicated(acceptors$seq_ID)])
table(grepl("random",potential_introns$seq_ID)[!duplicated(potential_introns$seq_ID)])

# 97 vs 39 (monoexons vs random regions)
values_per_feature <- donors %>% mutate(type = ifelse(grepl("random",seq_ID),"random","monoexon")) %>% group_by(seq_ID,type) %>% summarise(donors_per_monoexon = n())

values_per_feature %>% group_by(type) %>% summarise(mean_donors_per_feat = mean(donors_per_monoexon))

values_per_feature <- full_join(values_per_feature,acceptors %>% mutate(type = ifelse(grepl("random",seq_ID),"random","monoexon")) %>% group_by(seq_ID,type) %>% summarise(acceptors_per_monoexon = n()))

values_per_feature %>% group_by(type) %>% summarise(mean_donors_per_feat =
                                                      mean(donors_per_monoexon),
                                                    mean_acceptors_per_feat =
                                                      mean(acceptors_per_monoexon))
values_per_feature <- values_per_feature %>%
  mutate(donors_per_monoexon = replace_na(donors_per_monoexon,0))

values_per_feature <- full_join(values_per_feature,
                                potential_introns %>%
                                  mutate(type =
                                           ifelse(grepl("random",seq_ID),
                                                  "random",
                                                  "monoexon")) %>%
                                  group_by(seq_ID,type) %>%
                                  summarise(introns_per_monoexon = n()))

values_per_feature$introns_per_monoexon[is.na(values_per_feature$introns_per_monoexon)] <- 0
# values_per_feature <- values_per_feature %>%
#   mutate(introns_per_monoexon = replace_na(introns_per_monoexon, 0))

values_per_feature %>% group_by(type) %>% summarise(mean_donors_per_feat =
                                                      mean(donors_per_monoexon),
                                                    mean_acceptors_per_feat =
                                                      mean(acceptors_per_monoexon),
                                                    mean_pot_introns_per_feat =
                                                      mean(introns_per_monoexon))

values_per_feature <- full_join(values_per_feature,
                                donors %>%
                                  mutate(type =
                                           ifelse(grepl("random",seq_ID),
                                                           "random",
                                                           "monoexon"),
                                         width = extended_end - extended_start) %>%
                                  group_by(seq_ID,type) %>%
                                  summarise(width = width[1]))


missing_widths <-     acceptors %>%
                                  mutate(width = extended_end - extended_start) %>%
                                  group_by(seq_ID) %>%
                                  summarise(width = width[1])

values_per_feature$width[is.na(values_per_feature$width)] <-
  missing_widths$width[match(values_per_feature$seq_ID[is.na(values_per_feature$width)],
                             missing_widths$seq_ID)]
values_per_feature %>% group_by(type) %>% summarise(mean_width =
                                                     mean(width),
                                                   mean_donors_per_feat_width =
                                                     mean(donors_per_monoexon/width),
                                                   mean_acceptors_per_feat_width =
                                                     mean(acceptors_per_monoexon/width),
                                                   mean_pot_introns_per_feat_width =
                                                     mean(introns_per_monoexon/width))
# load maxEntScan results ----
ss3 <- read.table("/home/llorenzi/Documentos/test_outputs_100/test_outputs/monoexon_analysis_test_3ss.MaxEntScan.txt")
ss5 <- read.table("/home/llorenzi/Documentos/test_outputs_100/test_outputs/monoexon_analysis_test_5ss.MaxEntScan.txt")
#

table(ss3$V1%in%acceptors$acceptor_site_id)
table(ss5$V1%in%donors$donor_site_id)

potential_introns$donor_strength <- ss5$V3[match(potential_introns$donor_site_id,
                                                 ss5$V1)]
potential_introns$acceptor_strength <- ss3$V3[match(potential_introns$acceptor_site_id,
                                                    ss3$V1)]
potential_introns$combined_strength <- potential_introns$donor_strength + potential_introns$acceptor_strength
summary(potential_introns$combined_strength)

potential_introns <- potential_introns %>% mutate(type = ifelse(grepl("random",seq_ID),"random","monoexon"))

ggplot(potential_introns, aes(x=type, y=combined_strength)) + geom_boxplot()
library(ggpubr)

# Create the boxplot with statistical comparison
ggplot(potential_introns, aes(x = type, y = combined_strength)) +
  geom_boxplot() +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(c("monoexon", "random")))

# Subset data by 'monoexon' and 'random' types
monoexon_data <- potential_introns$combined_strength[potential_introns$type == "monoexon"]
random_data <- potential_introns$combined_strength[potential_introns$type == "random"]

test_result <- wilcox.test(monoexon_data, random_data)

# Compare max intron strength per sequence ----
max_potential_introns <- potential_introns %>% group_by(seq_ID,type) %>%
  summarise(max_combined_strength = max(combined_strength))

ggplot(max_potential_introns,aes(x=type, y = max_combined_strength)) + geom_boxplot()+
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(c("monoexon", "random")))

monoexon_data <- max_potential_introns$max_combined_strength[max_potential_introns$type == "monoexon"]
random_data <- max_potential_introns$max_combined_strength[max_potential_introns$type == "random"]

test_result <- wilcox.test(monoexon_data, random_data)

# Add annotation ----
# Split the seq_ID column
strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', potential_introns$seq_ID[1:10]), ' ')

potential_introns <- potential_introns %>%
  mutate(
    split_seq_ID = sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', seq_ID)
  ) %>%
  separate(
    col = split_seq_ID,
    into = c("transcript_id", "gene_name"),
    sep = " "
  )

dbListTables(db)
dbListFields(db, "transcript_core_data")
transcript_core_data <- dbReadTable(db,"transcript_core_data")
potential_introns$biotype <- transcript_core_data$gene_biotype[match(potential_introns$transcript_id,
                                                                     transcript_core_data$transcript_id)]

table(potential_introns$biotype[!duplicated(potential_introns$transcript_id)])
# filter potential intronic sites ----
# keep only potential introns with combined strength >5
table(potential_introns$combined_strength>5)

maxEntScan_all_genes <- dbReadTable(db, "gene_splice_sites_maxMaxEntScan")
metadata <- dbReadTable(db,"metadata")
maxEntScan_all_genes$combined_strength <- maxEntScan_all_genes$max_ss5_MaxEntScan +
  maxEntScan_all_genes$max_ss3_MaxEntScan
# Things to check:
# strength by monoexon biotype
# strength compared to real introns in each biotype (a sample of 100)
# check how many of the extended sequences match the consensus sequence
# analyse only splicing signals inside monoexons (i.e without extension)
# check my notes for the script and try to add extension as a command line input
# select only monoexons and potNovel with additional evidence level

# select randomly 100 introns from maxEntScan_all_genes
set.seed(1234)
rand100_real_introns <- maxEntScan_all_genes[sample(nrow(maxEntScan_all_genes),100),]

rand100_real_introns$biotype <- transcript_core_data$gene_biotype[match(rand100_real_introns$gene_name,
                                                                        transcript_core_data$gene_name)]

table(rand100_real_introns$biotype)

max_introns <- rand100_real_introns %>% select(seq_ID=gene_name,
                                               type=biotype,
                                               max_combined_strength=combined_strength)

merged_data <- rbind(max_potential_introns, max_introns)
ggplot(merged_data,aes(x=type, y = max_combined_strength)) + geom_boxplot()+
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(c("monoexon", "random"),
                                                                                    c("potNovel","monoexon"),
                                                                                    c("lncRNA","random"),
                                                                                    c("potNovel","random")))

# add biotype to monoexons
max_potential_introns$biotype <- potential_introns$biotype[match(max_potential_introns$seq_ID,
                                                                 potential_introns$seq_ID)]

max_potential_introns <- max_potential_introns %>% mutate(type = ifelse(type != "random",
                                                                        paste0("mono_",biotype),
                                                                        "random"))

merged_data <- rbind( max_potential_introns %>% select(seq_ID,type,max_combined_strength),
                      max_introns)

ggplot(merged_data,aes(x=type, y = max_combined_strength)) + geom_boxplot()+
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(c("mono_potNovel", "random"),
                                                                                   c("mono_lncRNA","random"),
                                                                                   c("potNovel","random")))

# Important things to check:
# 1. GC content ask cgpt
# 2. Normalize by length
# 3. Do it with all sequences (but only biotypes of interest)
# 4. Check if random sequences have potential issues...(do not think so)
# 5. Generate bed file with random sequences, useful for other analyses, e.g calculate "promoters"
# and TF binding
# 6. Repeat analyses with only stringent set of monoexons
