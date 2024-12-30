############################################
##
## Purpose of script: Set up a SQL database to organize gene data
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-12-11
##
## Email: lucialorenzi90@gmail.com
##
#  Notes ---------------------------
##
##
##
#  Setup ---------------------------

options(scipen = 999)
# Load packages---------------------------
# load up the packages we will need:  (uncomment as required)
require(tidyverse)
library(trastools)
# require(data.table)
# library(rtracklayer)
# library(bedtoolsr)
# library(ggplot2)

# Install db packages ----
install.packages("DBI")
install.packages("RSQLite")
library(DBI)
library(RSQLite)
library(dbplyr)
library(rtracklayer)
# Create database ----
# Create dir to store db
dir.create("outputs/dbs")
# Specify the full path for the database file
db <- dbConnect(RSQLite::SQLite(), "/home/llorenzi/Rprojects/StemLinc_analyses/outputs/dbs/StemLinc.db")


# Load data---------------------------

# original GTF:
gtf_path <- "data/raw/LSK_StemLinc.combined.sorted.gtf"
GTF <- readGFF(gtf_path)
# gene_level data, table created from annotated track (gffcompare), script ~/Rprojects/StemLinc_analyses/scripts/gene_level_filter_annotated_tracking.R
raw_gene_level_info <- read.table('/home/llorenzi/Rprojects/StemLinc_analyses/outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.tsv',header = T)
table(raw_gene_level_info$biotype)
table(is.na(raw_gene_level_info$biotype))

# transcript level data, sources: StringTie and gffcompare
raw_transcript_level_info <- read.table("outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.tsv", header = T)

# annotation of the reference transcriptome used, a merge between GENCODE and RefSeq
ref_annot_path="~/Rprojects/StemLinc_analyses/data/references/merged_refs_annotation/annotated_tracking_file.updated_gene_names.txt"
ref_annot=read.table(ref_annot_path,header = T)

# Parse transcipt level info ----
colnames(raw_transcript_level_info)[1:2] <- c("transcript_id","StringTie_gene_id")
colnames(raw_transcript_level_info)[5] <- c("classcode")

# Add meaningful reference ids from merged reference annotation table
raw_transcript_level_info$gencode_transcript <- ref_annot$gencode_transcript[match(raw_transcript_level_info$ref_transcript_id,
                                                                                   ref_annot$V1)]

raw_transcript_level_info$refseq_transcript <- ref_annot$refseq_transcript[match(raw_transcript_level_info$ref_transcript_id,
                                                                                   ref_annot$V1)]

raw_transcript_level_info$gencode_gene <- ref_annot$gencode_gene[match(raw_transcript_level_info$ref_transcript_id,
                                                                                   ref_annot$V1)]

raw_transcript_level_info$refseq_gene <- ref_annot$refseq_gene[match(raw_transcript_level_info$ref_transcript_id,
                                                                                 ref_annot$V1)]

raw_transcript_level_info$gencode_gene_name <- ref_annot$gencode_gene_name[match(raw_transcript_level_info$ref_transcript_id,
                                                                       ref_annot$V1)]
# Add gene name that matches the gene level annotation
raw_transcript_level_info <- raw_transcript_level_info %>% mutate(gene_name = ifelse(overlap_Ref,
                                                                       ref_gene_name,StringTie_gene_id))

# Add biotype "potNovel" for potential novel genes
raw_transcript_level_info <- raw_transcript_level_info %>% mutate(gene_biotype = ifelse(overlap_Ref,
                                                                                     ref_biotype,"potNovel"))


# Change GTF gene_name by gene_name in transcript annotation
GTF$gene_name <- raw_transcript_level_info$gene_name[match(GTF$transcript_id,
                                                           raw_transcript_level_info$transcript_id)]
# Make sample related info col names more descriptive
colnames(raw_transcript_level_info)[6:8]=c("rep1_info","rep2_info","rep3_info")
colnames(raw_transcript_level_info)[14:16] <- c("rep1_tpm","rep2_tpm","rep3_tpm")


# Parse the gene level info ----
# core data: "gene_name"     "chr"           "start"         "end"           "strand"        "N_total_exons"       "gene_biotype"

# I have to calculate the total exons per gene from the gtf
epg <- exons_per_gene(GTF%>% select(-gene_id),gene_col = "gene_name")

raw_gene_level_info <- left_join(raw_gene_level_info,
                                 epg)


# Create and populate gene level tables ----
# Create the gene_core_data table in the database
dbExecute(db, "
    CREATE TABLE gene_core_data (
        gene_name TEXT PRIMARY KEY,
        chr TEXT,
        start INTEGER,
        end INTEGER,
        strand TEXT,
        width INTEGER,
        N_total_exons INTEGER,
        biotype TEXT
    );
")

raw_gene_level_info$biotype[is.na(raw_gene_level_info$biotype)] <- "potNovel"

# Generate the gene_core_data table
gene_core_data <- raw_gene_level_info %>%
  dplyr::select(
    gene_name,
    chr = seqnames,
    start,
    end,
    strand,
    width,
    N_total_exons = N_exons,
    biotype,
    exonic_type
  )


# Write the data to the database
dbWriteTable(db, "gene_core_data", gene_core_data, append = TRUE, row.names = FALSE)


# add also the exonic type
exonic_type <- GTF %>% filter(type=="exon") %>% group_by(transcript_id, gene_name) %>%
  summarise(N_exons=n())

exonic_type <- exonic_type %>% group_by(gene_name) %>%
  summarise(exonic_type=ifelse(any(N_exons>1),"multiexonic","monoexonic"))

dbExecute(db, "ALTER TABLE gene_core_data ADD COLUMN exonic_type TEXT;")

#
gene_core_data <- dbReadTable(db, "gene_core_data")
gene_core_data$exonic_type <- exonic_type$exonic_type[match(gene_core_data$gene_name,
                                                            exonic_type$gene_name)]

dbWriteTable(db, "gene_core_data", gene_core_data, overwrite = TRUE, row.names = FALSE)


# Basic gene annotation
# from the table add Nsamps as N_samples, best_cc as best_classcode, max_mean_tpm as max_mean_tpm_LSK,
# pass_filter
# Add reference gene ids and names

gene_level_ref_annot <- raw_transcript_level_info %>% group_by(gene_name) %>%
  summarise(StringTie_gene_id=paste(unique(StringTie_gene_id),collapse = ","),
            gencode_gene=paste(unique(gencode_gene[gencode_gene!="-"]),collapse = ","),
            refseq_gene=paste(unique(refseq_gene[refseq_gene!="-"]),collapse = ","),
            gencode_gene_name=paste(unique(gencode_gene_name[gencode_gene_name!="-"]),collapse = ",")
            )

gene_level_ref_annot <- left_join(gene_level_ref_annot, raw_gene_level_info %>%
                                    dplyr::select(gene_name,
                                                  best_classcode = best_cc,
                                                  ))

#

# Create gene_level_ref_annotation table
dbExecute(db, "
    CREATE TABLE gene_level_ref_annotation (
        gene_name TEXT PRIMARY KEY,
        StringTie_gene_id TEXT,
        gencode_gene TEXT,
        refseq_gene TEXT,
        gencode_gene_name TEXT,
        best_classcode TEXT,
        FOREIGN KEY (gene_name) REFERENCES gene_core_data (gene_name)

    );
")


dbWriteTable(db, "gene_level_ref_annotation", gene_level_ref_annotation, append = TRUE, row.names = FALSE)

# Create gene_level_filter table
dbExecute(db, "
    CREATE TABLE gene_level_filter (
        gene_name TEXT PRIMARY KEY,
        N_replicates INTEGER,
        max_mean_tpm_LSK REAL,
        pass_filter BOOLEAN,
        FOREIGN KEY (gene_name) REFERENCES gene_core_data (gene_name)
    );
")


gene_level_filter <- raw_gene_level_info %>%
  dplyr::select(
    gene_name,
    N_replicates = Nsamps,
    max_mean_tpm_LSK = max_mean_tpm,
    pass_filter
  )

dbWriteTable(db, "gene_level_filter", gene_level_filter, append = TRUE, row.names = FALSE)

# Create and populate transcript tables ----

# create core transcript data
dbExecute(db, "
    CREATE TABLE transcript_core_data (
        transcript_id TEXT PRIMARY KEY,
        gene_name TEXT,
        chr TEXT,
        start INTEGER,
        end INTEGER,
        strand TEXT,
        N_exons INTEGER,
        gene_biotype TEXT,
        FOREIGN KEY (gene_name) REFERENCES gene_core_data (gene_name)
    );
")




# Prepare data for the transcript_core_data table
transcript_core_data <- raw_transcript_level_info %>%
  select(
    transcript_id,
    gene_name ,
    chr = seqid,
    start,
    end,
    strand,
    N_exons = Nexons,
    gene_biotype
  )


# Populate the transcript_core_data table
dbWriteTable(db, "transcript_core_data", transcript_core_data, append = TRUE, row.names = FALSE)


# Create the transcript_annotation table
dbExecute(db, "
    CREATE TABLE transcript_annotation (
        transcript_id TEXT PRIMARY KEY,
        StringTie_gene_id TEXT,
        classcode TEXT,
        rep1_info TEXT,
        rep2_info TEXT,
        rep3_info TEXT,
        overlap_Ref TEXT,
        Nsamps INTEGER,
        rep1_tpm REAL,
        rep2_tpm REAL,
        rep3_tpm REAL,
        mean_tpm REAL,
        gencode_transcript TEXT,
        refseq_transcript TEXT,
        gencode_gene TEXT,
        refseq_gene TEXT,
        gencode_gene_name TEXT,
        FOREIGN KEY (transcript_id) REFERENCES transcript_data (transcript_id)
    );
")


# Prepare data for transcript_annotation table
transcript_annotation <- raw_transcript_level_info %>%
  select(
    transcript_id,
    StringTie_gene_id,
    classcode,
    rep1_info,
    rep2_info,
    rep3_info,
    overlap_Ref,
    Nsamps,
    rep1_tpm,
    rep2_tpm,
    rep3_tpm,
    mean_tpm,
    gencode_transcript,
    refseq_transcript,
    gencode_gene,
    refseq_gene,
    gencode_gene_name
  )

# Populate the transcript_annotation table
dbWriteTable(db, "transcript_annotation", transcript_annotation, append = TRUE, row.names = FALSE)


# List all tables in the database
tables <- dbListTables(db)

# Loop through each table and list its fields
for (table in tables) {
  cat("Table:", table, "\n")
  fields <- dbGetQuery(db, paste0("PRAGMA table_info(", table, ");"))
  print(fields[, c("name", "type")]) # Display field names and types
  cat("\n")
}


# Confirm that all desired info are in created tables:
table(colnames(raw_transcript_level_info)[!colnames(raw_transcript_level_info)%in%c(colnames(transcript_core_data),
                                                                                    colnames(transcript_annotation))])

# Create metadata table ----
# Create the metadata table
dbExecute(db, "
    CREATE TABLE metadata (
        table_name TEXT PRIMARY KEY,
        script_name TEXT,
        level TEXT,
        filter_status TEXT,
        analysis_version TEXT,
        timestamp TEXT,
        source_data TEXT,
        orig_scripts TEXT,
        description TEXT
    );
")


# Populate the metadata table
# Populate the metadata table
# Populate the metadata table with updated filter status
metadata <- data.frame(
  table_name = c("gene_core_data", "gene_level_filter", "gene_level_ref_annotation",
                 "transcript_core_data", "transcript_annotation"),
  script_name = "set_SQL.R",
  level = c("gene", "gene", "gene", "transcript", "transcript"),
  filter_status = c("none", "contains pass_filter column", "none", "none", "none"),
  analysis_version = c("v1.0", "v1.0", "v1.0", "v1.0", "v1.0"),
  timestamp = Sys.time(),
  source_data = c(
    "annotated_tracking_file_gene_level:outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.tsv, GTF:data/raw/LSK_StemLinc.combined.sorted.gtf",
    "annotated_tracking_file_gene_level:outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.tsv",
    "annotated_tracking_transcript_level:outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.tsv, annotated_tracking_file_gene_level:outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.tsv",
    "annotated_tracking_transcript_level:outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.tsv, annotation_of_merged_refs:data/references/merged_refs_annotation/annotated_tracking_file.updated_gene_names.txt",
    "annotated_tracking_transcript_level:outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.tsv, annotation_of_merged_refs:data/references/merged_refs_annotation/annotated_tracking_file.updated_gene_names.txt"
  ),
  orig_scripts = c(
    "StringTie and gffCompare (hpc) and gene_level_filter_annotated_tracking.R",
    "gene_level_filter_annotated_tracking.R",
    "StringTie and gffCompare (hpc), Transcriptome_characterization_report.Rmd and gene_level_filter_annotated_tracking.R",
    "StringTie and gffCompare (hpc) and Transcriptome_characterization_report.Rmd",
    "StringTie and gffCompare (hpc) and Transcriptome_characterization_report.Rmd"
  ),
  description = c(
    "Core gene-level data with basic attributes like coordinates, strand, and biotype.",
    "Gene-level data including replicate counts, max TPM, and pass filter status.",
    "Gene-level reference annotations including Gencode and RefSeq IDs and gffCompare classcodees.",
    "Core transcript-level data including transcript ID, coordinates, and gene biotype.",
    "Transcript-level annotations with TPM values and overlaps with reference."
  )
)

# Write metadata to the database
dbWriteTable(db, "metadata", metadata, append = TRUE, row.names = FALSE)


# Add gene classification info ----
gene_classif <- read.table("/home/llorenzi/Rprojects/StemLinc_analyses/outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.tsv", header = T)

head(gene_classif)
closest_PCG=read.table("outputs/closest_PCGs.txt",header = T)

closest_PCG <- closest_PCG %>% group_by(gene_name= idQ) %>%
  arrange(abs(dist)) %>%
  summarise(distance_closest_PCG = dist[1])


# Select the new fields
gene_classification_to_PCG <- gene_classif %>%
  select(
    gene_name, best_classif_to_PCG,
    best_closest_PCG, classif_to_PCG, closest_PCG
  )

gene_classification_to_PCG <- left_join(gene_classification_to_PCG, closest_PCG)

dbExecute(db, "
    CREATE TABLE gene_classification_to_PCG (
        gene_name TEXT PRIMARY KEY,
        best_classif_to_PCG TEXT,
        best_closest_PCG TEXT,
        classif_to_PCG TEXT,
        closest_PCG TEXT,
        distance_closest_PCG INTEGER,
        FOREIGN KEY (gene_name) REFERENCES gene_core_data (gene_name)
    );
")

# Write the new data to the database
dbWriteTable(db, "gene_classification_to_PCG", gene_classification_to_PCG, append = TRUE, row.names = FALSE)


# Create metadata for the gene_level_classification table
new_metadata <- data.frame(
  table_name = "gene_classification_to_PCG",
  script_name = "set_SQL.R",
  level = "gene",
  filter_status = "3 replicates gene level, >= 0.2 TPM",
  analysis_version = "v1.0",
  timestamp = Sys.time(),
  source_data = "gene_classif table from classification to PCG analysis",
  orig_scripts = "gene_classification_to_PCG.R",
  description = "Gene-level classification relative to closest PCG."
)

# Write the metadata to the metadata table
dbWriteTable(db, "metadata", new_metadata, append = TRUE, row.names = FALSE)





