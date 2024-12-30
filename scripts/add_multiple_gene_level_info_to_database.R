############################################
##
## Purpose of script: Add info from multiple sources to tables in database
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-12-16
##
## Email: lucialorenzi90@gmail.com
##
#  Notes ---------------------------
##
##
##
#  Setup ---------------------------

options(scipen = 999)
library(rtracklayer)
library(tidyverse)
library(purrr)
library(DBI)
library(RSQLite)
library(dbplyr)

db <- dbConnect(RSQLite::SQLite(), "/home/llorenzi/Rprojects/StemLinc_analyses/outputs/dbs/StemLinc.db")

# Data collapsed in "outputs/gene_level_info_all_evidence_oct24.csv" ----
multiple_info <- read.csv("outputs/gene_level_info_all_evidence_oct24.csv")

colnames(multiple_info)

# Identify all tables
all_tables <- dbListTables(db)

# Identify tables that are gene-level (e.g., by name or metadata table if available)
gene_level_tables <- grep("gene",all_tables, value = T)

# Get the fields of each gene-level table
gene_level_fields <- lapply(gene_level_tables, function(table) {
  dbGetQuery(db, paste0("PRAGMA table_info(", table, ");"))$name
})

# Combine all fields into a single vector
all_gene_fields <- unique(unlist(gene_level_fields))
print(all_gene_fields)

# Expression data ----
cols=c(1,21:24,41:46)
colnames(multiple_info)[cols]
gene_expression_data <- multiple_info[,cols]
head(gene_expression_data)

dbExecute(db, "
    CREATE TABLE gene_expression_data (
        gene_name TEXT PRIMARY KEY,
        Nhigher_than_hist INTEGER,
        mean_vst_LSK_StL REAL,
        mean_vst_HSPC_polyA REAL,
        FC_vst_LSK_StL_vs_HSPC_polyA REAL,
        HSPC_vs_diff_logFC REAL,
        HSPC_vs_diff TEXT,
        HSPC_vs_prog_logFC REAL,
        HSPC_vs_prog TEXT,
        prog_vs_diff_logFC REAL,
        prog_vs_diff TEXT,
        FOREIGN KEY (gene_name) REFERENCES gene_core_data(gene_name)
    );
")

colnames(gene_expression_data) <- gsub("\\.","_",colnames(gene_expression_data))

dbWriteTable(
  conn = db,
  name = "gene_expression_data",
  value = gene_expression_data,
  overwrite=T,  # Append data to the existing table
  row.names = FALSE
)

metadata <- data.frame(
  table_name = "gene_expression_data",
  script_name = "add_multiple_gene_level_info_to_database.R",  # Replace with your script name
  level = "gene",
  filter_status = "3 replicates gene level, >= 0.2 TPM, only biotypes of interest",
  analysis_version = "v1.0",
  timestamp = Sys.time(),
  source_data = "outputs/gene_level_info_all_evidence_oct24.csv",
  orig_scripts = "multiple: including DGEA.R, collect_all_evidence.R",  # Replace with your actual script(s)
  description = "Gene expression data including mean expression values, fold changes, and differential expression results for LSK, HSPC, and differentiation stages. Linked to gene_core_data by gene_name."
)

dbWriteTable(db, "metadata", metadata, append = TRUE, row.names = FALSE)

# Correlation with signature genes ----
cols=c(1,18,25:28)
gene_correlation_with_PCG_signature <- multiple_info[,cols]
head(gene_correlation_with_PCG_signature)

dbExecute(db, "
    CREATE TABLE gene_correlation_with_PCG_signature (
        gene_name TEXT PRIMARY KEY,
        mean_corr_signature REAL,
        max_Signature_gene TEXT,
        max_corr_signature REAL,
        maxcorrPCG TEXT,
        max_corr REAL,
        FOREIGN KEY (gene_name) REFERENCES gene_core_data(gene_name)
    );
")

dbWriteTable(
  conn = db,
  name = "gene_correlation_with_PCG_signature",
  value = gene_correlation_with_PCG_signature,
  append = TRUE,  # Append data into the existing table
  row.names = FALSE
)

# Preview the table
dbGetQuery(db, "SELECT * FROM gene_correlation_with_PCG_signature LIMIT 10;")

metadata <- data.frame(
  table_name = "gene_correlation_with_PCG_signature",
  script_name = "add_multiple_gene_level_info_to_database.R",
  level = "gene",
  filter_status = "3 replicates gene level, >= 0.2 TPM, only biotypes of interest",
  analysis_version = "v1.0",
  timestamp = Sys.time(),
  source_data = "outputs/gene_level_info_all_evidence_oct24.csv",
  orig_scripts = "correlation_with_signature.R, collect_all_evidence.R" ,
  description = "Gene correlation with signature genes and all PCGs, including mean and maximum correlations and associated PCGs."
)

dbWriteTable(db, "metadata", metadata, append = TRUE, row.names = FALSE)

# Overlap chromatin marks ----
colnames(multiple_info)
cols=c(1,29:40)
colnames(multiple_info)[cols]
gene_overlap_with_chromatin_marks <- multiple_info[,cols]
head(gene_overlap_with_chromatin_marks)

dbExecute(db, "
    CREATE TABLE gene_overlap_with_chromatin_marks (
        gene_name TEXT PRIMARY KEY,
        evidence_level TEXT,
        Enhancer BOOLEAN,
        CAGE_within_100bp BOOLEAN,
        polyAsite_within_100bp BOOLEAN,
        H3K4me3_at_promoter BOOLEAN,
        H3K4me1_at_promoter BOOLEAN,
        H3K27ac_at_promoter BOOLEAN,
        H3K36me3_at_geneBody BOOLEAN,
        Enhancer_Atlas_at_promoter BOOLEAN,
        FANTOM_enhancers_at_promoter BOOLEAN,
        Enhancer_from_marks_at_promoter BOOLEAN,
        H3K27me3_at_promoter BOOLEAN,
        FOREIGN KEY (gene_name) REFERENCES gene_core_data(gene_name)
    );
")

# Add the table to the database
dbWriteTable(
  conn = db,
  name = "gene_overlap_with_chromatin_marks",
  value = gene_overlap_with_chromatin_marks,
  append = TRUE,
  row.names = FALSE
)

# Preview the first rows
dbGetQuery(db, "SELECT * FROM gene_overlap_with_chromatin_marks LIMIT 10;")

metadata <- data.frame(
  table_name = "gene_overlap_with_chromatin_marks",
  script_name = "add_multiple_gene_level_info_to_database.R",
  level = "gene",
  filter_status = "3 replicates gene level, >= 0.2 TPM, only biotypes of interest",
  analysis_version = "v1.0",
  timestamp = Sys.time(),
  source_data = "outputs/gene_level_info_all_evidence_oct24.csv",
  orig_scripts = "multiple overlap scripts and collect_all_evidence.R",  # Replace with your script
  description = "Gene overlap with chromatin marks, CAGE peaks, polyA site and enhancers."
)

# Append the metadata to the database
dbWriteTable(db, "metadata", metadata, append = TRUE, row.names = FALSE)

# Overlap with transcription factors ----
cols=c(1,47,56)
colnames(multiple_info)[cols]
gene_overlap_TFs_cistromes <- multiple_info[,cols]
dbExecute(db, "
    CREATE TABLE gene_overlap_TFs_cistromes_rellevA (
        gene_name TEXT PRIMARY KEY,
        TF_rellevA TEXT,
        N_TFs INTEGER,
        FOREIGN KEY (gene_name) REFERENCES gene_core_data(gene_name)
    );
")

colnames(gene_overlap_TFs_cistromes)[2] <- "TF_rellevA"
# Add the table to the database
dbWriteTable(
  conn = db,
  name = "gene_overlap_TFs_cistromes_rellevA",
  value = gene_overlap_TFs_cistromes,
  append = TRUE,
  row.names = FALSE
)
# Preview the first rows
dbGetQuery(db, "SELECT * FROM gene_overlap_TFs_cistromes_rellevA LIMIT 10;")

metadata <- data.frame(
  table_name = "gene_overlap_TFs_cistromes_rellevA",
  script_name = "add_multiple_gene_level_info_to_database.R",
  level = "gene",
  filter_status = "3 replicates gene level, >= 0.2 TPM, only biotypes of interest",
  analysis_version = "v1.0",
  timestamp = Sys.time(),
  source_data = "outputs/gene_level_info_all_evidence_oct24.csv",
  orig_scripts = "overlap_TFs_cistromes.R and collect_all_evidence.R",
  description = "Gene overlap with transcription factors from Cistromes doi: 10.1186/s13104-018-3856-x"
)

# Append the metadata to the database
dbWriteTable(db, "metadata", metadata, append = TRUE, row.names = FALSE)


# Overlap with previous studies Delas_Luo ----
cols=c(1,48:54)
colnames(multiple_info)[cols]
gene_overlap_Delas_Luo <- multiple_info[,cols]
gene_overlap_Delas_Luo <- gene_overlap_Delas_Luo %>% filter(rowSums(is.na(gene_overlap_Delas_Luo))<
                                                             ncol(gene_overlap_Delas_Luo) - 1)

colnames(gene_overlap_Delas_Luo)[c(2,4,5,7,8)] <- c("Delas_best_gene","Delas_alt_name","Delas_orientation",
                                                  "Luo_best_gene","Luo_enrichment")

head(gene_overlap_Delas_Luo)

dbExecute(db, "
    CREATE TABLE gene_overlap_Delas_Luo (
        gene_name TEXT PRIMARY KEY,
        Delas_best_gene TEXT,
        Delas_best_cc TEXT,
        Delas_alt_name TEXT,
        Delas_orientation TEXT,
        Delas_enrichment TEXT,
        Luo_best_gene TEXT,
        Luo_enrichment TEXT,
        FOREIGN KEY (gene_name) REFERENCES gene_core_data(gene_name)
    );
")

dbWriteTable(
  conn = db,
  name = "gene_overlap_Delas_Luo",
  value = gene_overlap_Delas_Luo,
  append = TRUE,
  row.names = FALSE
)

# Preview the table
dbGetQuery(db, "SELECT * FROM gene_overlap_Delas_Luo LIMIT 10;")

metadata <- data.frame(
  table_name = "gene_overlap_Delas_Luo",
  script_name = "add_multiple_gene_level_info_to_database.R",
  level = "gene",
  filter_status = "3 replicates gene level, >= 0.2 TPM, only biotypes of interest, only genes with overlap with Delas and/or Luo",
  analysis_version = "v1.0",
  timestamp = Sys.time(),
  source_data = "outputs/gene_level_info_all_evidence_oct24.csv",
  orig_scripts = "overlap_Delas_Luo.R, collect_all_evidence.R",  # Replace with your actual script
  description = "Gene overlap data with Delas and Luo studies, including enrichment in cell groups."
)

# Append metadata to the database
dbWriteTable(db, "metadata", metadata, append = TRUE, row.names = FALSE)


# Results from ABC ----
#check script "analyse_ABC_results.R"

# MaxEntScan splicing sites analyses ----
#check "analyse_splicing_motifs.R"

# Overlap repeats ----
