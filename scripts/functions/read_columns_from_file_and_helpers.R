# Define a global variable for column variants
column_id_variants <- list(
  transcript_id = c("transcriptID", "transcript_name", "transcript"),
  gene_id = c("geneID", "gene"),
  gene_name = c("gene_name", "geneSymbol", "gene"),
  biotype = c("biotype", "gene_biotype", "type"),
  seqid = c("chr", "SEQNAMES")
)

# Function to set or modify column_id_variants
# Define or modify column_id_variants in the global environment
set_column_id_variants <- function(new_variants) {
  # Ensure global access to column_id_variants
  if (!exists("column_id_variants", envir = .GlobalEnv)) {
    column_id_variants <<- list()
  }

  # Update or append to the column_id_variants list
  for (name in names(new_variants)) {
    if (name %in% names(column_id_variants)) {
      # Append new variants if they already exist in the global list
      column_id_variants[[name]] <<- unique(c(column_id_variants[[name]], new_variants[[name]]))
    } else {
      # Create a new entry if it doesn't exist
      column_id_variants[[name]] <<- new_variants[[name]]
    }
  }
}


library(dplyr)
library(readr)

check_file_header <- function(file_path, delimiter="\t"){
 names(read_delim(file_path, delim = delimiter, n_max = 1, col_names = TRUE))
}
# Function to read the specified columns based on variants
read_columns_from_file <- function(file_path, columns_to_read = c("transcript_id", "gene_id", "gene_name", "biotype"), delimiter = "\t") {

  # Create an empty list to hold final column names
  final_columns <- vector("list", length(columns_to_read))
  names(final_columns) <- columns_to_read

  # Read the file header to match column names
  file_header <- names(read_delim(file_path, delim = delimiter, n_max = 1, col_names = TRUE))

  # Find the actual column names in the file based on column_id_variants
  for (col in columns_to_read) {
    # Try to find the first matching variant in the file
    variants <- column_id_variants[[col]]
    if (is.null(variants)) {
      variants <- col  # If no variants are set, default to the column itself
    }
    matched_column <- intersect(variants, file_header)

    if (length(matched_column) > 0) {
      final_columns[[col]] <- matched_column[1]  # Use the first match found
    } else {
      message(paste("Warning: No matching column found for", col))
    }
  }

  # Select columns and return as a data frame
  selected_columns <- unlist(final_columns, use.names = FALSE)
  selected_columns <- selected_columns[selected_columns != ""]

  if (length(selected_columns) == 0) {
    stop("No valid columns were found based on the provided variants.")
  }

  # Read only the selected columns from the file
  data <- read_delim(file_path, delim = delimiter, col_names = TRUE, col_types = cols_only(!!!setNames(rep("c", length(selected_columns)), selected_columns)))

  return(data)
}

# # try
# new_variants=list("transcript_id"=c("V1","Tr"),
#                   "gene_id"=c("V2","geneID"))
# set_column_id_variants(new_variants)
# column_id_variants
# file_path="outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.filtered.20240919_140928.tsv"
# tracking_basic <- read_columns_from_file(file_path,columns_to_read = c("transcript_id","gene_id","V4"))
# tracking_basic <- read_columns_from_file(file_path,columns_to_read = c("transcript_id","gene_id","V4"),
#                                          set_given_names = T)
# tracking_basic <- read_columns_from_file(file_path,columns_to_read = c("transcript_ID","gene_id","V4"),
#                                          set_given_names = T)
#
# tracking_basic <- read_columns_from_file(file_path,columns_to_read = c("TR","gene_id","V4"))
