# Load necessary libraries and sources
files.sources = list.files("scripts/functions/", full.names = T)
sapply(files.sources, source)
library(tidyverse)
library(trastools)

# Define function for writing BED files with filtering options
write_bed_files <- function(data, prefix, class_col, exonic_col, output_dir, score_col, name_col, class_codes = NULL, exonic_types = NULL) {
  unique_classes <- unique(data[[class_col]])
  if (!is.null(class_codes)) unique_classes <- class_codes

  for (cl in unique_classes) {
    bed <- extract_bed(data %>% filter(!!sym(class_col) == cl), name = name_col, score = score_col)
    write_bed(bed, paste0(output_dir, "/", prefix, ".", cl, "_transfrags.bed"))

    if (!is.null(exonic_col)) {
      for (et in unique(data[[exonic_col]])) {
        bed <- extract_bed(data %>% filter(!!sym(class_col) == cl & !!sym(exonic_col) == et), name = name_col, score = score_col)
        write_bed(bed, paste0(output_dir, "/", prefix, ".", cl, "_transfrags.", et, ".bed"))
      }
    }
  }
}

# Define function for promoters
write_promoter_bed_files <- function(data, prefix, class_col, exonic_col, output_dir, score_col, name_col, before = 500, after = 250) {
  unique_classes <- unique(data[[class_col]])

  for (cl in unique_classes) {
    bed <- extract_bed(data %>% filter(!!sym(class_col) == cl), name = name_col, score = score_col)
    bed <- get_promoters_from_bed(bed, before = before, after = after)
    write_bed(bed, paste0(output_dir, "/", prefix, ".", cl, "_promoters.bed"))

    if (!is.null(exonic_col)) {
      for (et in unique(data[[exonic_col]])) {
        bed <- extract_bed(data %>% filter(!!sym(class_col) == cl & !!sym(exonic_col) == et), name = name_col, score = score_col)
        bed <- get_promoters_from_bed(bed, before = before, after = after)
        write_bed(bed, paste0(output_dir, "/", prefix, ".", cl, "_promoters.", et, ".bed"))
      }
    }
  }
}

# Generalized function for processing class codes and exonic types
process_data_to_bed <- function(data, prefix, class_col, exonic_col, score_col, name_col, output_dir, class_codes = NULL) {
  write_bed_files(data, prefix, class_col, exonic_col, output_dir, score_col, name_col, class_codes = class_codes)
}

# Paths
tr_path <- "outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.filtered.20240919_140928.tsv"
gl_path <- "outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240919_140928.tsv"

# Read data
tr <- read.table(tr_path, header = T)
gl <- read.table(gl_path, header = T)

# Process gene level data
gl <- gl %>% filter(pass_filter)
genomic_ranges <- get_genomic_range_by_gene(gtf_df = tr, by = "gene_name")
gl <- left_join(gl, genomic_ranges)
exonic_type <- tr %>% group_by(gene_name) %>% summarise(exonic_type = ifelse(any(Nexons > 1), "multiexonic", "monoexonic"))
gl <- left_join(gl, exonic_type)
max_Nexons <- tr %>% group_by(gene_name) %>% summarise(max_Nexons = max(Nexons))
gl <- left_join(gl, max_Nexons)

# Set output prefixes
tr_out_pref <- gsub(".tsv", "", basename(tr_path))
gl_out_pref <- gsub(".tsv", "", basename(gl_path))

# Process transcript level data
tr <- tr %>% filter(transfrag_class %in% c("protein_coding", "pseudogene", "TEC", "lncRNA", "i", "u", "x", "r"))
tr$exonic_type <- ifelse(tr$Nexons == 1, "monoexonic", "multiexonic")

# Output directories
output_dir_tr <- "outputs/bed_files/tr_level"
output_dir_gl <- "outputs/bed_files/gl_level"
dir.create(output_dir_tr, showWarnings = FALSE)
dir.create(output_dir_gl, showWarnings = FALSE)

# Write transcript level BED files
process_data_to_bed(tr, tr_out_pref, "transfrag_class", "exonic_type", "Nexons", "V1", output_dir_tr)

# Write promoter regions for transcripts
write_promoter_bed_files(tr, tr_out_pref, "transfrag_class", "exonic_type", paste0(output_dir_tr, "/promoter_regions_500_250"), "Nexons", "V1")

# Gene level processing
gl <- gl %>% filter(gene_class %in% c("protein_coding", "pseudogene", "TEC", "lncRNA", "i", "u", "x", "r"))

# Write gene level BED files
process_data_to_bed(gl, gl_out_pref, "gene_class", "exonic_type", "max_Nexons", "gene_name", output_dir_gl)

# Write promoter regions for genes
write_promoter_bed_files(gl, gl_out_pref, "gene_class", "exonic_type", paste0(output_dir_gl, "/promoter_regions_500_250"), "max_Nexons", "gene_name")
