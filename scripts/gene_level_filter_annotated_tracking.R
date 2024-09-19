library(trastools)
library(optparse)
library(tidyverse)

# Define the command-line arguments
option_list <- list(
  make_option(c("-f", "--file_path"), type="character", default=NULL,
              help="Path to the input file (mandatory)", metavar="character"),
  make_option(c("-n", "--N_samps"), type="integer", default=3,
              help="Number of samples in which each gene is found (default: %default)", metavar="integer"),
  make_option(c("-t", "--tpm_cutoff"), type="double", default=0.2,
              help="TPM cutoff value to filter genes (>=tpm_cutoff) (default: %default)", metavar="double"),
  make_option(c("-s", "--suffix"), type="character", default=format(Sys.time(), "%Y%m%d_%H%M%S"),
              help="Suffix for output (default: timestamp %default)", metavar="character")
)

# Create an option parser and parse the arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if the mandatory file_path argument is provided
if (is.null(opt$file_path)) {
  stop("Error: file_path is a mandatory argument. Use --file_path to specify the path.")
}

# Print the parsed arguments
cat("Parsed Arguments:\n")
cat("file_path:", opt$file_path, "\n")
cat("N_samps:", opt$N_samps, "\n")
cat("tpm_cutoff:", opt$tpm_cutoff, "\n")
cat("suffix:", opt$suffix, "\n")

ordered_classcodes <- c("=", "j", "k", "c", "m", "n", "e", "o", "p", "s", "x", "i", "y", "u", "r", ".")
ref_annot_path="~/Rprojects/StemLinc_analyses/data/references/merged_refs_annotation/annotated_tracking_file.updated_gene_names.txt"
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
  genes2keep <- gene_level_info %>% filter(Nsamps==Nsamps_per_gene&max_mean_tpm>=maxTPM) %>% pull(gene_name)
  gene_level_info$pass_filter <- gene_level_info$gene_name%in%genes2keep
  filtered_tracking <- annotated_tracking %>% filter(gene_name %in% genes2keep)
  return(list(filtered_tracking=filtered_tracking,
              gene_level_info=gene_level_info))
}


file_path=opt$file_path
N_samps=opt$N_samps
tpm_cutoff=opt$tpm_cutoff
suffix=opt$suffix
# Read data
annotated_tracking <- read.table(file_path,header = T)

# Filter the data based on the number of assembled samples and TPM cutoff
filt_list <- filter_transcriptome(annotated_tracking,Nsamps_per_gene = N_samps,maxTPM = tpm_cutoff)

# Write data
filtered_output_path <- paste0(gsub("tsv","",file_path),"filtered.",suffix,".tsv")
gene_level_info_output_path <- paste0(gsub("tsv","",file_path),"gene_level_info.",suffix,".tsv")


cat("Writing filtered tracking output to:", filtered_output_path, "\n")
write.table(filt_list$filtered_tracking, filtered_output_path, quote = F,row.names = F, sep = "\t")

cat("Writing filtered tracking output to:", gene_level_info_output_path, "\n")
write.table(filt_list$gene_level_info, gene_level_info_output_path, quote = F,row.names = F, sep = "\t")
