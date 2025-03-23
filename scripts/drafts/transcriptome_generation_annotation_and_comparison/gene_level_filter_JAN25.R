library(trastools)
library(optparse)
library(tidyverse)
library(tools)
library(rtracklayer)
ref_annot_path="~/Rprojects/StemLinc_analyses/data/references/merged_refs_annotation/annotated_tracking_file.updated_gene_names.20250121.txt"

# Functions ----
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

# Define the command-line arguments ----
option_list <- list(
  make_option(c("-a", "--annot_track_path"), type="character", default=NULL,
              help="Path to the input annotation file, an annotated tracking file (mandatory)", metavar="character"),
  make_option(c("-g","--gtf_path"),type="character",default = NULL,
              help="Path to the input gtf file (mandatory)", metavar="character"),
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
if (is.null(opt$annot_track_path)|is.null(opt$gtf_path)) {
  stop("Error: annot_track_path and gtf_path are mandatory arguments. Use --annot_track_path and --gtf_path to specify the paths.")
}

# Print the parsed arguments
cat("Parsed Arguments:\n")
cat("annot_track_path:", opt$annot_track_path, "\n")
cat("gtf_path:", opt$gtf_path,"\n")
cat("N_samps:", opt$N_samps, "\n")
cat("tpm_cutoff:", opt$tpm_cutoff, "\n")
cat("suffix:", opt$suffix, "\n")

annot_track_path=opt$annot_track_path
gtf_path=opt$gtf_path
N_samps=opt$N_samps
tpm_cutoff=opt$tpm_cutoff
suffix=opt$suffix

# Load data ----
annotated_tracking <- read.table(annot_track_path,header = T)
ref_annot=read.table(ref_annot_path,header = T)
gtf=readGFF(gtf_path)

# Annotate transfrags against Ref ----
## TEMPORAL: First correct ref annotation columns (error in previous step- need to fix) ----
annotated_tracking$ref_gene_name <-
  ref_annot$gene_name[match(annotated_tracking$ref_transcript_id,
                            ref_annot$gffC_transcript_id)]

annotated_tracking$ref_biotype <-
  ref_annot$simplified_gene_biotype[match(annotated_tracking$ref_transcript_id,
                                          ref_annot$gffC_transcript_id)]

## Assign gene_names to overlapping classcodes, XLOC to non-overlapping ----
annotated_tracking <- annotated_tracking %>%
  mutate(gene_name = ifelse(V4%in%trastools::overlapping_class_codes,
                            ref_gene_name,V2))
## Assign strand of reference transcript ----
annotated_tracking$ref_strand <-
  ref_annot$strand[match(annotated_tracking$ref_transcript_id,
                         ref_annot$gffC_transcript_id)]

#annotated_tracking$width <- annotated_tracking$end - annotated_tracking$start +1

## Assign biotype to potential novel genes ----
annotated_tracking <- annotated_tracking %>%
  mutate(biotype=ifelse(overlap_Ref,ref_biotype,"potNovel"))
cat("Number of initial loci after annotation:","\n")
cat("total:")
cat(length(unique(annotated_tracking$gene_name)),"\n")
cat("annotated:")
cat(length(unique(annotated_tracking$gene_name[!annotated_tracking$overlap_Ref])),"\n")
cat("potNovel:")
cat(length(unique(annotated_tracking$gene_name[annotated_tracking$overlap_Ref])),"\n")
cat("Number of initial transfrags","\n")
cat("total:","\n")
cat(nrow(annotated_tracking),"\n")
cat("overlap ref?:","\n")
cat(table(annotated_tracking$overlap_Ref),"\n")

# Filter 1: remove spurious classcodes - tr level ----
spurious_ccs <- c("e","p","o","s")
filtered_annot <- annotated_tracking %>% filter(!V4%in%spurious_ccs)
cat("Number of transfrags pass 1st filter:","\n")
cat(nrow(filtered_annot),"\n")
cat("Number of loci pass 1st filter:","\n")
cat(length(unique(filtered_annot$gene_name)),"\n")

# Filter 2: remove intronic in sense to PCGs ----
intronic_sense <- filtered_annot %>%
  filter(V4=="i"&ref_biotype=="protein_coding"&ref_strand==strand)
filtered_annot <- filtered_annot %>% filter(!V1%in%intronic_sense$V1)
cat("Number of transfrags pass 2nd filter:","\n")
cat(nrow(filtered_annot),"\n")
cat("Number of loci pass 2nd filter:","\n")
cat(length(unique(filtered_annot$gene_name)),"\n")
# Calculate number of replicates at loci level ----
N_samps_per_gene=get_n_samples_per_gene_from_tracking_gene_name(filtered_annot)

# Filter 3: remove loci supported by less than N reps - loci level----
filtered_annot <- filtered_annot %>% filter(gene_name%in%names(N_samps_per_gene)[N_samps_per_gene==N_samps])
cat("Number of transfrags pass 3rd filter:","\n")
cat(nrow(filtered_annot),"\n")
cat("Number of loci pass 3rd filter:","\n")
cat(length(unique(filtered_annot$gene_name)),"\n")


# Extract max mean_TPM  pre gene ----
max_TPM <- get_max_per_group_by_arrange(filtered_annot,
                                        group_var = gene_name,
                                        arr_var = mean_tpm,
                                        outname = "max_mean_tpm")
# Filter 4: remove loci with max mean_TPM < 0.2 - loci level----
filtered_annot <- filtered_annot %>%
  filter(gene_name%in%max_TPM$gene_name[max_TPM$max_mean_tpm>=tpm_cutoff])

print(filtered_annot %>% group_by(gene_name) %>%
  summarise(biotype=unique(biotype),
            exonic_type=ifelse(any(Nexons>1),"multiexonic","monoexonic")) %>%
  group_by(biotype,exonic_type) %>% summarise(n()))

# Add colnames to filtered annot ----
# Annotate and filter GTF ----
gtf <- gtf%>%filter(transcript_id%in%filtered_annot$V1)
gtf$gene_name <- filtered_annot$gene_name[match(gtf$transcript_id,filtered_annot$V1)]
gtf$old_gene_id <- gtf$gene_id
gtf$gene_id <- gtf$gene_name

# Write outputs ----
## Annotated tracking with gene_names and biotype----
annot_track_out <- paste(file_path_sans_ext(annot_track_path),"gene_names",suffix,file_ext(annot_track_path),sep = ".")
write.table(annotated_tracking,annot_track_out,quote = F,row.names = F,sep = "\t")
## Filtered annotated tracking ----
filtered_annot_out <- paste(file_path_sans_ext(annot_track_path),"filtered",paste0("good_ccs_",N_samps,"reps_",tpm_cutoff,"_tpm"),suffix,file_ext(annot_track_path),sep = ".")
write.table(filtered_annot,filtered_annot_out,quote = F,row.names = F,sep = "\t")
## Filtered annotated GTF ----
gtf_out <- paste(paste0("outputs/gtf/",basename(file_path_sans_ext(gtf_path))),"filtered",paste0("good_ccs_",N_samps,"reps_",tpm_cutoff,"_tpm"),suffix,file_ext(gtf_path),sep = ".")
export(gtf,gtf_out,format = "gtf")
