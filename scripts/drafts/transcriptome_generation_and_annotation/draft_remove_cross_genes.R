library(GenomicRanges)
library(rtracklayer)
library(dplyr)

# 1. Load GTF File
annot_track_path <- "outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.tsv"
gtf_file <- '/home/llorenzi/Rprojects/StemLinc_analyses/data/raw/LSK_StemLinc.combined.gtf'  # Replace with your actual GTF file path
gtf <- import(gtf_file)

ref_gtf_file <- '/home/llorenzi/Rprojects/StemLinc_analyses/data/references/merged_refs_annotation/merged_refs.combined.updated_gene_names.20250121.gene_names.gtf'
ref_gtf <- import(ref_gtf_file)

ref_transcripts <- ref_gtf[ref_gtf$type == "transcript"]
# 3. Generate Gene Coordinates
# Group transcripts by gene_id and calculate min(start) and max(end) per gene
ref_gene_ranges <- ref_transcripts %>%
  as.data.frame() %>%
  group_by(gene_id, seqnames, strand) %>%
  summarise(start = min(start), end = max(end), .groups = "drop")

# 4. Convert DataFrame to GRanges Object
ref_genes <- GRanges(
  seqnames = ref_gene_ranges$seqnames,
  ranges = IRanges(start = ref_gene_ranges$start, end = ref_gene_ranges$end),
  strand = ref_gene_ranges$strand
)
ref_genes$gene_names=ref_gene_ranges$gene_id
# 2. Extract gene and transcript features
transcripts <- gtf[gtf$type == "transcript"]

# 3. Find Transcripts Overlapping Multiple Genes in Sense Direction
#   - Overlapping means: transcript has exons overlapping >1 gene
#   - Sense direction means transcript and genes are on the same strand

overlap_counts <- countOverlaps(transcripts, ref_genes, ignore.strand = FALSE)  # Keep strand information
head(transcripts$transcript_id,20)
head(overlap_counts,20)
cbind(head(transcripts$transcript_id,20),
      head(overlap_counts,20))
View(cbind(transcripts$transcript_id,overlap_counts))
# 4. Filter Out "Cross-Gene" Transcripts
clean_transcripts <- transcripts[overlap_counts == 1]  # Keep only transcripts overlapping 1 gene

annotated_tracking <- read.table(annot_track_path,header = T)

# 5. Keep Only Related Exons and CDS Features
# Find transcript IDs to keep
valid_transcript_ids <- unique(mcols(clean_transcripts)$transcript_id)

# Filter GTF to keep only valid transcripts, their exons, and CDS
gtf_cleaned <- gtf[mcols(gtf)$transcript_id %in% valid_transcript_ids]

# 6. Save the Cleaned GTF
output_file <- "cleaned_output.gtf"
export(gtf_cleaned, con = output_file, format = "gtf")

cat("Filtered GTF saved as:", output_file, "\n")
