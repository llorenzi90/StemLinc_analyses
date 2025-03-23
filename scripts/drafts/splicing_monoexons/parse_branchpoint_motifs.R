# Load libraries ----
library(tidyverse)
library(glue)
library(ggseqlogo)

# Load branchpoint data ----
branchpoint_table <- readxl::read_xlsx("data/public_data/splicing_branchpoint/Allison_et_al_2017/Supplemental_Table_S3.xlsx")
table(branchpoint_table$`Motif Model`)

length(unique(branchpoint_table$`Branchpoint Sequence - ( ) around bulge,     * - BP nucleotide`))
length(unique(branchpoint_table$`Branchpoint Sequence - ( ) around bulge,     * - BP nucleotide`[branchpoint_table$`Motif Model`=="canonical"]))

branchpoint_table %>% group_by(`Motif Model`) %>%
  summarise(total=n(),
    unique=length(unique(`Branchpoint Sequence - ( ) around bulge,     * - BP nucleotide`)))

branchpoint_table <-  branchpoint_table %>% mutate(Sequence = gsub("\\(|\\)|\\*","",`Branchpoint Sequence - ( ) around bulge,     * - BP nucleotide`))

# Generate sequence logos for each Motif model ----
sequences <- branchpoint_table$Sequence
names(sequences) <- branchpoint_table$`Motif Model`

table(nchar(branchpoint_table$Sequence))
length_distribution <- branchpoint_table %>% mutate(seq_len=nchar(Sequence)) %>%
                                                      group_by(`Motif Model`, seq_len) %>%
                                                    summarize(count = n(), .groups = "drop")

for (n in unique(names(sequences))) {
  g=ggseqlogo(sequences[names(sequences)==n])
  print(g + ggtitle(n))
}

pdf(file = "outputs/motifs/splicing_motifs/Allison_branchpoint_logos.pdf",onefile = T)
for (n in unique(names(sequences))) {
  g=ggseqlogo(sequences[names(sequences)==n])
  print(g + ggtitle(n))
}
dev.off()

pdf(file = "outputs/motifs/splicing_motifs/Allison_branchpoint_logos.unique_seqs.pdf",onefile = T)
for (n in unique(names(sequences))) {
  tmp_seqs=sequences[names(sequences)==n]
  g=ggseqlogo(tmp_seqs[!duplicated(tmp_seqs)])
  print(g + ggtitle(n))
}
dev.off()
# Generate meme formatted motif file ----
## Step 1: Input data for each motif model ----
all_motifs <- branchpoint_table %>%
  rename(
    MotifModel = `Motif Model`,
    BranchpointSequence = Sequence
  ) %>%
  group_by(MotifModel) %>%
  summarize(
    sequences = list(BranchpointSequence),
    .groups = "drop"
  )

unique_motifs <- branchpoint_table %>%
  rename(
    MotifModel = `Motif Model`,
    BranchpointSequence = Sequence
  ) %>%
  group_by(MotifModel) %>%
  summarize(
    unique_sequences = list(unique(BranchpointSequence)),
    .groups = "drop"
  )


## Step 2: Create MEME file content ----
create_meme_motif <- function(motif_name, sequences) {
  num_seqs <- length(sequences)
  seq_length <- nchar(sequences[1])  # Assumes all sequences have the same length

  # Build MEME motif header
  header <- glue("
MOTIF {motif_name}
letter-probability matrix: alength= 4 w= {seq_length} nsites= {num_seqs} E= 0.0
")

  # Create letter-probability matrix
  matrix <- sapply(1:seq_length, function(pos) {
    base_counts <- table(substring(sequences, pos, pos))
    all_bases <- c("A", "C", "G", "T")
    freqs <- sapply(all_bases, function(base) {
      count <- base_counts[base]
      ifelse(is.na(count), 0, count)  # Replace NA with 0
    })
    freqs / num_seqs
  })

  matrix_content <- apply(matrix, 2, function(col) paste(sprintf("%.3f", col), collapse = " "))

  # Combine header and matrix
  glue("{header}\n{paste(matrix_content, collapse = '\\n')}")
}

## Step 3: Generate motifs for all categories (unique and duplicated) ----
meme_motifs_unique <- unique_motifs %>%
  mutate(motif_content = purrr::map2(
    MotifModel, unique_sequences, create_meme_motif
  ))

meme_motifs_all <- all_motifs %>%
  mutate(motif_content = purrr::map2(
    MotifModel, sequences, create_meme_motif
  ))

## Step 4: Write to MEME file ----

out_list <- list(output_unique=list(meme_motifs_unique,
                                    "outputs/motifs/splicing_motifs/branchpoint_motifs.unique_sequences.meme"),
                 output_all=list(meme_motifs_all,
                                    "outputs/motifs/splicing_motifs/branchpoint_motifs.all_sequences.meme")
)

sapply(names(out_list),function(o){
  meme_motifs <- out_list[[o]][[1]]
  output_file <- out_list[[o]][[2]]
  write_lines(
    c(
      "MEME version 5.0",
      "ALPHABET= ACGT",
      "strands: + -",
      "Background letter frequencies
    A 0.25 C 0.25 G 0.25 T 0.25",
      unlist(meme_motifs$motif_content)
    ),
    output_file
  )

  cat(glue("Motif file saved to: {output_file}\n"))

})
