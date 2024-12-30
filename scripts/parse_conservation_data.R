library(rtracklayer)
library(tidyverse)
library(purrr)
library(DBI)
library(RSQLite)
library(dbplyr)
phastcons_files <- list.files("outputs/bed_files/",pattern = "mm39.phastCons35way",full.names = T)

phastcons <- lapply(phastcons_files,read.table )


phastcons <- lapply(phastcons, function(cons_res){
  colnames(cons_res) <- c("name","size","covered_size","sum_cov","mean_noncov0","mean_cov")
  return(cons_res)
}
  )
names(phastcons) <- c("merged_exons","promoters","exons","introns")

# merged exons conservation ----
# merged exons conservation gives a more neat representation
# of the mean exonic sequence conservation of a gene
phastcons$merged_exons <- separate(phastcons$merged_exons,col = 1,into = c("gene_name","merged_exon_number"),
                                   sep = "_(?!.*_)")

summary_merged_exons_conservation <- phastcons$merged_exons %>% group_by(gene_name) %>%
  arrange(-mean_cov) %>% summarise(phastCons_most_conserved_merged_exon = mean_cov[1],
                                   phastCons_mean_merged_exons = mean(mean_cov))

# promoter conservation ----
# match conservation info with corresponding gene and transcript names
# load GTF and transcript annotation
gtf_path <- "data/raw/LSK_StemLinc.combined.sorted.gtf"
GTF <- readGFF(gtf_path)
phastcons$promoters <- left_join(phastcons$promoters, GTF %>% filter(type=="transcript") %>%
                                   dplyr::select(transcript_id,gene_name), by=c(name="transcript_id"))

summary_promoters_conservation <- phastcons$promoters %>% group_by(gene_name) %>%
  arrange(-mean_cov) %>% summarise(phastCons_most_conserved_promoter = mean_cov[1],
                                   most_conserved_promoter_transcript_id = paste(name[mean_cov==max(mean_cov)],
                                                                                 collapse = ","),
                                   phastCons_mean_promoters = mean(mean_cov))


# exon conservation ----
# match conservation info with corresponding gene and transcript names
# load GTF and transcript annotation

exons_gtf <- GTF %>% filter(type=="exon")
exons_gtf$exon_id <- paste0(exons_gtf$seqid,":",
                            exons_gtf$start - 1,"-" ,
                            exons_gtf$end,":",
                              exons_gtf$strand)

exons_conservation <- left_join(exons_gtf %>% dplyr::select(exon_id, transcript_id,gene_name),
                                phastcons$exons, by=c(exon_id="name"))

summary_exons_conservation <- exons_conservation %>% group_by(gene_name) %>% arrange(-mean_cov) %>%
  summarise(phastCons_most_conserved_exon = mean_cov[1],
            most_conserved_exon_id = paste(unique(exon_id[mean_cov==max(mean_cov)]),
                                           collapse = ","),
            phastCons_mean_exons = mean(mean_cov[!duplicated(exon_id)]))

# intron conservation ----
head(phastcons$introns)

intron_bed <- read.table("outputs/bed_files/LSK_StemLinc.combined.sorted.introns.bed")
intron_bed$intron_id <- paste0(intron_bed$V1,":",
                               intron_bed$V2,"-",
                               intron_bed$V3,":",
                               intron_bed$V6)
phastcons$introns <- left_join(phastcons$introns,intron_bed %>% dplyr::select(V4, intron_id),
                               by=c(name="V4"))

phastcons$introns <- separate(phastcons$introns,col = 1,into = c("transcript_name","intron_number"),
                              sep = "_(?!.*_)", remove = F)


phastcons$introns <- left_join(phastcons$introns, GTF %>% filter(type=="transcript") %>%
                                 dplyr::select(transcript_id,gene_name), by=c(transcript_name="transcript_id"))

summary_introns_conservation <- phastcons$introns %>% group_by(gene_name) %>%
  arrange(-mean_cov) %>% summarise(phastCons_most_conserved_intron = mean_cov[1],
                                   most_conserved_intron_id = paste(unique(intron_id[mean_cov==max(mean_cov)]),
                                                                    collapse = ","),
                                   phastCons_mean_introns = mean(mean_cov[!duplicated(intron_id)]))


# List of tables to merge
tables <- list(summary_exons_conservation,summary_introns_conservation,summary_merged_exons_conservation,summary_promoters_conservation)

# Merge all tables by "gene_name"
merged_table <- purrr::reduce(tables, full_join, by = "gene_name")

# Write table in database ----

# Specify the full path for the database file
db <- dbConnect(RSQLite::SQLite(), "/home/llorenzi/Rprojects/StemLinc_analyses/outputs/dbs/StemLinc.db")

dbExecute(db, "
    CREATE TABLE gene_conservation_phastCons35way (
        gene_name TEXT PRIMARY KEY,
        phastCons_most_conserved_exon REAL,
        most_conserved_exon_id TEXT,
        phastCons_mean_exons REAL,
        phastCons_most_conserved_intron REAL,
        most_conserved_intron_id TEXT,
        phastCons_mean_introns REAL,
        phastCons_most_conserved_merged_exon REAL,
        phastCons_mean_merged_exons REAL,
        phastCons_most_conserved_promoter REAL,
        most_conserved_promoter_transcript_id TEXT,
        phastCons_mean_promoters REAL,
        FOREIGN KEY (gene_name) REFERENCES gene_core_data(gene_name)
    );
")

dbWriteTable(
  conn = db,
  name = "gene_conservation_phastCons35way",
  value = merged_table,
  append = TRUE,  # Append data to the existing table
  row.names = FALSE
)

dbGetQuery(db, "SELECT * FROM gene_conservation_phastCons35way LIMIT 10;")

# Create metadata for the new table
new_metadata <- data.frame(
  table_name = "gene_conservation_phastCons35way",
  script_name = "parse_conservation_data.R",  # Replace with the actual script name
  level = "gene",
  filter_status = "none",
  analysis_version = "v1.0",
  timestamp = Sys.time(),
  source_data = "merged table from multiple phastCons conservation analysis (exons, merged exons, promoters, introns)",
  orig_scripts = "multiple: bigWigAverageOverBed was used with mm39.phastCons35way.bw over multiple bed files",  # Replace with your actual script(s)
  description = "Gene conservation data including exon, intron, merged exon, and promoter conservation metrics based on phastCons 35-way scores."
)

# Write the metadata to the database
dbWriteTable(db, "metadata", new_metadata, append = TRUE, row.names = FALSE)
