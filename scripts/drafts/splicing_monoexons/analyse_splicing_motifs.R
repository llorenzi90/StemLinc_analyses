# Load libraries and set paths to data ----
library(rtracklayer)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
source("scripts/source_all_functions.R")
#
  # gene annotation: biotypes,
  # gtf, intronic
  # intronic info with splicing sites
  # MaxEntScan results

# Define paths
gene_annot_path <- "outputs/gene_level_info_all_evidence_oct24.csv"
transcript_info_path <- "outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.gene_classif.tsv"
gtf_path <- "data/raw/LSK_StemLinc.combined.gtf"
intron_info_path <- "outputs/motifs/intronic_motifs/LSK_StemLinc.combined_intron_info/introns_analysis.csv"
ss5_path <- "outputs/motifs/intronic_motifs/LSK_StemLinc.combined_intron_info/MaxEntScan_5ss.txt"
ss3_path <- "outputs/motifs/intronic_motifs/LSK_StemLinc.combined_intron_info/MaxEntScan_3ss.txt"

# Define out dir for plots ----
outdir_plots="outputs/motifs/intronic_motifs/plots"
dir.create(outdir_plots)

# Load tables ----
gene_level_info <- read.csv(gene_annot_path)
transcript_level_info <- read.table(transcript_info_path,header = T)
gtf <- readGFF(gtf_path)
intron_info <- read.csv(intron_info_path)
MaxEnt_5 <- read.table(ss5_path)
colnames(MaxEnt_5) <- c("intron_id","ss5_seq","ss5_strength")
MaxEnt_3 <- read.table(ss3_path)
colnames(MaxEnt_3) <- c("intron_id","ss3_seq","ss3_strength")

colnames(transcript_level_info)[1]="transcript_id"
table(intron_info$transcript_id%in%transcript_level_info$transcript_id)
table(intron_info$transcript_id%in%gtf$transcript_id)

# Add MaxEntScores
intron_info <- intron_info %>% mutate(intron_id=paste(transcript_id,"intron",intron_number,sep = "_"))

intron_info <- left_join(intron_info,MaxEnt_5)
intron_info <- left_join(intron_info,MaxEnt_3)

# Filter to retain only transcripts that pass expression filter:
intron_info <- intron_info %>%filter(transcript_id%in%transcript_level_info$transcript_id)

# Add gene name and biotype
intron_info <- left_join(intron_info, transcript_level_info %>% dplyr::select(transcript_id, gene_name, gene_biotype))

# mutate donor and acceptor sites to upper case
intron_info <- intron_info %>% mutate(orig_donor_site=donor_site,
                                      orig_acceptor_site=acceptor_site_2,
                                      orig_acceptor_site3=acceptor_site_3,
                                      donor_site=toupper(donor_site),
                                      acceptor_site_2=toupper(acceptor_site_2),
                                      acceptor_site_3=toupper(acceptor_site_3))

# mutate intron number ----
intron_info <- intron_info %>% mutate(N_introns=ifelse(intron_number<10,intron_number,">10"))
# Filter to retain only genes in filtered set
table(gene_level_info$gene_name[gene_level_info$exonic_type=="multiexonic"]%in%intron_info$gene_name)

intron_info <- intron_info %>%filter(gene_name %in% gene_level_info$gene_name)

# Intron length distribution ----
g=ggplot(intron_info, aes(gene_biotype,y=intron_length, fill=gene_biotype)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_manual(values = nejm_pal) +
  theme_minimal()

g
svg(paste0(outdir_plots,"/intron_length_distribution_per_biotype.boxplot.svg"),height = 8, width = 10 )
print(g)
dev.off()


g=ggplot(intron_info, aes(x=intron_length, col=gene_biotype)) +
  geom_density() +
  scale_x_log10() +
  scale_color_manual(values = nejm_pal)+
  theme_minimal()

svg(paste0(outdir_plots,"/intron_length_distribution_per_biotype.density.svg"),height = 8, width = 10 )
print(g)
dev.off()


# Splice site strength ----

g=ggplot(intron_info, aes(gene_biotype,y=ss5_strength, fill=gene_biotype)) +
  geom_boxplot()  +
  scale_fill_manual(values = nejm_pal) +
  theme_minimal()+
  ggtitle("5' Splice site strength MaxEntScan")

g

g <- ggplot(intron_info, aes(x = gene_biotype, y = ss5_strength, fill = gene_biotype)) +
  geom_boxplot() +
  scale_fill_manual(values = nejm_pal) +
  theme_minimal() +
  ggtitle("5' Splice site strength MaxEntScan") +
  stat_compare_means(method = "kruskal.test", label.y = max(intron_info$ss5_strength, na.rm = TRUE) * 1.1) +  # Overall test
  stat_compare_means(
    comparisons = list(c("protein_coding", "lncRNA"), c("lncRNA", "potNovel"),
                       c("protein_coding","potNovel")),  # Pairwise comparisons
    method = "wilcox.test",
    label = "p.signif"
  )

svg(paste0(outdir_plots,"/5ss_strength.boxplot.svg"),height = 8, width = 10 )
print(g)
dev.off()

png(paste0(outdir_plots,"/5ss_strength.boxplot.png"),height = 800, width = 1000,res=150 )
print(g)
dev.off()

ggplot(intron_info, aes(gene_biotype,y=ss3_strength, fill=gene_biotype)) +
  geom_boxplot()  + scale_fill_manual(values = nejm_pal)


g <- ggplot(intron_info, aes(x = gene_biotype, y = ss3_strength, fill = gene_biotype)) +
  geom_boxplot() +
  scale_fill_manual(values = nejm_pal) +
  theme_minimal() +
  ggtitle("3' Splice site strength MaxEntScan") +
  stat_compare_means(method = "kruskal.test",
                     label.y = max(intron_info$ss3_strength, na.rm = TRUE) * 1.1) +  # Overall test
  stat_compare_means(
    comparisons = list(c("protein_coding", "lncRNA"), c("lncRNA", "potNovel"),
                       c("protein_coding","potNovel")),  # Pairwise comparisons
    method = "wilcox.test",
    label = "p.signif"
  )

svg(paste0(outdir_plots,"/3ss_strength.boxplot.svg"),height = 8, width = 10 )
print(g)
dev.off()

png(paste0(outdir_plots,"/3ss_strength.boxplot.png"),height = 800, width = 1000,res=150 )
print(g)
dev.off()
# strength per intron number
ggplot(intron_info, aes(x=N_introns,y=ss5_strength))+geom_boxplot()


ggplot(intron_info, aes(x=N_introns,y=ss3_strength))+geom_boxplot()

# Summary of splice sites strength at gene level ----
ss3_strengthGL <- intron_info %>% group_by(gene_name) %>%
  arrange(-ss3_strength) %>% summarise(max_ss3_MaxEntScan=ss3_strength[1],
                                       intron_id_maxss3= intron_id[1])
ss5_strengthGL <- intron_info %>% group_by(gene_name) %>%
  arrange(-ss5_strength) %>% summarise(max_ss5_MaxEntScan=ss5_strength[1],
                                       intron_id_maxss5 = intron_id[1])

gene_splice_sites_maxMaxEntScan <- left_join(ss5_strengthGL, ss3_strengthGL)

# Write into the database ----
library(RSQLite)
library(DBI)
db <- dbConnect(RSQLite::SQLite(), "outputs/dbs/StemLinc.db")

dbExecute(db, "
    CREATE TABLE gene_splice_sites_maxMaxEntScan (
        gene_name TEXT PRIMARY KEY,
        max_ss5_MaxEntScan REAL,
        intron_id_maxss5 TEXT,
        max_ss3_MaxEntScan REAL,
        intron_id_maxss3 TEXT,
        FOREIGN KEY (gene_name) REFERENCES gene_core_data(gene_name)
    );
")

dbWriteTable(
  conn = db,
  name = "gene_splice_sites_maxMaxEntScan",
  value = gene_splice_sites_maxMaxEntScan,
  append = TRUE,
  row.names = FALSE
)

# Preview the first rows
dbGetQuery(db, "SELECT * FROM gene_splice_sites_maxMaxEntScan LIMIT 10;")

metadata <- data.frame(
  table_name = "gene_splice_sites_maxMaxEntScan",
  script_name = "analyse_splicing_motifs.R",  # Replace with your script name
  level = "gene",
  filter_status = "3 replicates gene level, >= 0.2 TPM, only biotypes of interest, multiexonc genes (at least 1 intron)",
  analysis_version = "v1.0",
  timestamp = Sys.time(),
  source_data = "MaxEntScan results",
  orig_scripts = "run_MaxEntScan.sh",  # Replace with your actual script name
  description = "As a measure of splicing strength at gene level, for each gene, the introns with maximum MaxEntScan scores for splice sites (5' ss and 3' ss) are retrieved together with their scores"
)

# Append the metadata to the database
dbWriteTable(db, "metadata", metadata, append = TRUE, row.names = FALSE)


# Check frequences of donor and acceptor sites
table(toupper(intron_info$donor_site),intron_info$gene_biotype)
table(toupper(intron_info$acceptor_site_2),intron_info$gene_biotype)
table(toupper(intron_info$acceptor_site_3),intron_info$gene_biotype)
table(toupper(intron_info$donor_site),toupper(intron_info$acceptor_site_2))

# Intron characteristics
intron_info %>% group_by(gene_biotype) %>%
  summarise(N_genes=length(unique(gene_name)),
            N_transcripts=length(unique(transcript_id)),
            N_introns=n(),
            mean_introns_per_transcript=N_introns/N_transcripts
            )

intron_length_summary <- intron_info %>%
  group_by(gene_biotype) %>%
  summarise(
    mean_length = mean(intron_length, na.rm = TRUE),
    median_length = median(intron_length, na.rm = TRUE),
    min_length = min(intron_length, na.rm = TRUE),
    max_length = max(intron_length, na.rm = TRUE),
    sd_length = sd(intron_length, na.rm = TRUE),
    n_introns = n()  # Total number of introns
  )

# View the summary table
print(intron_length_summary)
# Frequences of donor and accpetor sites ----

## Individual sites ----
# List of site types to process
splice_sites <- c("donor_site", "acceptor_site_2", "acceptor_site_3")

# Define a function for heatmap creation
generate_heatmap <- function(site_type) {
  # Calculate frequencies and totals for the given site type
  site_frequencies <- intron_info %>%
    group_by(gene_biotype, .data[[site_type]]) %>%  # Dynamically reference the column
    summarise(frequency = n(), .groups = "drop") %>%
    group_by(gene_biotype) %>%
    mutate(total_introns = sum(frequency),
           scaled_frequency = (frequency / total_introns) * 100)

  # Create matrices for heatmap
  scaled_matrix <- site_frequencies %>%
    dplyr::select(gene_biotype, .data[[site_type]], scaled_frequency) %>%
    pivot_wider(names_from = .data[[site_type]], values_from = scaled_frequency, values_fill = 0) %>%
    column_to_rownames(var = "gene_biotype")

  count_matrix <- site_frequencies %>%
    dplyr::select(gene_biotype, .data[[site_type]], frequency) %>%
    pivot_wider(names_from = .data[[site_type]], values_from = frequency, values_fill = 0) %>%
    column_to_rownames(var = "gene_biotype")

  # Function to add text annotations
  text_annotation <- function(j, i, x, y, width, height, fill) {
    grid.text(count_matrix[i, j], x, y, gp = gpar(fontsize = 8, col = "black"))
  }

  # Create heatmap
  heatmap <- Heatmap(
    scaled_matrix,
    name = paste0("Scaled Frequency (%) - ", site_type),
    col = colorRamp2(c(0, max(scaled_matrix)), c("white", nejm_pal[4])),
    cell_fun = text_annotation,
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_title = "Gene Biotype",
    column_title = paste("Splice Site -", site_type)
  )

  # Save heatmap to file
  output_file <- paste0(outdir_plots,"/",site_type, "_heatmap.svg")
  svg(output_file, width = 10, height = 8)
  draw(heatmap)
  dev.off()

  cat("Heatmap saved to:", output_file, "\n")
}

# Loop over splice sites and generate heatmaps
for (site in splice_sites) {
  generate_heatmap(site)
}

## Donor x acceptor sites per biotype ----
# Function to generate heatmap for donor_site x acceptor_site_2 for a specific biotype
generate_combined_heatmap <- function(biotype, outdir_plots) {
  # Filter data for the given biotype
  biotype_data <- intron_info %>%
    filter(gene_biotype == biotype) %>%
    group_by(donor_site, acceptor_site_2) %>%
    summarise(frequency = n(), .groups = "drop") %>%
    mutate(total_combination = sum(frequency),
           scaled_frequency = (frequency / total_combination) * 100)

  # Create matrices for heatmap
  scaled_matrix <- biotype_data %>%
    dplyr::select(donor_site, acceptor_site_2, scaled_frequency) %>%
    pivot_wider(names_from = acceptor_site_2, values_from = scaled_frequency, values_fill = 0) %>%
    column_to_rownames(var = "donor_site")

  count_matrix <- biotype_data %>%
    dplyr::select(donor_site, acceptor_site_2, frequency) %>%
    pivot_wider(names_from = acceptor_site_2, values_from = frequency, values_fill = 0) %>%
    column_to_rownames(var = "donor_site")

  # Function to add text annotations
  text_annotation <- function(j, i, x, y, width, height, fill) {
    grid.text(count_matrix[i, j], x, y, gp = gpar(fontsize = 8, col = "black"))
  }

  # Create heatmap
  heatmap <- Heatmap(
    scaled_matrix,
    name = paste0("Scaled Frequency (%) - ", biotype),
    col = colorRamp2(c(0, max(scaled_matrix)), c("white", "#20854EFF")),
    cell_fun = text_annotation,
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_title = "Donor Site",
    column_title = "Acceptor Site (2 nt)"
  )

  # Save heatmap to file
  output_file <- paste0(outdir_plots, "/", biotype, "_donor_acceptor_heatmap.svg")
  svg(output_file, width = 10, height = 8)
  draw(heatmap)
  dev.off()

  cat("Heatmap saved to:", output_file, "\n")
}

# Get unique biotypes and generate heatmaps for each
unique_biotypes <- unique(intron_info$gene_biotype)

for (biotype in unique_biotypes) {
  generate_combined_heatmap(biotype, outdir_plots)
}


# Branchpoint analysis ----
## Count number of branchpoints per intron ----
intron_info$branchpoint_count <- sapply(strsplit(intron_info$branchpoints, ","), length)

intron_info%>% group_by(gene_biotype) %>%
  summarise(has_BP=sum(branchpoint_count>0)/n(),
            mean_BP_per_600bp3ss=mean(branchpoint_count))

summary(intron_info$intron_length[intron_info$branchpoint_count==0])

# Summary statistics
summary(intron_info$branchpoint_count)

# Plot branchpoint count distribution
g=ggplot(intron_info, aes(x = branchpoint_count)) +
  geom_histogram(binwidth = 1, fill = "blue", alpha = 0.7) +
  labs(title = "Branchpoint Count Distribution",
       x = "Number of Branchpoints", y = "Frequency")
g


svg(paste0(outdir_plots,"/branchpoint_distribution_all_introns.svg"),height = 8, width = 10 )
print(g)
dev.off()

svg(paste0(outdir_plots,"/branchpoint_distribution_introns_per_biotype.svg"),height = 8, width = 10 )
print(g+ facet_wrap(~gene_biotype, scales = "free_y"))
dev.off()

g=ggplot(intron_info, aes(x = branchpoint_count, y=intron_length)) +
  geom_point() + scale_y_log10()+
  labs(title = "Branchpoint Count vs Intron Length",
       x = "Number of Branchpoints", y = "Intron Length")
g


g=ggplot(intron_info, aes(x = branchpoint_count, y=intron_length)) +
  geom_point(alpha=0.6) + scale_y_log10()+
  labs(title = "Branchpoint Count vs Intron Length",
       x = "Number of Branchpoints", y = "Intron Length")+
  facet_wrap(~gene_biotype)
g

tiff(paste0(outdir_plots,"/branchpoint_vs_intron_length_per_biotype.tiff"),height = 800, width = 1000, res = 300 )
print(g)
dev.off()

## Distance to branchpoint ----

# Select closest, second and third closest
intron_info <- intron_info %>%
  mutate(
    branchpoints_split = strsplit(branchpoints, ","),
    closest_branchpoint = as.character(sapply(branchpoints_split, function(bp) tail(bp, 1))),       # Last element
    second_closest_branchpoint = as.character(sapply(branchpoints_split, function(bp) tail(bp, 2)[1])), # Second to last element
    third_closest_branchpoint = as.character(sapply(branchpoints_split, function(bp) tail(bp, 3)[1])),# Third to last element
    dist_split = strsplit(distances_to_3ss, ","),
    closest_dist = as.integer(sapply(dist_split,  function(bp) tail(bp, 1))),
    second_closest_dist = as.integer(sapply(dist_split, function(bp) tail(bp, 2)[1])), # Second to last element
    third_closest_dist = as.integer(sapply(dist_split, function(bp) tail(bp, 3)[1]))
  ) %>%
  dplyr::select(-branchpoints_split,-dist_split)  # Remove intermediate list column


# Plot distance distribution

g=ggplot(intron_info, aes(x = gene_biotype, y = closest_dist, fill = gene_biotype)) +
  geom_boxplot() + theme_minimal() + scale_fill_manual(values = nejm_pal)+
  labs(title = "Closest Branchpoint Distance to 3'SS", x = "Biotype", y = "Distance to 3ss (bp)")

g

svg(paste0(outdir_plots,"/closest_branchpoint_distance.boxplot.svg"),height = 8, width = 10 )
print(g)
dev.off()

tiff(paste0(outdir_plots,"/closest_branchpoint_distance.boxplot.tiff"),height = 800, width = 1000, res=150 )
print(g)
dev.off()

png(paste0(outdir_plots,"/closest_branchpoint_distance.boxplot.png"),height = 800, width = 1000, res=150 )
print(g)
dev.off()

intron_info %>% group_by(gene_biotype) %>%
  summarise(mean_dist_closest=mean(closest_dist,na.rm=T),
            median_dist_closest=median(closest_dist,na.rm=T),
            mean_dist_2nd_closest=mean(second_closest_dist, na.rm=T),
            median_dist_2nd_closest=mean(second_closest_dist,na.rm=T),
            mean_dist_3nd_closest=mean(third_closest_dist, na.rm=T),
            median_dist_3nd_closest=mean(third_closest_dist,na.rm=T))

intron_info %>% group_by(gene_biotype) %>%
  summarise(
            mean_dist_2nd_closest=mean(second_closest_dist, na.rm=T),
            median_dist_2nd_closest=median(second_closest_dist,na.rm=T),
            mean_dist_3nd_closest=mean(third_closest_dist, na.rm=T),
            median_dist_3nd_closest=median(third_closest_dist,na.rm=T))


g=ggplot(intron_info, aes(x = gene_biotype, y = second_closest_dist, fill = gene_biotype)) +
  geom_boxplot() + theme_minimal() + scale_fill_manual(values = nejm_pal)+
  labs(title = "2nd closest Branchpoint Distance to 3'SS", x = "Biotype", y = "Distance to 3ss (bp)")

g

svg(paste0(outdir_plots,"/second_closest_branchpoint_distance.boxplot.svg"),height = 8, width = 10 )
print(g)
dev.off()

tiff(paste0(outdir_plots,"/second_closest_branchpoint_distance.boxplot.tiff"),height = 800, width = 1000, res=150 )
print(g)
dev.off()

png(paste0(outdir_plots,"/second_closest_branchpoint_distance.boxplot.png"),height = 800, width = 1000, res=150 )
print(g)
dev.off()


g=ggplot(intron_info, aes(x = closest_dist)) +
  geom_histogram(binwidth = 1, fill = "blue", alpha = 0.7) +
  labs(title = "Closest Branchpoint distance to 3'ss",
       x = "Distance (bp)", y = "Frequency") + facet_wrap(~gene_biotype, scales = "free_y")
g

svg(paste0(outdir_plots,"/Closest_branchpoint_distance.histogram.svg"),height = 8, width = 10 )
print(g)
dev.off()

g + xlim(c(0,200))
svg(paste0(outdir_plots,"/Closest_branchpoint_distance.histogram.truncated200bp.svg"),height = 8, width = 10 )
print(g + xlim(c(0,200)))
dev.off()

intron_info <- intron_info %>% mutate(closest_dist2=ifelse(closest_dist>50,">50",closest_dist))

sort(table(intron_info$closest_dist2))

## Branchpoint sequence distribution ----
### closest branchpoint ----

g=ggplot(intron_info, aes(x=closest_branchpoint,fill=closest_branchpoint)) +
  geom_bar() + facet_wrap(~gene_biotype, scales="free_y") +
  theme(axis.text.x = element_text(angle=45, hjust=1))

g
svg(paste0(outdir_plots,"/Closest_branchpoint_sequence_frequency_per_biotype.svg"),height = 8, width = 10 )
print(g )
dev.off()

### Check distance distribution for TTTAT ----
g <- ggplot(intron_info %>% filter(closest_branchpoint=="TTTAT"),
            aes(x=closest_dist)) +geom_histogram() + xlab("Distance to 3'ss (bp)")+
  facet_wrap(~gene_biotype, scales="free_y") + ggtitle("Distance distribution of closest 'TTTAT' motif")

g

svg(paste0(outdir_plots,"/Closest_TTTAT_sequence_distance_per_biotype.svg"),height = 8, width = 10 )
print(g )
dev.off()

### second closest branchpoint sequence distribution ----
g=ggplot(intron_info %>%filter, aes(x=second_closest_branchpoint,fill=second_closest_branchpoint)) +
  geom_bar() + facet_wrap(~gene_biotype, scales="free_y") +
  theme(axis.text.x = element_text(angle=45, hjust=1))

g
svg(paste0(outdir_plots,"/Second_closest_branchpoint_sequence_frequency_per_biotype.svg"),height = 8, width = 10 )
print(g )
dev.off()

# Next analyses ----


## Sequence logos for closest branchpoints ----

## Sequence logos for all branchpoints ----

## Strength per branchpoint type ----

## Fraction of introns that have a branchpoint within X bases ----

## strength per exon number ----




# Branchpoints vs. intron length
ggplot(intron_info, aes(x = intron_length, y = branchpoint_count)) +
  geom_point(alpha = 0.6) +
  labs(title = "Branchpoints vs. Intron Length", x = "Intron Length (bp)", y = "Number of Branchpoints")




# Match intron id with gene, transcript id, biotype,

# Parameters of interest: gene biotype, splicing strength, intron length, exon number

# Things I want to analyse: intron length distribution, 5' and 3' ss strength, splicing acceptor and donor sites,
# distribution, logos?,
# number of branchpoints and distance to end,



