## Head -------------------------------------
##
##
## Purpose of script: check overlap with repeat regions
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-12-17
##
## Email: lucialorenzi90@gmail.com
##
## Notes ---------------------------
##
##
##
##
## Setup ---------------------------

options(scipen = 999)
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
library(rtracklayer)
library(trastools)
require(bedtoolsr)
library(DBI)
library(RSQLite)
library(AnnotationHub)
#repeats=read.table("/home/llorenzi/work_local/UCSC_RepeatMasker_mm39.bed")
db <- dbConnect(RSQLite::SQLite(), "outputs/dbs/StemLinc.db")
GTF <- readGFF("data/raw/LSK_StemLinc.combined.sorted.gtf")
# annotation of the reference transcriptome used, a merge between GENCODE and RefSeq

transcript_level_info <- dbReadTable(db,"transcript_core_data")
GTF$gene_id <- transcript_level_info$gene_name[match(GTF$transcript_id,
                                                       transcript_level_info$transcript_id)]
gene_level_info <- dbReadTable(db,"gene_core_data")



# Initialize AnnotationHub
ah <- AnnotationHub()

# Search for RepeatMasker data for the mm39 genome
query_mm39 <- query(ah, c("mm39", "RepeatMasker"))

# List available resources
query_mm39

# Retrieve the RepeatMasker data
rmsk_mm39 <- query_mm39[["AH99013"]]

# Inspect the GRanges object
rmsk_mm39

# Convert RepeatMasker GRanges to a BED-like data frame
rmsk_bed <- data.frame(
  chrom = as.character(seqnames(rmsk_mm39)),     # Chromosome
  start = start(rmsk_mm39) - 1,                 # 0-based start
  end = end(rmsk_mm39),                         # End
  name = rmsk_mm39$repName,                     # Repeat name
  score = 0,                                    # Default score
  strand = as.character(strand(rmsk_mm39)),     # Strand
  repClass = rmsk_mm39$repClass,                # Repeat class
  repFamily = rmsk_mm39$repFamily               # Repeat family
)

# Replace NA strands with "."
rmsk_bed$strand[is.na(rmsk_bed$strand)] <- "."



get_overlap_with_repeats_from_GTF=function(gtf){

  pgme=trastools::get_per_gene_merged_exons(gtf)
  pgme=pgme[,c(3:5,2,6,7)]
  pgme=pgme[order(pgme$seqnames,pgme$start),]
  exonic_length=pgme%>%group_by(group_name)%>%
    summarise(exonic_length=sum(width))
  pgme$start <- pgme$start - 1 # 0-based for bed format
  ol=bt.intersect(a = pgme,b=rmsk_bed,wo = T)
  ol=ol%>% rowwise()%>% mutate(ol_st=max(c(V2,V8)),ol_en=min(c(V3,V9)))
  colnames(ol) <- c("chr", "start","end","gene_name","merged_exon_width","strand",
                                        "rep_chr","rep_start","rep_end","rep_name","score","rep_strand",
                                        "rep_class","rep_family","ol_len","ol_start","ol_end")
  olsm=ol%>%group_by(gene_name)%>%
    summarise(repeats=paste(unique(rep_name[order(-ol_len)]),
                            collapse = ","),
              repeats_classes=paste(unique(rep_class[order(-ol_len)]),
                                    collapse = ","),
              t.overlap_length=sum(IRanges::width(IRanges::reduce(IRanges(
                start = ol_start,end = ol_end)))))
  olsm=left_join(olsm,exonic_length,by=c(gene_name="group_name"))
  olsm$repeat.fraction=olsm$t.overlap_length/olsm$exonic_length

  # strand-specific
  olsm_stranded <- ol%>%group_by(gene_name,rep_strand)%>%
    summarise(repeats=paste(unique(rep_name[order(-ol_len)]),
                            collapse = ","),
              repeats_classes=paste(unique(rep_class[order(-ol_len)]),
                                    collapse = ","),
              t.overlap_length=sum(IRanges::width(IRanges::reduce(IRanges(
                start = ol_start,end = ol_end)))))

  olsm_stranded=left_join(olsm_stranded,exonic_length,by=c(gene_name="group_name"))
  olsm_stranded$repeat.fraction=olsm_stranded$t.overlap_length/olsm_stranded$exonic_length


  return(list(olsm,olsm_stranded,ol, exonic_length))
}


#

# calculate repeat overlap ----
ol_repeats=get_overlap_with_repeats_from_GTF(GTF)

# calculate the total fraction of bases covered by repeats ----
# per biotype


exonic_length <- ol_repeats[[4]]
exonic_length=left_join(exonic_length,
                        gene_level_info%>%dplyr::select(gene_name,
                                                        biotype),
                        by=c(group_name="gene_name"))

olsm=ol_repeats[[1]]
olsm <- left_join(olsm, gene_level_info %>%dplyr::select(gene_name,biotype))

olsm_per_biotype <- olsm %>% group_by(biotype) %>% summarise(
  covered_len=sum(t.overlap_length)
)

exonic_length_biot=exonic_length%>%group_by(biotype) %>%
  summarise(total_exonic_length=sum(exonic_length))

olsm_per_biotype$total_exonic_length=exonic_length_biot$total_exonic_length[match(olsm_per_biotype$biotype,
                                                                       exonic_length_biot$biotype)]

olsm_per_biotype$fraction_covered=olsm_per_biotype$covered_len/olsm_per_biotype$total_exonic_length
olsm_per_biotype$percent_covered=olsm_per_biotype$covered_len/olsm_per_biotype$total_exonic_length*100
olsm_per_biotype

# stranded
olsm_stranded=ol_repeats[[2]]
olsm_stranded <- left_join(olsm_stranded, gene_level_info %>%dplyr::select(gene_name,biotype,strand))

olsm_stranded_per_biotype <- olsm_stranded %>% group_by(biotype, sense = rep_strand==strand) %>% summarise(
  covered_len=sum(t.overlap_length)
)

olsm_stranded_per_biotype$total_exonic_length=exonic_length_biot$total_exonic_length[match(olsm_stranded_per_biotype$biotype,
                                                                                  exonic_length_biot$biotype)]

olsm_stranded_per_biotype$fraction_covered=olsm_stranded_per_biotype$covered_len/olsm_stranded_per_biotype$total_exonic_length
olsm_stranded_per_biotype$percent_covered=olsm_stranded_per_biotype$covered_len/olsm_stranded_per_biotype$total_exonic_length*100
olsm_stranded_per_biotype


olsm_stranded <- mutate(olsm_stranded, orientation = ifelse(rep_strand==strand,"sense","antisense"))

olsm_stranded_summ <- pivot_wider(olsm_stranded, id_cols = gene_name,names_from = orientation,
                             values_from = repeat.fraction, names_prefix = "repeat.fraction_")
## write tables in db ----

### all orientations ----
head(olsm)

gene_overlap_RepeatMasker_any_orientation <- olsm %>% dplyr::select(gene_name,repeats,
                                                                    repeats_classes,total_overlap_length = t.overlap_length,
                                                                    total_exonic_length = exonic_length,
                                                                    fraction_covered_by_repeats = repeat.fraction)


# Create the table in the database
dbExecute(db, "
    CREATE TABLE gene_overlap_RepeatMasker_any_orientation (
        gene_name TEXT PRIMARY KEY,
        repeats TEXT,
        repeats_classes TEXT,
        total_overlap_length INTEGER,
        total_exonic_length INTEGER,
        fraction_covered_by_repeats REAL,
        FOREIGN KEY (gene_name) REFERENCES gene_core_data (gene_name)
    );
")

# Populate the table using the original data frame name
dbWriteTable(
  conn = db,
  name = "gene_overlap_RepeatMasker_any_orientation",
  value = gene_overlap_RepeatMasker_any_orientation,
  append = TRUE,
  row.names = FALSE
)

# Insert metadata for the new table
dbExecute(db, "
    INSERT INTO metadata (
        table_name,
        script_name,
        level,
        filter_status,
        analysis_version,
        timestamp,
        source_data,
        orig_scripts,
        description
    )
    VALUES (
        'gene_overlap_RepeatMasker_any_orientation',
        'overlap_RepeatMasker_all_genes.R',
        'gene',
        'Only genes with overlap with any repeat in any orientation',
        'v1.0',
        DATETIME('now'),
        'transcript_core_data and repeatmasker_annotations_mm39',
        'overlap_repeat_analysis.R',
        'Overlap of genes with RepeatMasker annotations, any orientation.'
    );
")


### sense orientation ----

gene_overlap_RepeatMasker_sense_orientation <-  olsm_stranded %>% filter(orientation=="sense") %>%
  dplyr::select(gene_name,gene_strand = strand,  sense_repeats = repeats,sense_repeats_classes = repeats_classes,
                total_sense_overlap_length = t.overlap_length,
                total_exonic_length = exonic_length, fraction_covered_by_sense_repeats = repeat.fraction)

# Create the table in the database
dbExecute(db, "
    CREATE TABLE gene_overlap_RepeatMasker_sense_orientation (
        gene_name TEXT PRIMARY KEY,
        gene_strand TEXT,
        sense_repeats TEXT,
        sense_repeats_classes TEXT,
        total_sense_overlap_length INTEGER,
        total_exonic_length INTEGER,
        fraction_covered_by_sense_repeats REAL,
        FOREIGN KEY (gene_name) REFERENCES gene_core_data (gene_name)
    );
")


# Populate the table using the original data frame name
dbWriteTable(
  conn = db,
  name = "gene_overlap_RepeatMasker_sense_orientation",
  value = gene_overlap_RepeatMasker_sense_orientation,
  append = TRUE,
  row.names = FALSE
)



# Insert metadata for the new table
dbExecute(db, "
    INSERT INTO metadata (
        table_name,
        script_name,
        level,
        filter_status,
        analysis_version,
        timestamp,
        source_data,
        orig_scripts,
        description
    )
    VALUES (
        'gene_overlap_RepeatMasker_sense_orientation',
        'overlap_RepeatMasker_all_genes.R',
        'gene',
        'Only genes with overlap with any repeat in sense orientation',
        'v1.0',
        DATETIME('now'),
        'transcript_core_data and repeatmasker_annotations_mm39',
        'overlap_repeat_analysis.R',
        'Overlap of genes with RepeatMasker annotations, sense orientation.'
    );
")


### antisense orientation ----

gene_overlap_RepeatMasker_antisense_orientation <-  olsm_stranded %>% filter(orientation=="antisense") %>%
  dplyr::select(gene_name,gene_strand = strand,  antisense_repeats = repeats,antisense_repeats_classes = repeats_classes,
                total_antisense_overlap_length = t.overlap_length,
                total_exonic_length = exonic_length, fraction_covered_by_antisense_repeats = repeat.fraction)

# Create the table in the database
dbExecute(db, "
    CREATE TABLE gene_overlap_RepeatMasker_antisense_orientation (
        gene_name TEXT PRIMARY KEY,
        gene_strand TEXT,
        antisense_repeats TEXT,
        antisense_repeats_classes TEXT,
        total_antisense_overlap_length INTEGER,
        total_exonic_length INTEGER,
        fraction_covered_by_antisense_repeats REAL,
        FOREIGN KEY (gene_name) REFERENCES gene_core_data (gene_name)
    );
")


# Populate the table using the original data frame name

# first remove the * strand genes
gene_overlap_RepeatMasker_antisense_orientation <- gene_overlap_RepeatMasker_antisense_orientation%>%filter(gene_strand!="*")
dbWriteTable(
  conn = db,
  name = "gene_overlap_RepeatMasker_antisense_orientation",
  value = gene_overlap_RepeatMasker_antisense_orientation,
  append = TRUE,
  row.names = FALSE
)



# Insert metadata for the new table
dbExecute(db, "
    INSERT INTO metadata (
        table_name,
        script_name,
        level,
        filter_status,
        analysis_version,
        timestamp,
        source_data,
        orig_scripts,
        description
    )
    VALUES (
        'gene_overlap_RepeatMasker_antisense_orientation',
        'overlap_RepeatMasker_all_genes.R',
        'gene',
        'Only genes with overlap with any repeat in antisense orientation',
        'v1.0',
        DATETIME('now'),
        'transcript_core_data and repeatmasker_annotations_mm39',
        'overlap_repeat_analysis.R',
        'Overlap of genes with RepeatMasker annotations, antisense orientation.'
    );
")


#### still to process this part ----


overlap_repeat_classes <- ol_repeats[[2]]
colnames(overlap_repeat_classes) <- c("chr", "start","end","gene_name","merged_exon_width","strand",
                                      "rep_chr","rep_start","rep_end","rep_name","score","rep_strand",
                                      "rep_class","rep_family","ol_start","ol_end")
gene_level_info=left_join(gene_level_info,
                          ol_repeats%>%dplyr::select(gene_id,repeats,repeat.fraction))
gene_level_info$repeat.fraction[is.na(gene_level_info$repeat.fraction)]=0

# plots ----

# distribution of repeats in each biotype:
# percent of total bases covered by each repeat in each biotype ----
frac_per_biot_rep=olsm%>% group_by(V11,biotype)%>%
  summarise(total_overlap_length=sum(t.overlap_length))

frac_per_biot_rep$total_exonic_length=exonic_length_biot$total_exonic_length[match(
  frac_per_biot_rep$biotype,exonic_length_biot$biotype
)]


frac_per_biot_rep$percent_repeat=frac_per_biot_rep$total_overlap_length/frac_per_biot_rep$total_exonic_length*100

frac_per_biot_rep_wide=pivot_wider(frac_per_biot_rep,
                                   id_cols = biotype,
                                   names_from = V11,
                                   values_from = percent_repeat)
frac_per_biot_rep_wide[is.na(frac_per_biot_rep_wide)] <- 0

frac_per_biot_rep_wide_long=pivot_longer(frac_per_biot_rep_wide,
                                         cols = 2:ncol(frac_per_biot_rep_wide),
                                         names_to = "repeat_type",
                                         values_to = "percent_length_covered")

library(ggpubr)
ggplot(frac_per_biot_rep_wide_long,aes(x=biotype,y=percent_length_covered)) +
  geom_col() + facet_wrap( ~ repeat_type)

# remove low repeats and add all repeats ----
frac_per_biot_rep_wide=left_join(frac_per_biot_rep_wide,
                                 exonic_length_biot)
colnames(frac_per_biot_rep_wide)
colnames(frac_per_biot_rep_wide)[20]="all_repeats"

sel_cols=c("biotype","DNA","LINE","LTR","Simple_repeat","SINE","all_repeats")

tmpdat=frac_per_biot_rep_wide %>% dplyr::select(all_of(sel_cols))
frac_per_biot_rep_wide_long=pivot_longer(tmpdat,
                                         cols = 2:ncol(tmpdat),
                                         names_to = "repeat_type",
                                         values_to = "percent_length_covered")

### plot2 ----
g=ggplot(frac_per_biot_rep_wide_long,aes(x=biotype,y=percent_length_covered,fill=biotype)) +
  geom_col() + facet_wrap( ~ repeat_type) + scale_fill_nejm()

wdir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/analyses/repeat_overlap/plots/"


save_tiff_svg(plo = g,outdir = wdir,filename = "2_percent_length_covered_per_biotype_and_repeat_highclassesANDall_repeats")
# fraction of all loci with overlap with each repeat class ----
# get number of genes in each biotype
Ngenes_biot=exonic_length%>%group_by(simpl_biotype)%>%
  summarise(N_total_genes=n())
colnames(Ngenes_biot)[1]="biotype"

Ngenes_any_repeat=olsm%>%group_by(biotype)%>%
  summarise(N_genes_with_repeats=length(unique(V4)))

# N genes per repeat class per biotype:
Ngenes_rep_class_biot=olsm%>%group_by(V11,biotype)%>%
  summarise(Ngenes=length(unique(V4)))

Ngenes_rep_class_biot=left_join(Ngenes_rep_class_biot,
                                Ngenes_biot)

Ngenes_rep_class_biot$percent_Genes=Ngenes_rep_class_biot$Ngenes/Ngenes_rep_class_biot$N_total_genes*100

colnames(Ngenes_rep_class_biot)[1]="repeat_type"

g=ggplot(Ngenes_rep_class_biot,aes(x=biotype,y=percent_Genes,fill=biotype)) +
  geom_col() + facet_wrap( ~ repeat_type) + scale_fill_nejm()
g

### plot 3 ----
save_tiff_svg(plo = g,outdir = wdir,filename = "3_percent_Genes_overlap_repeat_class_per_biotype")

#
widerdat=pivot_wider(Ngenes_rep_class_biot,id_cols = biotype,
                     names_from = repeat_type,
                     values_from = percent_Genes)

widerdat[is.na(widerdat)] <- 0

Ngenes_biot=left_join(Ngenes_biot,Ngenes_any_repeat)
Ngenes_biot$any_repeat=Ngenes_biot$N_genes_with_repeats/Ngenes_biot$N_total_genes*100

widerdat=left_join(widerdat,Ngenes_biot)
sel_cols=c("biotype","DNA","LINE","Low_complexity","LTR","Simple_repeat","SINE","any_repeat")
tmpdat=widerdat%>%dplyr::select(all_of(sel_cols))
longdat=pivot_longer(tmpdat,
                     cols = 2:ncol(tmpdat),
                     names_to = "repeat_type",
                     values_to = "percent_Genes")
### plot 4----
g=ggplot(longdat,aes(x=biotype,y=percent_Genes,fill=biotype)) +
  geom_col() + facet_wrap( ~ repeat_type) + scale_fill_nejm()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
g
save_tiff_svg(plo = g,outdir = wdir,filename = "4_percent_Genes_overlap_repeat_class_per_biotype_highclasses_and_anyRepeat",h = 8,w = 10)
# do the same comparing mono and multi-exonic ----
