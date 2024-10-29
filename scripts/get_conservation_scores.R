## Notes ---------------------------
## Which measure to use?
## The PhastCons score is a probability that each nucleotide
## belongs to a conserved element, whereas abs(phyloP) is the -log(p-value)
## under a null hypothesis of neutral evolution, and a negative sign indicates
## faster-than expected evolution, while positive values imply conservation.May 27, 2010
## I can try and compare both
##
## Setup ---------------------------

options(scipen = 999)
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
library(rtracklayer)
source("scripts/source_all_functions.R")
## Load data---------------------------
gtf_path="data/raw/LSK_StemLinc.combined.gtf"
annot_path="outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.gene_classif.tsv"

out_pref=gsub(".gene_classif|.tsv","",basename(annot_path))

# Conservation at merged exons ----
### Generate a bed file of merged exonic regions for each gene ----
GTF=readGFF(gtf_path)
annot=read.table(annot_path,header = T)

GTF <- GTF%>%filter(transcript_id%in%annot$V1) %>%
  mutate(gene_id=annot$gene_name[match(transcript_id,
                                       annot$V1)])
gtf_exons=GTF %>% filter(type=="exon")

txdb <-makeGRangesFromDataFrame(as.data.frame(gtf_exons),keep.extra.columns = T)

per_gene_merged_exons <- reduce(split(txdb,txdb$gene_id))
per_gene_merged_exons <- as.data.frame(per_gene_merged_exons)

per_gene_merged_exons$exon_name=get_unique_names(per_gene_merged_exons$group_name)

merged_exons_bed <- extract_bed(per_gene_merged_exons,name = "exon_name")

### Write bed ----
out_bed_path=paste0("outputs/bed_files/",
                    out_pref,
                    ".merged_exons.bed")
write.table(merged_exons_bed,
            out_bed_path,
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

### bigWigAverageOverBed ----
#  system run UCSC's  bigWigAverageOverBed
# usage:
#  bigWigAverageOverBed in.bw in.bed out.tab
bigWig_path="~/Descargas/mm39.phastCons35way.bw"
bw_pref=gsub(".bw","",basename(bigWig_path))
out_tab_path=gsub(".bed$",
                  paste0(".",bw_pref,".tsv"),out_bed_path)

system(paste("~/Descargas/bigWigAverageOverBed",
             bigWig_path,
             out_bed_path,
             out_tab_path))
### Read results ----
cons_res <- read.table(out_tab_path)
colnames(cons_res) <- c("name","size","covered_size","sum_cov","mean_noncov0","mean_cov")

cons_res$gene_name <- per_gene_merged_exons$group_name[match(cons_res$name,
                                                             per_gene_merged_exons$exon_name)]

# Conservation at promoter regions ----
### Generate a bed file of transcript promoter regions ----
bed <- extract_bed(GTF%>%filter(type=="transcript"),name = "transcript_id")
proms_bed <- get_promoters_from_bed(bed,before = 500, after = 250)

### Write bed ----
out_bed_path=paste0("outputs/bed_files/",
                    out_pref,
                    ".promoters.bed")
write.table(proms_bed,
            out_bed_path,
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

### bigWigAverageOverBed ----
#  system run UCSC's  bigWigAverageOverBed
# usage:
#  bigWigAverageOverBed in.bw in.bed out.tab
bigWig_path="~/Descargas/mm39.phastCons35way.bw"
bw_pref=gsub(".bw","",basename(bigWig_path))
out_tab_path=gsub(".bed$",paste0(".",bw_pref,".tsv"),out_bed_path)

system(paste("~/Descargas/bigWigAverageOverBed",
             bigWig_path,
             out_bed_path,
             out_tab_path))
### Read results ----
cons_res <- read.table(out_tab_path)
colnames(cons_res) <- c("name","size","covered_size","sum_cov","mean_noncov0","mean_cov")

cons_res$gene_name <- GTF$gene_id[match(cons_res$name,
                                          GTF$transcript_id)]

# Conservation at slices of 200 bp in exonic regions ----
