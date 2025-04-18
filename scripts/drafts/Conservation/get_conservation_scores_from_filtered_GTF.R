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

bigWigAOB <- function(bW_path=bigWig_path,bed_path){
  bw_pref=gsub(".bw","",basename(bW_path))
  out_tab_path=gsub(".bed$",paste0(".",bw_pref,".tsv"),bed_path)

  system(paste("scripts/public_scripts/bigWigAverageOverBed",
               bW_path,
               bed_path,
               out_tab_path))
  ### Read results ----
  cons_res <- read.table(out_tab_path)
  colnames(cons_res) <- c("name","size","covered_size","sum_cov","mean_noncov0","mean_cov")

  cons_res$gene_name <- GTF$gene_id[match(cons_res$name,
                                          GTF$transcript_id)]
  write.table(cons_res,out_tab_path,row.names = F,quote = F,sep = "\t")

}
args <- commandArgs(trailingOnly = TRUE)

# Check if an input file is provided
if (length(args) == 0) {
  stop("Please provide the input files as a command-line argument.")
}

gtf_path=args[1] # "data/raw/LSK_StemLinc.combined.gtf"
#annot_path=args[2] #"outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.gene_classif.tsv"

out_pref=gsub(".gtf","",basename(gtf_path))

# Conservation at merged exons ----
### Generate a bed file of merged exonic regions for each gene ----
GTF=readGFF(gtf_path)
#annot=read.table(annot_path,header = T)

GTF <- GTF%>%  mutate(gene_id=gene_name)
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
bigWigAOB(bed_path = out_bed_path)

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
bigWig_path="data/references/conservation_bigWigs/mm39.phastCons35way.bw"
bw_pref=gsub(".bw","",basename(bigWig_path))
out_tab_path=gsub(".bed$",paste0(".",bw_pref,".tsv"),out_bed_path)

system(paste("scripts/public_scripts/bigWigAverageOverBed",
             bigWig_path,
             out_bed_path,
             out_tab_path))
### Read results ----
cons_res <- read.table(out_tab_path)
colnames(cons_res) <- c("name","size","covered_size","sum_cov","mean_noncov0","mean_cov")

cons_res$gene_name <- GTF$gene_id[match(cons_res$name,
                                          GTF$transcript_id)]
write.table(cons_res,out_tab_path,row.names = F,quote = F,sep = "\t")

# Conservation at exons and introns ----
unique_exons_bed <- paste0("outputs/bed_files/",gsub(".gtf",".exons.bed",basename(gtf_path)))

system(paste("scripts/unique_exons_from_gtf.sh", gtf_path, unique_exons_bed))
unique_exons <- read.table(unique_exons_bed)
unique_exons$V5 <- 1000
write.table(unique_exons,unique_exons_bed,quote = F,col.names = F,row.names = F,sep = "\t")
bigWigAOB(bed_path = unique_exons_bed)

introns_bed <- paste0("outputs/bed_files/",gsub(".gtf",".introns.bed",basename(gtf_path)))

system(paste("python scripts/introns_from_gtf.py",gtf_path,introns_bed))
bigWigAOB(bed_path = introns_bed)
