options(scipen = 999)
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
require(rtracklayer)
require(bedtoolsr)
require(trastools)

chr_order=paste0("chr",c(seq(1:22),"M","X","Y"))
get_near_peaks_unstranded <- function(bed,mark,w=1000,cols="2,3,8,9"){
  win=bt.window(a = bed,mark,w = w)
  win.ov=bt.overlap(win,cols=cols)
  return(win.ov)
} # this function finds overlap with the mark bed in a window of 1000 by default
# then it computes the amount of overlap (positive values) or distance (negative values)
get_near_peaks_stranded <- function(bed,mark,cols="2,3,8,9",lwin=1000,rwin=0,sm=T){
  win=bt.window(a = bed,mark,l = lwin,r = rwin,sw = T,sm=sm)
  win.ov=bt.overlap(win,cols=cols)
  return(win.ov)
} # this function finds overlap with the mark in a strand specific manner, e.g
# useful to compute overlap relative to TSS (e.g CAGE)
# or end of feature (e.g polyA signal), in case of TF, we want the window
# relative to the TSS but we do not care about
# strand, then sm=F

get_closest_peaks_stranded <- function(bed,mark,s=T){
  mark <- mark %>% arrange(match(V1,chr_order),V2)
  closest_mark <- bt.closest(bed, mark, D = "a",s = s ) # this function already retrieves the distance relative to transcript strand
  return(closest_mark)
} # this function is more general, it finds the closest mark whatsoever, even if it is far far away...
# With this approach we can plot distribution of distance for all features and later filter based
# on a user defined proximity cutoff


# documentation:
# https://bedtools.readthedocs.io/en/latest/content/tools/window.html
# https://bedtools.readthedocs.io/en/latest/content/tools/overlap.html
# https://bedtools.readthedocs.io/en/latest/content/tools/closest.html
# params:
#         sw Define -l and -r based on strand. For example if used, -l 500 for a negative-stranded feature will add 500 bp downstream. - Default = disabled.
#         sm  Only report hits in B that overlap A on the _same_ strand. - By default, overlaps are reported without respect to strand.
#         cols  start1,end1,start2,end2


## Load data---------------------------
source("scripts/source_all_functions.R")
outdir="outputs/overlap_marks/"
dir.create(outdir)

gene_level_data_path <- "outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.tsv"
check_file_header(gene_level_data_path)
gene_level_info <- read.table(gene_level_data_path,header = T)
gene_level_bed <- extract_bed(gene_level_info)

transcript_level_data_path <- "outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.gene_classif.tsv"
transcript_level_data <- read.table(transcript_level_data_path,header = T)
transcript_level_bed <- extract_bed(transcript_level_data,name = "V1",score = "gene_name")
TSS_bed <- transcript_level_bed %>% mutate(end=start+1)
TSS_bed <- TSS_bed %>% filter(!duplicated(TSS_bed%>%select(seqid,start,end,strand)))

gl_out_pref=gsub(".tsv","",basename(gene_level_data_path))
tl_out_pref=gsub(".tsv","",basename(transcript_level_data_path))
