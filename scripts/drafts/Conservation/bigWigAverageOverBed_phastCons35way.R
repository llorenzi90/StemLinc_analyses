### bigWigAverageOverBed ----
#  system run UCSC's  bigWigAverageOverBed
# usage:
#  bigWigAverageOverBed in.bw in.bed out.tab
input_bed=commandArgs(trailingOnly = T)[1]
bigWig_path="data/references/conservation_bigWigs/mm39.phastCons35way.bw"
bw_pref=gsub(".bw","",basename(bigWig_path))
out_tab_path=gsub(".bed$",
                  paste0(".",bw_pref,".tsv"),input_bed)

system(paste("scripts/public_scripts/bigWigAverageOverBed",
             bigWig_path,
             input_bed,
             out_tab_path))

# run from project directory ~/Rprojects/StemLinc_analyses
