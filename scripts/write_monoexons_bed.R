library(rtracklayer)
library(DBI)
library(tidyverse)
db <- dbConnect(RSQLite::SQLite(), "outputs/dbs/StemLinc.db")


gtf=readGFF("data/raw/LSK_StemLinc.combined.sorted.gtf")

DBI::dbListTables(db)

transcript_annot <- dbReadTable(db,"transcript_annotation")

transcript_core <- dbReadTable(db,"transcript_core_data")

all_monoexons_bed <- transcript_core %>% filter(N_exons==1) %>%
                                 select(chr,start,end,transcript_id,gene_name,N_exons,strand)


all_monoexons_bed <- all_monoexons_bed %>% mutate(start = start - 1,
  name=paste(transcript_id,gene_name,sep = "_")) %>%
  select(chr, start, end, name, N_exons, strand)

write.table(all_monoexons_bed, "outputs/bed_files/tr_level/LSK_StemLinc.combined.sorted.all_monoexons.bed",
            sep = "\t", quote = F,col.names = F,row.names = F)




chr1_monoexons <- all_monoexons_bed %>% filter(chr =="chr1")
set.seed(1234)
chr1_monoexons_100_random <- chr1_monoexons[sample(1:nrow(chr1_monoexons),100),]
write.table(chr1_monoexons_100_random, "outputs/bed_files/tr_level/LSK_StemLinc.combined.sorted.chr1_monoexons.100_random.bed",
            sep = "\t", quote = F,col.names = F,row.names = F)

chr1_monoexons_20_random <- chr1_monoexons[sample(1:nrow(chr1_monoexons),20),]
write.table(chr1_monoexons_20_random, "outputs/bed_files/tr_level/LSK_StemLinc.combined.sorted.chr1_monoexons.20_random.bed",
            sep = "\t", quote = F,col.names = F,row.names = F)
export(gtf %>% filter(seqid == "chr1"),
       "/home/llorenzi/Documentos/test_data/gtf/LSK_StemLinc.combined.sorted.chr1.gtf", format = "GTF")


