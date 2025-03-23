# create bed from annotated filtered tracking
# classes
#        biotypes: lncRNA, protein_coding, pseudogenes, potNovel
#        exonic type: monoexonic / multiexonic

# transcripts/genes

# potNovel divided by classcode
files.sources = list.files("scripts/functions/",full.names = T)
sapply(files.sources, source)
library(tidyverse)
library(trastools)

tr_path="outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.filtered.20240919_140928.tsv"

tr=read.table(tr_path,header = T)

gl_path="outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240919_140928.tsv"
gl=read.table(gl_path,header = T)

# output prefix
tr_out_pref=gsub(".tsv","",basename(tr_path))
gl_out_pref=gsub(".tsv","",basename(gl_path))


# pre-process gene level data
gl=gl %>% filter(pass_filter)
genomic_ranges <- get_genomic_range_by_gene(gtf_df = tr,by = "gene_name")
gl <- left_join(gl,genomic_ranges)
exonic_type <- tr%>% group_by(gene_name)%>%summarise(exonic_type=ifelse(any(Nexons>1),"multiexonic","monoexonic"))
gl <- left_join(gl, exonic_type)
max_Nexons <- tr %>% group_by(gene_name) %>%summarise(max_Nexons=max(Nexons))
gl <- left_join(gl,max_Nexons)



# transcript level ----
# keep only transcripts of interest
tr <- tr %>% filter(transfrag_class %in% c("protein_coding","pseudogene","TEC","lncRNA","i","u","x","r"))
tr$exonic_type=ifelse(tr$Nexons==1,"monoexonic","multiexonic")
# loop across transfrag class and write corresponding bed
for (cl in unique(tr$transfrag_class)) {
  bed <- extract_bed(tr%>%filter(transfrag_class==cl),name = "V1",score = "Nexons")
  write_bed(bed, paste0("outputs/bed_files/tr_level/",tr_out_pref,".",cl,"_transfrags.bed"))
  for (et in unique(tr$exonic_type)) {
    bed <- extract_bed(tr%>%filter(transfrag_class==cl&exonic_type==et),name = "V1",score = "Nexons")
    write_bed(bed, paste0("outputs/bed_files/tr_level/",tr_out_pref,".",cl,"_transfrags.",et,".bed"))
  }
}

## collapse all classcodes into potNovel ----
bed <- extract_bed(tr%>%filter(transfrag_class%in%non_overlapping_class_codes),name = "V1",score = "Nexons")
write_bed(bed, paste0("outputs/bed_files/tr_level/",tr_out_pref,".potNovel_transfrags.bed"))
for (et in unique(tr$exonic_type)) {
  bed <- extract_bed(tr%>%filter(transfrag_class%in%non_overlapping_class_codes&exonic_type==et),
                     name = "V1",score = "Nexons")
  write_bed(bed, paste0("outputs/bed_files/tr_level/",tr_out_pref,".potNovel_transfrags.",et,".bed"))
}


## collapse intergenic classcodes ----
bed <- extract_bed(tr%>%filter(transfrag_class%in%c("u","r")),name = "V1",score = "Nexons")
write_bed(bed, paste0("outputs/bed_files/tr_level/",tr_out_pref,".intergenic_transfrags.bed"))
for (et in unique(tr$exonic_type)) {
  bed <- extract_bed(tr%>%filter(transfrag_class%in%c("u","r")&exonic_type==et),
                     name = "V1",score = "Nexons")
  write_bed(bed, paste0("outputs/bed_files/tr_level/",tr_out_pref,".intergenic_transfrags.",et,".bed"))
}

## promoter regions ----
dir.create("outputs/bed_files/tr_level/promoter_regions_500_250")
for (cl in unique(tr$transfrag_class)) {
  bed <- extract_bed(tr%>%filter(transfrag_class==cl),name = "V1",score = "Nexons")
  bed <- get_promoters_from_bed(bed, before = 500, after = 250)
  write_bed(bed, paste0("outputs/bed_files/tr_level/promoter_regions_500_250/",tr_out_pref,".",cl,"_promoters.bed"))
  for (et in unique(tr$exonic_type)) {
    bed <- extract_bed(tr%>%filter(transfrag_class==cl&exonic_type==et),name = "V1",score = "Nexons")
    bed <- get_promoters_from_bed(bed, before = 500, after = 250)
    write_bed(bed, paste0("outputs/bed_files/tr_level/promoter_regions_500_250/",tr_out_pref,".",cl,"_promoters.",et,".bed"))
  }
}

### collapse all classcodes into potNovel ----
bed <- get_promoters_bed_from_annotated_tracking(tr%>%filter(transfrag_class%in%non_overlapping_class_codes),name = "V1",score = "Nexons")
write_bed(bed, paste0("outputs/bed_files/tr_level/promoter_regions_500_250/",tr_out_pref,".potNovel_promoters.bed"))
for (et in unique(tr$exonic_type)) {
  bed <- get_promoters_bed_from_annotated_tracking(tr%>%filter(transfrag_class%in%non_overlapping_class_codes&exonic_type==et),
                     name = "V1",score = "Nexons")
  write_bed(bed, paste0("outputs/bed_files/tr_level/promoter_regions_500_250/",tr_out_pref,".potNovel_promoters.",et,".bed"))
}


### collapse intergenic classcodes ----
bed <- get_promoters_bed_from_annotated_tracking(tr%>%filter(transfrag_class%in%c("u","r")),name = "V1",score = "Nexons")
write_bed(bed, paste0("outputs/bed_files/tr_level/promoter_regions_500_250/",tr_out_pref,".intergenic_promoters.bed"))
for (et in unique(tr$exonic_type)) {
  bed <- get_promoters_bed_from_annotated_tracking(tr%>%filter(transfrag_class%in%c("u","r")&exonic_type==et),
                     name = "V1",score = "Nexons")
  write_bed(bed, paste0("outputs/bed_files/tr_level/promoter_regions_500_250/",tr_out_pref,".intergenic_promoters.",et,".bed"))
}



# gene level ----
table(gl$gene_class)
gl <- gl %>% filter(gene_class %in% c("protein_coding","pseudogene","TEC","lncRNA","i","u","x","r"))


# loop across gene class and write corresponding bed
for (cl in unique(gl$gene_class)) {
  bed <- extract_bed(gl%>%filter(gene_class==cl),name = "gene_name",score = "max_Nexons")
  write_bed(bed, paste0("outputs/bed_files/gl_level/",gl_out_pref,".",cl,"_genes.bed"))
  for (et in unique(gl$exonic_type)) {
    bed <- extract_bed(gl%>%filter(gene_class==cl&exonic_type==et),name = "gene_name",score = "max_Nexons")
    write_bed(bed, paste0("outputs/bed_files/gl_level/",gl_out_pref,".",cl,"_genes.",et,".bed"))
  }
}

## collapse all classcodes into potNovel ----
bed <- extract_bed(gl%>%filter(gene_class%in%non_overlapping_class_codes),name = "gene_name",score = "max_Nexons")
write_bed(bed, paste0("outputs/bed_files/gl_level/",gl_out_pref,".potNovel_genes.bed"))
for (et in unique(gl$exonic_type)) {
  bed <- extract_bed(gl%>%filter(gene_class%in%non_overlapping_class_codes&exonic_type==et),
                     name = "gene_name",score = "max_Nexons")
  write_bed(bed, paste0("outputs/bed_files/gl_level/",gl_out_pref,".potNovel_genes.",et,".bed"))
}


## collapse intergenic classcodes ----
bed <- extract_bed(gl%>%filter(gene_class%in%c("u","r")),name = "gene_name",score = "max_Nexons")
write_bed(bed, paste0("outputs/bed_files/gl_level/",gl_out_pref,".intergenic_genes.bed"))
for (et in unique(gl$exonic_type)) {
  bed <- extract_bed(gl%>%filter(gene_class%in%c("u","r")&exonic_type==et),
                     name = "gene_name",score = "max_Nexons")
  write_bed(bed, paste0("outputs/bed_files/gl_level/",gl_out_pref,".intergenic_genes.",et,".bed"))
}


## promoter regions ----
dir.create("outputs/bed_files/gl_level/promoter_regions_500_250")
for (cl in unique(gl$gene_class)) {
  bed <- get_promoters_bed_from_annotated_tracking(gl%>%filter(gene_class==cl),name = "gene_name",score = "max_Nexons")
  write_bed(bed, paste0("outputs/bed_files/gl_level/promoter_regions_500_250/",gl_out_pref,".",cl,"_promoters.bed"))
  for (et in unique(gl$exonic_type)) {
    bed <- get_promoters_bed_from_annotated_tracking(gl%>%filter(gene_class==cl&exonic_type==et),name = "gene_name",score = "max_Nexons")
    write_bed(bed, paste0("outputs/bed_files/gl_level/promoter_regions_500_250/",gl_out_pref,".",cl,"_promoters.",et,".bed"))
  }
}

### collapse all classcodes into potNovel ----
bed <- get_promoters_bed_from_annotated_tracking(gl%>%filter(gene_class%in%non_overlapping_class_codes),name = "gene_name",score = "max_Nexons")
write_bed(bed, paste0("outputs/bed_files/gl_level/promoter_regions_500_250/",gl_out_pref,".potNovel_promoters.bed"))
for (et in unique(gl$exonic_type)) {
  bed <- get_promoters_bed_from_annotated_tracking(gl%>%filter(gene_class%in%non_overlapping_class_codes&exonic_type==et),
                     name = "gene_name",score = "max_Nexons")
  write_bed(bed, paste0("outputs/bed_files/gl_level/promoter_regions_500_250/",gl_out_pref,".potNovel_promoters.",et,".bed"))
}


### collapse intergenic classcodes ----
bed <- get_promoters_bed_from_annotated_tracking(gl%>%filter(gene_class%in%c("u","r")),name = "gene_name",score = "max_Nexons")
write_bed(bed, paste0("outputs/bed_files/gl_level/promoter_regions_500_250/",gl_out_pref,".intergenic_promoters.bed"))
for (et in unique(gl$exonic_type)) {
  bed <- get_promoters_bed_from_annotated_tracking(gl%>%filter(gene_class%in%c("u","r")&exonic_type==et),
                     name = "gene_name",score = "max_Nexons")
  write_bed(bed, paste0("outputs/bed_files/gl_level/promoter_regions_500_250/",gl_out_pref,".intergenic_promoters.",et,".bed"))
}
