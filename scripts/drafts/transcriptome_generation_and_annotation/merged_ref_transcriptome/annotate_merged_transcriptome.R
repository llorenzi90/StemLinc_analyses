library(trastools)
library(tidyverse)
library(rtracklayer)
library(ggvenn)
merged_tracking_path <- "outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/StringTie_merge_per_sample_transcriptomes/gffcompare/StemLinc_vs_Ref.tracking"
merged_tmap_path <- "outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/StringTie_merge_per_sample_transcriptomes/gtfs/StemLinc_vs_Ref.StemLinc_filtered_merged_LSK_Macro_T.cell.gtf.tmap"
ref_annot_path <- "data/references/merged_refs_annotation/annotated_tracking_file.updated_gene_names.20250121.txt"
tracking <- read.table(merged_tracking_path)
tmap <- read.delim(merged_tmap_path)
merged_gtf_path <- "outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/StringTie_merge_per_sample_transcriptomes/gtfs/StemLinc_filtered_merged_LSK_Macro_T.cell.gtf"
path_track_merged_as_ref_vs_samples <- "outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/StringTie_merge_per_sample_transcriptomes/gffcompare/StemLinc_merged_gtfs_vs_meta.tracking"
ref_gtf_path <- "data/references/merged_refs_annotation/merged_refs.combined.updated_gene_names.20250121.gtf"
# Annotate tmap with ref annot ----
ref_annot <- read.delim(ref_annot_path)

tmap$ref_gene_name <- ref_annot$gene_name[match(tmap$ref_id,
                                                ref_annot$gffC_transcript_id)]
tmap$ref_biotype <- ref_annot$simplified_gene_biotype[match(tmap$ref_id,
                                                 ref_annot$gffC_transcript_id)]
tmap <- tmap %>% mutate(overlapRef = class_code %in% trastools::overlapping_class_codes,
                        gene_name = ifelse(overlapRef,ref_gene_name,qry_gene_id),
                        biotype = ifelse(overlapRef,ref_biotype,"potNovel"))
table(tmap$biotype,tmap$num_exons)

biots=c("protein_coding","lncRNA","potNovel","TEC","pseudogene")
length(unique(tmap$gene_name[tmap$biotype%in%biots]))


# Add coordinates ----
gtf <- readGFF(merged_gtf_path)
# Annotate tmap with each per-sample transcriptome ----
tracking <- read.table(path_track_merged_as_ref_vs_samples)
tracking <- separate(tracking,col = 3,into = c("ST_gene","transcript_id"),remove = F,sep = "\\|")

missed_genes <- unique(tmap$qry_gene_id)[!unique(tmap$qry_gene_id)%in%tracking$ST_gene]
View(tmap%>%filter(qry_gene_id%in%missed_genes))
View(gtf%>%filter(gene_id%in%missed_genes))
# All these are genes with "*" strand
table(gtf$strand=="*")
# It is OK to discard those

# Summarise tracking ----
# for each merged transcript I want to know: in which samples it was found
# best classcode, samples of best classcode, transcripts (ids) merged per sample, N transcripts merged in total, N samples merged
# first parse the tracking:
tracking <- trastools::split_samples_info(tracking,cols = 7:9,qnames = c("LSK","macrophage","T.cell"),remove = F )
tracking <- tracking %>%filter(V4%in%trastools::overlapping_class_codes)
ordered_classcodes <- c("=", "j", "k", "c", "m", "n", "e", "o", "p", "s", "x", "i", "y", "u", "r", ".")

tracking_summ_per_cc <- tracking %>% group_by(transcript_id,V4) %>% arrange(match(V4,ordered_classcodes)) %>%
  summarise(LSK=any(LSK_gene!="-"),
            macrophage=any(macrophage_gene!="-"),
            T.cell=any(T.cell_gene!="-"),
            LSK_tr=paste(unique(na.omit(LSK_transcript)),collapse = ","),
            macrophage_tr=paste(unique(na.omit(macrophage_transcript)),collapse = ","),
            T.cell_tr=paste(unique(na.omit(T.cell_transcript)),collapse = ","),
            N_merged_transcripts=sum(LSK_gene!="-") + sum(macrophage_gene!="-") + sum(T.cell_gene!="-")
            )
#
tracking_summ <- tracking_summ_per_cc %>% group_by(transcript_id) %>%
  summarise(best_cc=V4[1],
            all_ccs=paste(unique(V4),collapse = ","),
            LSK_best_cc=any(LSK[V4==best_cc]),
            macrophage_best_cc=any(macrophage[V4==best_cc]),
            T.cell_best_cc=any(T.cell[V4==best_cc]),
            LSK=any(LSK),
            macrophage=any(macrophage),
            T.cell=any(T.cell),
            LSK_tr=paste(LSK_tr[LSK_tr!=""],collapse = ","),
            macrophage_tr=paste(macrophage_tr[macrophage_tr!=""],collapse = ","),
            T.cell_tr=paste(T.cell_tr[T.cell_tr!=""],collapse = ","),
            N_merged_transcripts_best_cc=sum(N_merged_transcripts[V4==best_cc]),
            total_merged_transcripts=sum(N_merged_transcripts)
            )

table(tracking_summ$LSK,tracking_summ$macrophage,tracking_summ$T.cell)

# later I can further refine if I want more info...
# add transcript coordinates
tracking_summ <- left_join(tracking_summ,gtf%>%filter(type=="transcript") %>%dplyr::select(transcript_id=transcript_id,seqid,start,end,strand))
# add transcript info, including gene_name, nexons, reference info
tracking_summ <- left_join(tracking_summ, tmap %>% dplyr::select(transcript_id=qry_id, gene_name, biotype, ST_gene=qry_gene_id, ref_cc=class_code,num_exons,ref_gene_name,ref_biotype))

# Add number of exons in ref
ref_gtf <- readGFF(ref_gtf_path)
ref_gtf_exons <- ref_gtf %>% filter(type=="exon") %>% group_by(transcript_id) %>%summarise(N_exons=n())
tracking_summ <- left_join(tracking_summ, tmap %>% dplyr::select(transcript_id=qry_id,ref_id))
tracking_summ$ref_num_exons <- ref_gtf_exons$N_exons[match(tracking_summ$ref_id,ref_gtf_exons$transcript_id)]
# NOTE: with stringite --merge some transcripts get too long, that's why there are many "k" transcripts, for example StemLinc.10024.1

# Summarise at gene level ----
gtf <- gtf %>% filter(transcript_id%in%tracking_summ$transcript_id)
gtf$gene_id <- tracking_summ$gene_name[match(gtf$transcript_id,tracking_summ$transcript_id)]
gene_Exonic_length <- trastools::compute_gene_non_redundant_exonic_length(gtf)
gtf_exons <- gtf%>%filter(type=="exon") %>% mutate(exon_id=paste0(seqid,":",start,"-",end,":",strand))
total_exons_per_gene <- gtf_exons %>% group_by(gene_name=gene_id) %>% summarise(unique_exons=length(unique(exon_id)))
gene_level_annot <- tracking_summ %>% group_by(gene_name) %>%  arrange(match(ref_cc,ordered_classcodes)) %>%
  summarise(chr=seqid[1],
            start = min(start),
            end=max(end),
            strand=strand[1],
            best_cc=ref_cc[1],
            LSK=any(LSK),
            macrophage=any(macrophage),
            T.cell=any(T.cell),
            biotype=biotype[1],
            exonic_type=ifelse(any(num_exons>1),"multiexonic","monoexonic"),
            max_num_exons=max(num_exons),
            max_ref_num_exons=ifelse(best_cc%in%overlapping_class_codes,max(ref_num_exons),NA)

            )

gene_level_annot <- left_join(gene_level_annot, gene_Exonic_length %>% select(gene_name=gene_id,exonic_length))
gene_level_annot <- left_join(gene_level_annot, total_exons_per_gene )

table(gene_level_annot$biotype,gene_level_annot$exonic_type,gene_level_annot$LSK)
table(gene_level_annot$biotype,gene_level_annot$exonic_type,gene_level_annot$macrophage)
table(gene_level_annot$biotype,gene_level_annot$exonic_type,gene_level_annot$T.cell)



plot_three_venn <- function(traGL,
                            samp1_var,
                            samp2_var,
                            samp3_var,
                            feat_var,
                            filter_cond=rep(T,nrow(traGL)), # defaults the entire data
                            sampnames=c(sample_names,"Ref"),
                            title=""){
  traGL <- traGL %>% filter(filter_cond)
  if(is.logical(traGL%>%pull({{samp1_var}}))){
    vennlist <- list(traGL %>% filter({{samp1_var}}) %>% pull({{feat_var}}) ,
                     traGL %>% filter({{samp2_var}}) %>% pull({{feat_var}}),
                     traGL %>% filter({{samp3_var}}) %>% pull({{feat_var}}))
  }else{
    vennlist <- list(traGL %>% filter({{samp1_var}}!="-") %>% pull({{feat_var}}) ,
                     traGL %>% filter({{samp2_var}}!="-") %>% pull({{feat_var}}),
                     traGL %>% filter({{samp3_var}}!="-") %>% pull({{feat_var}}))
  }

  names(vennlist)=sampnames

  g=ggvenn(
    vennlist,
    fill_color = c("#0072B5FF","#FFDC91FF","#20854EFF"),
    stroke_size = 0.5, set_name_size = 4
  ) +ggtitle(title)
  print(g)
  #return(list(vennlist,g))
} # makes venn diagram of three sets of transcripts

gene_level_info_venn <- gene_level_annot %>%select(gene_name,biotype,LSK,macrophage,T.cell)

plot_three_venn(traGL = gene_level_info_venn,filter_cond = gene_level_info_venn$biotype=="potNovel",
                LSK,macrophage,T.cell,gene_name,sampnames = c("LSK","macrophage","T.cell"),
                title = "potNovel")

for (biot in biots) {
  plot_three_venn(traGL = gene_level_info_venn,filter_cond = gene_level_info_venn$biotype==biot,
                  LSK,macrophage,T.cell,gene_name,sampnames = c("LSK","macrophage","T.cell"),title = biot)


}


# write tables ----
## transcript level annot ----
dir.create("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation")
write.table(tracking_summ,
            "outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation/all_merged_transcripts_annotation.tsv",
            quote = F,row.names = F,sep = "\t")

write.table(tracking_summ %>% filter(biotype%in%biots),
            "outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation/all_merged_transcripts_annotation.sel_biotypes.tsv",
            quote = F,row.names = F,sep = "\t")

## gene_level annot ----
write.table(gene_level_annot,"outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation/all_merged_genes_annotation.tsv",
            quote = F,row.names = F,sep = "\t")
write.table(gene_level_annot%>%filter(biotype%in%biots),"outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation/all_merged_genes_annotation.sel_biotypes.tsv",
            quote = F,row.names = F,sep = "\t")
## gtf with gene_names ----
export(gtf,"outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/StringTie_merge_per_sample_transcriptomes/gtfs/StemLinc_filtered_merged_LSK_Macro_T.cell.gene_names.gtf",
       format = "gtf")
# write plots ----
pdf("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation/Venn_diagrams_samples_per_biotype_gene_level.pdf")
for (biot in biots) {
  plot_three_venn(traGL = gene_level_info_venn,filter_cond = gene_level_info_venn$biotype==biot,
                  LSK,macrophage,T.cell,gene_name,sampnames = c("LSK","macrophage","T.cell"),title = biot)


}
dev.off()
