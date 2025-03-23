combined_transcriptome=read.table("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/filtered_LSK_T-cell_macrophage.combined_annotated_tracking.gene_names.tsv",header = T)
table(combined_transcriptome$transfrag_class,combined_transcriptome$Nsamps)

gl=read.table("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/filtered_LSK_T-cell_macrophage.gene_level_filtered.20250117.tsv",header = T)

table(gl$biotype,gl$Nsamps)

library(rtracklayer)
library(tidyverse)
gtf=readGFF("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/LSK_T-cell_macrophage.combined.updated_names.sel_biotypes.gtf")
expression=read.table("~/Documentos/all_samples_analyses/kallisto/filtered_assembly_StemLinc_samples_kallisto.TPM.tsv",header = T)
expression$transcript_id=gsub("TCONS","STTRA",expression$transcript_id)
expression$gene_name=combined_transcriptome$gene_name[match(expression$transcript_id,
                                                            combined_transcriptome$V1)]

expression <- expression[!is.na(expression$gene_name),]

View(gl%>%filter(biotype=="potNovel",exonic_type=="monoexonic",best_cc=="u"))
View(gl%>%filter(biotype=="potNovel",exonic_type=="multiexonic",best_cc=="u",Nsamps==3))

colSums(expression[expression$gene_name=="XLOC_009890",2:10])
expression$mean_expression <- apply(expression[])

LSK_tracking=read.table("outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.tsv",header = T)

conservation <- read.delim("outputs/bed_files/LSK_T-cell_macrophage.combined.updated_names.sel_biotypes.merged_exons.mm39.phastCons35way.tsv")

correlation <- read.delim("outputs/correlation_analysis/corr_ncRNAs_vs_allPCG.txt")

ref_annot <- read.delim("~/Documentos/references/annotated_tracking_file.updated_gene_names.20250121.txt")
ref_annot_gl <- trastools::get_genomic_range_by_gene(ref_annot,by = "gene_name")
ref_gtf <- readGFF("~/Documentos/references/merged_refs.combined.gtf")
ref_gtf$gene_name <- ref_annot$gene_name[match(ref_gtf$transcript_id,
                                               ref_annot$gffC_transcript_id)]

ref_annot_Nexons <- ref_gtf %>% filter(type=="exon") %>% group_by(transcript_id) %>%
  summarise(Nexons=n())

write.table(ref_annot_gl,"~/Documentos/references/gene_level_mereged_ref.updated_gene_names.20250121.bed",
            quote = F,col.names = F,row.names = F,sep = "\t")

# multiexons consistent between reps
LSK_tracking_non_ol_ref <- LSK_tracking %>% filter(!overlap_Ref)
table(LSK_tracking_non_ol_ref$Nsamps>1,LSK_tracking_non_ol_ref$Nexons>1)
View(LSK_tracking_non_ol_ref%>%filter(Nsamps>1,Nexons>1))
LSK_tracking_non_ol_ref_multiexon_cons <- LSK_tracking_non_ol_ref %>%filter(Nsamps>1 & Nexons>1)
length(unique(LSK_tracking_non_ol_ref_multiexon_cons$V2))

merged_LSK_Reps <- readGFF("~/Documentos/merged_LSK_reps.m200.T0.1.gtf")
View(merged_LSK_Reps%>% filter(seqid=="chr5",start>=75889731, end<=75891362))
merged_LSK_Reps <- readGFF("~/Documentos/merged_LSK_reps.default.gtf")
View(merged_LSK_Reps%>% filter(seqid=="chr5",start>=75889731, end<=75891362))
merged_LSK_Reps <- readGFF("~/Documentos/merged_LSK_reps.i.gtf")
View(merged_LSK_Reps%>% filter(seqid=="chr5",start>=75889731, end<=75891362))
