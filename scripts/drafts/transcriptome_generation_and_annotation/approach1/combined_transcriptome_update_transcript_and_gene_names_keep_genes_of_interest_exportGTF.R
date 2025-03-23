library(rtracklayer)
annotated_tracking <- read.table("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/filtered_LSK_T-cell_macrophage.combined_annotated_tracking.gene_names.tsv")
gtf_path <- "data/hpc_data/filtered_LSK_T-cell_macrophage.combined.gtf"
# Write GTF with updated gene_names and transcript names ----
GTF=readGFF(gtf_path)
GTF <- GTF %>% mutate(transcript_id=gsub("TCONS","STTRA",transcript_id),
                      gene_id=gsub("XLOC","STLOC",gene_id))
GTF$gene_name <- annotated_tracking$gene_name[match(GTF$transcript_id,
                                                    annotated_tracking$V1)]
# keep only genes of interest:
GTF <- GTF %>% filter(gene_name %in% gene_level_info$gene_name)
export(GTF,"outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/LSK_T-cell_macrophage.combined.updated_names.sel_biotypes.gtf", format = "gtf")
