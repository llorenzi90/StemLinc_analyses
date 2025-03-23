library(rtracklayer)
ref_gtf=readGFF('/home/llorenzi/Rprojects/StemLinc_analyses/data/references/merged_refs_annotation/merged_refs.combined.updated_gene_names.20250121.gtf')
ref_gtf$gene_id=ref_gtf$gene_name

export(ref_gtf,
       '/home/llorenzi/Rprojects/StemLinc_analyses/data/references/merged_refs_annotation/merged_refs.combined.updated_gene_names.20250121.gene_names.gtf',
       format = "gtf")

canonical_chrs=paste0("chr",c(1:19,"M","Y","X"))
gtf_withref=readGFF("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/StringTie_merge_per_sample_transcriptomes/gtfs/StemLinc_filtered_merged_LSK_Macro_T.cell.withRef.gtf")
export(gtf_withref%>%filter(seqid%in%canonical_chrs),"outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/StringTie_merge_per_sample_transcriptomes/gtfs/StemLinc_filtered_merged_LSK_Macro_T.cell.withRef.gtf")
