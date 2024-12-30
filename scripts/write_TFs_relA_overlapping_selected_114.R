TF_ol=readRDS("outputs/overlap_marks/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.gene_classif.TFs_overlap_500_250.rds")
TF_ol <- TF_ol %>%filter(V5%in%gene_level_info_sel$gene_name)
TF_ol <- TF_ol %>% filter(grepl("_A",V11))
TF_ol <- TF_ol %>% select(V7:V13)
write.table(TF_ol,"outputs/overlap_marks/TFs_relA_overlapping_selected_114.bed",quote = F,col.names = F,
            row.names = F,sep = "\t")
