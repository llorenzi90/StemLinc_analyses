gene_level_info_sel <- read.csv("outputs/gene_level_info_all_evidence_oct24.selected114.csv")

# select info ----
gene_level_info_sel_new_names <- gene_level_info_sel %>%
  dplyr::select(gene_name,new_gene_name,biotype,
                exonic_type,mean_corr_signature,max_Signature_gene,
                best_classif_to_PCG,classif,
                best_closest_PCG,dist_PCG,
                maxcorrPCG,TF_rel.lev_A,N_TFs,evidence_level,HSPC_vs_diff.logFC,HSPC_vs_prog.logFC,
                prog_vs_diff.logFC,Enhancer,Delas_enrichment,best_Luo_enrichment,seqnames,start,end,strand)

gene_level_info_sel_new_names <- gene_level_info_sel_new_names %>% arrange(desc(mean_corr_signature))
write.csv(gene_level_info_sel_new_names,"gene_level_info_114.csv")
