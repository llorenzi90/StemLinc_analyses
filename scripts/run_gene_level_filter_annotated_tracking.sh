# Run script to filter based on reproducibility and minimum TPM cutoff, using annotated tracking as input
Rscript scripts/gene_level_filter_annotated_tracking.R -f "outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.tsv" -n 3 -t 0.2 # timestamp: 20240930_173216

# T-cells
Rscript scripts/gene_level_filter_annotated_tracking.R -f "outputs/transcriptome_characterization/T_cell.combined/T_cell.combined_annotated_tracking.tsv" -n 3 -t 0.2 

Rscript scripts/classify_genes_relative_to_PCGs_v2.R "outputs/transcriptome_characterization/T_cell.combined/T_cell.combined_annotated_tracking.filtered.20250114_170649.tsv" "outputs/transcriptome_characterization/T_cell.combined/T_cell.combined_annotated_tracking.gene_level_info.20250114_170649.tsv"

Rscript scripts/filter_gtf_assign_gene_name_to_gene_id.R "data/raw/T_cell.combined.gtf" "outputs/transcriptome_characterization/T_cell.combined/T_cell.combined_annotated_tracking.filtered.20250114_170649.gene_classif.tsv" # second argument is the transcript level output from classify_genes_relative_to_PCGs...

# macrophages
Rscript scripts/gene_level_filter_annotated_tracking.R -f "outputs/transcriptome_characterization/Macro.combined/Macro.combined_annotated_tracking.tsv" -n 3 -t 0.2 

Rscript scripts/classify_genes_relative_to_PCGs_v2.R "outputs/transcriptome_characterization/Macro.combined/Macro.combined_annotated_tracking.filtered.20250114_170752.tsv" "outputs/transcriptome_characterization/Macro.combined/Macro.combined_annotated_tracking.gene_level_info.20250114_170752.tsv"

Rscript scripts/filter_gtf_assign_gene_name_to_gene_id.R "data/raw/Macro.combined.gtf" "outputs/transcriptome_characterization/Macro.combined/Macro.combined_annotated_tracking.filtered.20250114_170752.gene_classif.tsv"
