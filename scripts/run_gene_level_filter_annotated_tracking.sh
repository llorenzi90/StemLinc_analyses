# Run script to filter based on reproducibility and minimum TPM cutoff, using annotated tracking as input
Rscript scripts/gene_level_filter_annotated_tracking.R -f "outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.tsv" -n 3 -t 0.2 # timestamp: 20240930_173216
