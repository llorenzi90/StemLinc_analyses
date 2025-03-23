fasta=outputs/fasta/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.promoters.fa

./scripts/get_markov_model_for_fimo.sh $fasta

fasta=outputs/fasta/LSK_StemLinc.combined.introns.fa

./scripts/get_markov_model_for_fimo.sh $fasta
