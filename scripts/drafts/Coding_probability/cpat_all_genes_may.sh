cd '/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/'

# cpat.py -x Mouse_Hexamer.tsv --antisense -d Mouse_logitModel.RData --top-orf=5 -g input_fa.fa -o result_folder/prefix

cpat.py -x '/home/llorenzi/references/CPAT_files/Mouse_Hexamer.tsv' --antisense -d '/home/llorenzi/references/CPAT_files/Mouse_logitModel.RData' --top-orf=5 -g 'data/various_fasta/merged_PCG_lncRNA_pseudo_transcriptome.240510.fa' -o 'data/CPAT_results/merged_PCG_lncRNA_pseudo_transcriptome.240510'


