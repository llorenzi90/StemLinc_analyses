
# Check that all input sequences have the expected length and defined nucleotides

python scripts/validate_maxentscan.py outputs/motifs/intronic_motifs/LSK_StemLinc.combined_intron_info/5ss_sequences.fasta 9
python scripts/validate_maxentscan.py outputs/motifs/intronic_motifs/LSK_StemLinc.combined_intron_info/3ss_sequences.fasta 23


# Run MaxEntScan
cd scripts/MaxEntScan

perl score3.pl ../../outputs/motifs/intronic_motifs/LSK_StemLinc.combined_intron_info/3ss_sequences.clean.fasta > ../../outputs/motifs/intronic_motifs/LSK_StemLinc.combined_intron_info/MaxEntScan_3ss.txt


perl score5.pl ../../outputs/motifs/intronic_motifs/LSK_StemLinc.combined_intron_info/5ss_sequences.clean.fasta > ../../outputs/motifs/intronic_motifs/LSK_StemLinc.combined_intron_info/MaxEntScan_5ss.txt

# Extract fasta sequences names to match results to sequences
grep ">" 3ss_sequences.clean.fasta |sed "s/^>//" > 3ss_sequences.clean.names.txt
grep ">" 5ss_sequences.clean.fasta |sed "s/^>//" > 5ss_sequences.clean.names.txt


# Merge both tables

paste 3ss_sequences.clean.names.txt MaxEntScan_3ss.txt > tmp.txt
mv tmp.txt MaxEntScan_3ss.txt

paste 5ss_sequences.clean.names.txt MaxEntScan_5ss.txt > tmp.txt
mv tmp.txt MaxEntScan_5ss.txt
