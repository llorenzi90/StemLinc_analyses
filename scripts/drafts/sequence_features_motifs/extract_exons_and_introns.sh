cd /Rprojects/StemLinc_analyses/outputs/bed_files


# extract unique exons
~/Rprojects/StemLinc_analyses/scripts/unique_exons_from_gtf.sh ~/Rprojects/StemLinc_analyses/data/raw/LSK_StemLinc.combined.gtf

# extract introns
python ~/Rprojects/StemLinc_analyses/scripts/introns_from_gtf.py ~/Rprojects/StemLinc_analyses/data/raw/LSK_StemLinc.combined.gtf > ~/Rprojects/StemLinc_analyses/logs/introns_from_gtf.log 2>&1

