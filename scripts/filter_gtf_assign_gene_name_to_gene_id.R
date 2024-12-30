library(rtracklayer)
gtf_path="data/raw/LSK_StemLinc.combined.gtf"
annot_path="outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.filtered.20240930_173216.gene_classif.tsv"
timestamp <- gsub(".*?(([0-9]{8}_[0-9]{6})(_([0-9]{3}))?).*?$", "\\2.\\4", annot_path)

# read raw sample gtf
GTF=readGFF(gtf_path)
# read filtered annot to select filtered transcripts:
annot=read.table(annot_path,header = T)

# filter GTF to retain only transcripts in the filtered set
# and assign gene_names to gene_id from the annot file
GTF <- GTF%>%filter(transcript_id%in%annot$V1) %>%
  mutate(gene_id=annot$gene_name[match(transcript_id,
                                       annot$V1)])
out_path=gsub(".gtf",paste0(".",timestamp,"filtered.gene_name.gtf"),basename(gtf_path))
outdir="outputs/gtf/"
dir.create(outdir)
export(GTF,con = paste0(outdir,out_path),format = "gtf")
