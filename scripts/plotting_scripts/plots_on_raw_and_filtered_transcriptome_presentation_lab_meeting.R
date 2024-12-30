raw_tr_path="data/raw/LSK_StemLinc.tracking"
tracking=read.table(raw_tr_path)

tracking <- tracking %>% mutate(transcript_class=ifelse(V4=="=",
                                                        "exact_match",
                                                        ifelse(V4%in%c("c","e","j","k","m","n","o"),
                                                               "novel_isoform",
                                                               ifelse(V4=="i",
                                                                      "intronic",
                                                                      "new_transcript"))))

N_Transcripts_per_class <- tracking %>% group_by(transcript_class) %>% summarise(N_transcripts=n() )

pie(N_Transcripts_per_class$N_transcripts,
    labels = paste0(N_Transcripts_per_class$transcript_class,
                    " - " ,
                    round((N_Transcripts_per_class$N_transcripts/
                      sum(N_Transcripts_per_class$N_transcripts))*100,1),"%"),
    col = c("#20854EFF", "#EE4C97FF","#0072B5FF","#E18727FF"))

#"#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF" "#FFDC91FF" "#EE4C97FF"

raw_gene_path="outputs/transcriptome_characterization/LSK_StemLinc.combined/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.tsv"
gene_level_raw=read.table(raw_gene_path,header = T)

gene_level_raw <- gene_level_raw %>% mutate(class=ifelse(best_cc=="=",
                                                         "exact_match",
                                                         ifelse(best_cc%in%c("c","e","j","k","m","n","o"),
                                                                "novel_isoform",
                                                                ifelse(best_cc=="i",
                                                                       "intronic",
                                                                       "new_transcript"))))

tracking$Nsamps=get_n_samples_per_transcript_from_tracking(tracking,cols = c(5:7))

tracking_class_Nsamps=tracking %>% select(transcript_class,Nsamps)
tracking_class_Nsamps$level="transcript"
colnames(tracking_class_Nsamps)[1]="class"

gene_class_Nsamps=gene_level_raw %>% select(class,Nsamps)
gene_class_Nsamps$level="gene"
colnames(gene_class_Nsamps)

class_Nsamps=rbind(tracking_class_Nsamps,
                   gene_class_Nsamps)

class_Nsamps$level=factor(class_Nsamps$level,levels = c("transcript","gene"))

g=ggplot(class_Nsamps,aes(x=class,fill=as.factor(Nsamps))) + geom_bar() +
  scale_fill_manual(values = nejm_pal) + facet_wrap(~level) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

g

save_tiff_svg(g,outdir ,"raw_transcriptome_transcript_gene_level_vs_Nsamples",w=8, h=5)


class_Nsamps %>% filter(Nsamps==3) %>% group_by(level,class) %>% summarise(n())

table(gene_level_raw$gene_class)

gene_level_raw_sel_biotypes <- gene_level_raw %>% filter(gene_class%in%c(biots_of_interest,"i","x","u","r"))

gene_level_raw_sel_biotypes <- gene_level_raw_sel_biotypes %>% mutate(overlapRef=factor(ifelse(class%in%c("exact_match","novel_isoform"),
                                                                                        "overlapRef","potNovel"),
                                                                                        levels=c("potNovel",
                                                                                                 "overlapRef")))

g=ggplot(gene_level_raw_sel_biotypes,
         aes(x=reorder(gene_class,max_mean_tpm,median) , y=max_mean_tpm,
             fill=reorder(gene_class,max_mean_tpm,median))) + geom_boxplot() +
  scale_y_log10() + theme_minimal()+
  facet_wrap(~overlapRef,scales="free_x") +labs(x="gene class", y="mean TPM",
                                                fill="gene class")

g

save_tiff_svg(g,outdir,"expression_raw_gene_level_per_interesting_classes")

g=ggplot(gene_level_raw_sel_biotypes,
       aes(x=reorder(gene_class,max_mean_tpm,median) , y=max_mean_tpm,
           fill=reorder(gene_class,max_mean_tpm,median))) + geom_boxplot() +
  scale_y_log10() + theme_minimal() + geom_hline(yintercept = 0.2,linetype = "dashed") +
  facet_wrap(~overlapRef,scales="free_x") +labs(x="gene class", y="mean TPM",
                                                fill="gene class")

g

save_tiff_svg(g,outdir,"expression_raw_gene_level_per_interesting_classes_dashed_0.2")

#
gene_level_info %>% ggplot(aes(x=biotype)) + geom_bar()
gene_level_info %>% group_by(biotype) %>% summarise(n())


N_genes_per_class <- gene_level_info %>% group_by(biotype) %>% summarise(N_genes=n() )

pie(N_genes_per_class$N_genes,
    labels = paste0(N_genes_per_class$biotype," - ",N_genes_per_class$N_genes,
                    " - " ,
                    round((N_genes_per_class$N_genes/
                             sum(N_genes_per_class$N_genes))*100,1),"%"),
    col = c("#7876B1FF","#E18727FF", "#0072B5FF" ,"#20854EFF","#EE4C97FF"))


#"#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF" "#FFDC91FF" "#EE4C97FF"

# nc gene classif and exonic type ----

g1 <- ggplot(gene_level_info%>% filter(biotype!="protein_coding"), aes(x=biotype, fill = classif)) +
  geom_bar(position = "dodge") + scale_fill_manual(values = nejm_pal) + facet_wrap(~exonic_type)
g1

g2 <- ggplot(gene_level_info%>% filter(biotype!="protein_coding"), aes(x=biotype, fill = exonic_type)) +
  geom_bar() + scale_fill_manual(values = nejm_pal) + facet_wrap(~classif , scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g2


outdir="outputs/plots/"
save_tiff_svg(g1,outdir = outdir, filename = "non_coding_gene_classification.simple.exonic_type", h = 8, w = 10)
save_tiff_svg(g2,outdir = outdir, filename = "non_coding_genes_exonic_type_per_biotype_facet_classif", h = 8, w = 10)


g3 <- ggplot(gene_level_info, aes(x=biotype,fill=exonic_type)) +geom_bar(position = "fill") +
  scale_fill_manual(values = nejm_pal)
g3
save_tiff_svg(g3,outdir = outdir, filename = "filtered_genes_exonic_type_per_biotype", h = 8, w = 10)


g4 <- ggplot(gene_level_info%>% filter(biotype!="protein_coding"), aes(x=classif, fill = exonic_type)) +
  geom_bar() + scale_fill_manual(values = nejm_pal) + facet_wrap(~biotype , scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g4
save_tiff_svg(g4,outdir = outdir, filename = "non_coding_genes_exonic_type_per_classif_facet_biotype", h = 8, w = 10)



g5 <- ggplot(gene_level_info%>% filter(biotype!="protein_coding"), aes(x=classif, fill = exonic_type)) +
  geom_bar(position = "fill") + scale_fill_manual(values = nejm_pal) + facet_wrap(~biotype , scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g5
save_tiff_svg(g5,outdir = outdir, filename = "non_coding_genes_exonic_type_per_classif_facet_biotype_fill", h = 8, w = 10)
