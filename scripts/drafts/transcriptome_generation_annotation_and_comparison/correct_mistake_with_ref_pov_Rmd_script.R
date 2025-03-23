# first split the tracking ref col in ref_gene_id and ref_transcript_id
tracking=separate(tracking,col = 3,
                  into = c("ref_gene_id","ref_transcript_id"),
                  sep = "\\|")

### add ref gene name
tracking$ref_gene_name=ref_annot$gene_name[match(tracking$ref_transcript_id,
                                                 ref_annot$V1)]


tracking$ref_biotype=ref_annot$simplified_gene_biotype[match(tracking$ref_transcript_id,
                                                             ref_annot$V1)]

# define priority order of classcodes
order_cc=c("=","j","k","c","m","n","e","o","p","s","x","i","y","r","u",".")

# get info of all reference transcripts
ref_transcripts_info=tchars_refGTF$EPT

# add info of overlapping transfrags
ref_transcripts_info <- left_join(ref_transcripts_info,
                                  tracking %>% filter(overlap_Ref) %>% dplyr::select(V1,V4,ref_transcript_id,Nexons),
                                  by=c(transcript_id="ref_transcript_id"))

# sort by classcode priority to select the best matching transfrag
ref_transcripts_info <- ref_transcripts_info%>%arrange(match(V4,order_cc))

ref_transcripts_info <- ref_transcripts_info%>%filter(!duplicated(transcript_id))

# add gene biotype
ref_transcripts_info$simpl_gene_biotype=ref_annot$simplified_gene_biotype[
  match(ref_transcripts_info$gene_id,
        ref_annot$gene_name)]
colnames(ref_transcripts_info)[4:6]=c("best_transfrag","transfrag_classcode","N_exons_transfrag")

ref_genes_info_old=ref_transcripts_info%>%group_by(gene_id)%>%
  summarise(multiexonic=ifelse(any(N_exons!=1),"Ref_gene_multiexonic","Ref_gene_monoexonic"),
            best_transfrag_cc=transfrag_classcode[1],
            simpl_gene_biotype=simpl_gene_biotype[1],
            transfrag_exonic_type=ifelse(N_exons_transfrag[1]!=1,"multiexonic","monoexonic"))


ref_genes_info=ref_transcripts_info%>%group_by(gene_id)%>%
  summarise(multiexonic=ifelse(any(N_exons!=1),"Ref_gene_multiexonic","Ref_gene_monoexonic"),
            best_transfrag_cc=transfrag_classcode[1],
            simpl_gene_biotype=simpl_gene_biotype[1],
            transfrag_exonic_type=ifelse(!all(is.na(N_exons_transfrag)),
                                         ifelse(any(N_exons_transfrag!=1,na.rm = T),
                                                "multiexonic","monoexonic"),NA))

table(ref_genes_info$multiexonic)/nrow(ref_genes_info)

ref_transcripts_info$gene_multiexonic=ref_genes_info$multiexonic[
  match(ref_transcripts_info$gene_id,
        ref_genes_info$gene_id)
]


## plots
ggplot(ref_transcripts_info,aes(x=simpl_gene_biotype,fill=gene_multiexonic)) +
  geom_bar() +theme_classic() +ylab("# transcripts") + ggtitle("Biotype distribution of reference transcripts and exonic type") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(ref_transcripts_info%>%filter(gene_multiexonic=="Ref_gene_multiexonic"),aes(x=simpl_gene_biotype,fill=N_exons>1)) +
  geom_bar(position = "fill") +theme_classic() +ylab("# transcripts") + ggtitle("Fraction of monoexonic transcripts in multiexonic reference genes")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(ref_transcripts_info,aes(x=gene_multiexonic,fill=transfrag_classcode)) +
  geom_bar(position = "fill") +theme_classic() +ylab("# transcripts") + ggtitle("Classcode distribution of transfrags overlapping reference multiexonic and monoexonic genes") + xlab("Ref_exonic_type")

ggplot(ref_transcripts_info,aes(x=transfrag_classcode,fill=simpl_gene_biotype)) +
  geom_bar() +
  theme_classic() +ylab("# Ref transcripts") + ggtitle("Reference transcripts recovered or not (NA) by transfrags") + facet_wrap(~ gene_multiexonic)

# ggplot(ref_transcripts_info,aes(x=simpl_gene_biotype,fill=transfrag_classcode )) +
#   geom_bar() +
#   theme_classic() +ylab("# transcripts") + ggtitle("All reference genes") + facet_wrap(~ gene_multiexonic) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#


ggplot(ref_transcripts_info,aes(x=simpl_gene_biotype,fill=transfrag_classcode )) +
  geom_bar(position="fill") +
  theme_classic() +ylab("fraction transcripts") + ggtitle("Fraction of ref transcripts overlapped by transfrags") + facet_wrap(~ gene_multiexonic) + theme(axis.text.x = element_text(angle = 45, hjust = 1))


#same but at gene level

ggplot(ref_genes_info,aes(x=simpl_gene_biotype,fill=best_transfrag_cc )) +
  geom_bar() +
  theme_classic() +ylab("# genes") + ggtitle("Fraction of reference genes overlapped by transfrags") + facet_wrap(~ multiexonic)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(ref_genes_info,aes(x=simpl_gene_biotype,fill=best_transfrag_cc )) +
  geom_bar(position="fill") +
  theme_classic() +ylab("fraction of genes") + ggtitle("Fraction of reference genes overlapped by transfrags") + facet_wrap(~ multiexonic)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# further divide classcodes in good overlap and other/antisense overlap
good_cc=order_cc[1:5]
other_cc=order_cc[6:length(order_cc)]
ref_genes_info$type_overlap=ifelse(is.na(ref_genes_info$best_transfrag_cc),NA,
                                   ifelse(ref_genes_info$best_transfrag_cc%in%good_cc,"good_overlap","other_overlap"))

ggplot(ref_genes_info,aes(x=simpl_gene_biotype,fill=type_overlap )) +
  geom_bar(position="fill") +
  theme_classic() +ylab("fraction of genes") + ggtitle("Fraction of reference genes overlapped by transfrags") + facet_wrap(~ multiexonic)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(ref_genes_info,aes(x=best_transfrag_cc,fill=simpl_gene_biotype )) +
  geom_bar() +
  theme_classic() +ylab("# genes") + ggtitle("Number of ref gene biotypes in each classcode") + facet_wrap(~ multiexonic)

ggplot(ref_genes_info,aes(x=best_transfrag_cc,fill=transfrag_exonic_type )) +
  geom_bar() +
  theme_classic() +ylab("# genes") + ggtitle("Exonic type of transfrags overlapping monoexonic and multiexonic ref genes") + facet_wrap(~ multiexonic)

ggplot(ref_genes_info,aes(x=simpl_gene_biotype,fill=transfrag_exonic_type )) +
  geom_bar() +
  theme_classic() +ylab("# genes") + ggtitle("Reference genes overlapped by tranfrags") + facet_wrap(~ multiexonic)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(ref_genes_info,aes(x=simpl_gene_biotype,fill=transfrag_exonic_type )) +
  geom_bar(position = "fill") +
  theme_classic() +ylab("fraction of genes") + ggtitle("Fraction of reference genes overlapped by tranfrags") + facet_wrap(~ multiexonic)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
