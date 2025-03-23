# genes that do not have an assigned gene name in refseq
# but whose transcripts do have a gencode gene name (through
# gffcompare) get, in principle, the gencode_gene_name
# posibilities: its OK, their transcripts belong to a unique gffcompare id

gffcmp.gene_id_woRefseqID=merged_ref_annot$V2[is.na(merged_ref_annot$gene_name)]
length(unique(gencode_gene_names_noRefseq))
length(unique())

# get number of XLOC per gene_name
gencode_XLOC=merged_ref_annot%>%filter(gencode_gene!="-")%>%group_by(gencode_gene_name)%>%
  reframe(XLOC=unique(V2))
refseq_XLOC=merged_ref_annot%>%filter(refseq_gene!="-")%>%group_by(refseq_gene)%>%
  reframe(XLOC=unique(V2))
## # in "summarise at gene level" NOTE: in principle I should check if those transcripts without gene id in gencode
# are isoforms of genes with assigned id in gencode and in that case,
# if the gene names do match between refseq and gencode. For now I haven't done this
##
##
# setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/")
# source("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/scripts/StemLinc_code/functions_to_parse_gffCompare_tracking_files.R")


# NOTE: in principle I should check if those transcripts without gene id in gencode
# are isoforms of genes with assigned id in gencode and in that case,
# if the gene names do match between refseq and gencode. For now I haven't done this



# get genomic coordinates
refGTF$old_gene_id=refGTF$gene_id
refGTF$gene_name=annotated_track$gene_name[match(refGTF$transcript_id,
                                                 annotated_track$V1)]
refGTF$gene_id=annotated_track$gene_name[match(refGTF$transcript_id,
                                               annotated_track$V1)]
refGTF_genomic_regions=trastools::get_genomic_range_by_gene(refGTF)

annotated_track_gene_level <- left_join(annotated_track_gene_level,
                                        refGTF_genomic_regions,by=c(gene_name="gene_id"))

table(annotated_track_gene_level$Nchr)
# when summarising at gene level there are a few too long genes
# the ones I checked were because they were different loci annotated
# with the same gene_name

annotated_track_uniq_idname=annotated_track
annotated_track_uniq_idname$id_name=paste(annotated_track$gencode_gene,
                                          annotated_track$gene_name,sep = ";")
annotated_track_uniq_idname=annotated_track_uniq_idname%>%filter(!duplicated(id_name))
dup_gene_names=annotated_track_uniq_idname$gene_name[duplicated(annotated_track_uniq_idname$gene_name)]

gff_compare_genes_per_gene_name=refGTF%>%group_by(gene_name)%>%summarise(N_XLOC=length(unique(old_gene_id)))
# I could just make a reduce using the entire transcript coordinates for each gene
# and check if there are genes composed of non-overlapping transcripts
# from different annotations...

## to be continued ...




# for now, save what I've done so far
# write temp out files ----
export(refGTF,"~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/mouse/annotation/merged_refs_annotation/merged_refs.combined.with_gene_name.gtf",format = "gtf")
write.table(refGTF_genomic_regions,"~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/mouse/annotation/merged_refs_annotation/refGTF_genomic_regions.txt",sep="\t",quote = F,row.names = F)
write.table(annotated_track,"~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/mouse/annotation/merged_refs_annotation/annotated_tracking_file.txt",quote = F,row.names = F,sep="\t")
write.table(annotated_track_gene_level,"~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/mouse/annotation/merged_refs_annotation/annotated_tracking_file_gene_level.txt",quote = F,row.names = F,sep="\t")

EPT_refGTF=trastools::exons_per_transcript(refGTF)


merged_ref_annot_gene_level=merged_ref_annot %>% group_by(gene_name)%>%
  summarise(Nchr=length(unique(chr)),
            chrs=paste0(sort(unique(chr)),collapse = ","),
            N_gencode_biot=length(unique(gencode_biotype[!is.na(gencode_biotype)])),
            gencode_biots=paste0(sort(unique(gencode_biotype)),collapse = ","),
            N_refseq_biot=length(unique(refseq_biotype[!is.na(refseq_biotype)])),
            refseq_biots=paste0(sort(unique(refseq_biotype)),collapse = ","))

# special case B
#                                    B.1) Are there any gene_id within the same ref that has two different assigned names (A and B)?
#                                         if a matching transcript that in ref 1 has name A and in the other ref has B, then take B (maybe
#                                         A is misannotated, or majority vote)
#                                    B.2) if also in the other reference, do their transcripts match transcript (or gene name) from the other reference?
#                                    then could be truly different overlapping transcriptional units. For this, transcript level quantification might be useful
#                                    B.3) do they also have any of the same genes in the other reference?, majority vote
#                                    B.1) do they also have overlapping genes in the other reference but do not match transcripts?
#


setwd("/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/")

source("/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/scripts/StemLinc_code/compute_gene_non_redundant_exons_and_exonic_length.functions.R")
source("/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/scripts/StemLinc_code/functions_to_parse_gffCompare_tracking_files.R")
source("/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/scripts/StemLinc_code/get_genomic_regions_function.R")
source("/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/scripts/StemLinc_code/get_GTF_info_function.R")
# gtf_path="/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare/GENCODE_lncRNAs/StemLinc_vs_GENCODElncRNAs_blind.combined.gtf"
# gtf_path="/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/analyses/transcriptome_characterization_and_comparison/input_files/gtfs/LSK_StemLinc.combined.gtf"
# tracking_path="/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare_blind/LSK_StemLinc.blind.tracking"
# tracking_path="/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/analyses/transcriptome_characterization_and_comparison/input_files/tracking/LSK_StemLinc.tracking"
#
#
# ### gtf ----
# gtf=readGFF(gtf_path)
# ### tracking file ----
# tracking=read.table(tracking_path)

### ref annot ----
# refGTF$gene_id=merged_ref_annot$gene_name[match(refGTF$transcript_id,
#                                                 merged_ref_annot$V1)]

refseq=readGFF("/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/references/mouse/annotation/RefSeq/GCF_000001635.27_GRCm39_genomic.nocmm.chr.gff")


# I do have to check if the refseq gene names that do not match any gencode name
# are present in other genes
table(non_dup_XLOC$refseq[non_dup_XLOC$type_gene=="no_match"]%in%merged_ref_annot$gene_name)
table(non_dup_XLOC$refseq[non_dup_XLOC$type_gene=="no_match"]%in%merged_ref_annot$gencode_gene_name)
# none of them is, then I can just take the gencode gene_name for those genes
# in the same way, check if the refseq only genes match gencode genes
table(non_dup_XLOC$refseq[non_dup_XLOC$type_gene=="refseq"]%in%merged_ref_annot$gencode_gene_name)
# only 13


# rule: if for a given shared transcript gencode and refseq have names,
# take gencode one
# if for a refseq transcript in a same XLOC
# there is a name not present in the gencode names
# associated with that XLOC, check,
# if it has a different biotype, keep as a
# separate gene, if it shares biotype with one of the XLOC biotypes,
# assign the name


# # only refseq:
#
# only_refseq_XLOC=still_complex_XLOC_refseq_only_nomatch%>%filter(gencode_biotype=="")
#
# still_complex_XLOC_refseq_only_nomatch=still_complex_XLOC_refseq_only_nomatch%>%filter(gencode_biotype!="")

write.table(merged_ref_annot,"/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/references/mouse/annotation/merged_refs_annotation/annotated_tracking_file.updated_gene_names.txt",row.names = F,quote = F,sep = "\t")
merged_ref_annot1=read.table("~/Documentos/references/annotated_tracking_file.updated_gene_names.txt",sep = "\t",header = T)


merged_transcripts_per_gene <- get_per_gene_merged_transcripts(refGTF_new_gene_id)
length(unique(merged_transcripts_per_gene$group_name))
length(unique(refGTF_new_gene_id$gene_id))

in_bed <- merged_ref_annot%>%dplyr::select(seqid,start,end,V1,gene_name,strand)

in_bed_test=in_bed[order(in_bed[,1],in_bed[,2]),]
in_bed_test=in_bed_test[1:50,]
