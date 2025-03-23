options(scipen = 999)
## Load packages---------------------------
require(tidyverse)
require(data.table)
require(rtracklayer)
require(trastools)
library(bedtoolsr)

# New combined transcriptome filtering.
# Individual assemblies were generated with StringTie v.2.1.7
# e.g: bam="LSK_1.bam", guided by GENCODE VM31
# stringtie -p $PPN -v -G $ref_gtf -M 0.3 --rf -o $bam.strtie.guided.gtf $bam
# The individual assemblies were combined and compared with gffCompare and against the
# latest version of GENCODE:
# refgtf="/ijc/USERS/llorenzi/references/gencode.vM36.annotation.gtf"
# cat gtf_paths:
# /mnt/beegfs/llorenzi/jobs/LSK_RNA_seq/analyses/LSK_1/LSK_1Aligned.sortedByCoord.out.MarkedDups.minq1.primarychrs.bam.strtie.guided.gtf
# /mnt/beegfs/llorenzi/jobs/LSK_RNA_seq/analyses/LSK_2/LSK_2Aligned.sortedByCoord.out.MarkedDups.minq1.primarychrs.bam.strtie.guided.gtf
# /mnt/beegfs/llorenzi/jobs/LSK_RNA_seq/analyses/LSK_3/LSK_3Aligned.sortedByCoord.out.MarkedDups.minq1.primarychrs.bam.strtie.guided.gtf
# /mnt/beegfs/llorenzi/jobs/T_cell_macrophage_RNA_seq/analyses/Macro_R1/Macro_R1Aligned.sortedByCoord.out.MarkedDups.minq1.primarychrs.bam.strtie.guided.gtf
# /mnt/beegfs/llorenzi/jobs/T_cell_macrophage_RNA_seq/analyses/Macro_R2/Macro_R2Aligned.sortedByCoord.out.MarkedDups.minq1.primarychrs.bam.strtie.guided.gtf
# /mnt/beegfs/llorenzi/jobs/T_cell_macrophage_RNA_seq/analyses/Macro_R3/Macro_R3Aligned.sortedByCoord.out.MarkedDups.minq1.primarychrs.bam.strtie.guided.gtf
# /mnt/beegfs/llorenzi/jobs/T_cell_macrophage_RNA_seq/analyses/T_cells_R1/T_cells_R1Aligned.sortedByCoord.out.MarkedDups.minq1.primarychrs.bam.strtie.guided.gtf
# /mnt/beegfs/llorenzi/jobs/T_cell_macrophage_RNA_seq/analyses/T_cells_R2/T_cells_R2Aligned.sortedByCoord.out.MarkedDups.minq1.primarychrs.bam.strtie.guided.gtf
# /mnt/beegfs/llorenzi/jobs/T_cell_macrophage_RNA_seq/analyses/T_cells_R3/T_cells_R3Aligned.sortedByCoord.out.MarkedDups.minq1.primarychrs.bam.strtie.guided.gtf

# gffcompare -i gtf_paths.txt -r $refgtf -s $smaskedgenome -V -j $outdir/StemLinc_replicates.vs.VM36.new_junction.tsv -o $outdir/StemLinc_replicates.vs.VM36


# Load data ----
tracking_path <- "~/Documentos/all_samples_analyses/gffcompare_replicates.vs.VM36/StemLinc_replicates.vs.VM36.tracking"
gtf_path <- "~/Documentos/all_samples_analyses/gffcompare_replicates.vs.VM36/StemLinc_replicates.vs.VM36.combined.gtf"
ref_gtf_path <- "~/Documentos/references/gencode.vM36.annotation.gtf"
tracking <- read.table(tracking_path)
gtf <- readGFF(gtf_path)

# Extract transcript features ----
N_exons <- as.integer(trastools::get_n_exons_per_transcript(tracking))
N_libraries <- trastools::get_n_samples_per_transcript_from_tracking(tracking)
TPMs <- trastools::get_expression_values(tracking)
table(N_exons,N_libraries)

ordered_classcodes <- c("=", "j", "k", "c", "m", "n", "e", "o", "p", "s", "x", "i", "y", "u", "r", ".")

tracking$N_exons=N_exons
tracking$N_libraries=N_libraries
tracking$mean_tpm=TPMs$mean_tpm
tracking <- tracking %>% mutate(overlapRef=V4%in%overlapping_class_codes)
table(tracking$overlapRef)
rm(N_exons)
rm(N_libraries)
gc()
ref_gtf <- readGFF(ref_gtf_path)
tracking <- separate(tracking,col = 3,into = c("ref_gene_id","ref_transcript_id"),sep = "\\|")
tracking$ref_transcript_id[is.na(tracking$ref_transcript_id)] <- "-"
tracking$ref_biotype <- ref_gtf$gene_type[match(tracking$ref_transcript_id,ref_gtf$transcript_id)]

# simplify ref biotypes ----
ref_biotypes=unique(tracking$ref_biotype)
#simplify biotypes
simpl_biotypes <- list("protein_coding" = c("protein_coding",
                                            grep("_gene|_V_",
                                                 grep("pseudogene",
                                                      ref_biotypes,
                                                      invert = T,value = T),value = T)),
                       "lncRNA" = "lncRNA",
                       "pseudogene" = grep("pseudogene",
                                           ref_biotypes,
                                           value = T),
                       "TEC" = grep("TEC",
                                    ref_biotypes,
                                    value = T),
                       "miRNA" =grep("miRNA",
                                     ref_biotypes,
                                     value = T),
                       "sno-snRNA" = c("snoRNA"  ,
                                       "snRNA",
                                       "ncRNA;snoRNA"),
                       "rRNA" = grep("rRNA",
                                     ref_biotypes,
                                     value = T),
                       "tRNA" = grep("tRNA",
                                     ref_biotypes,
                                     value = T))
#"other" = )

table(ref_biotypes%in%unlist(simpl_biotypes))

ref_biotypes[!ref_biotypes%in%unlist(simpl_biotypes)]

simpl_biotypes <- c(simpl_biotypes,
                    list(other = unique(ref_biotypes)[!unique(ref_biotypes)%in%
                                                        unlist(simpl_biotypes)&!is.na(unique(ref_biotypes))]))
simpl_biotypes_table <- data.frame(orig.biot=unlist(simpl_biotypes),
                                   final.biot=rep(names(simpl_biotypes),sapply(simpl_biotypes,length)))

tracking$ref_biotype=simpl_biotypes_table$final.biot[match(tracking$ref_biotype,
                                                                               simpl_biotypes_table$orig.biot)]


table(tracking$ref_biotype)
# add transcript coordinates ----
tracking <- left_join(tracking, gtf%>%filter(type=="transcript") %>%
                        dplyr::select(V1=transcript_id,seqid,start,end,strand))

# add ref strand ----
tracking$ref_strand <- ref_gtf$strand[match(tracking$ref_transcript_id,ref_gtf$transcript_id)]


# Check numbers ----
tracking <- tracking %>% mutate(gene_id=ifelse(overlapRef,ref_gene_id,V2))

# Assembled annotated lncRNAs ----
Ref_genes <- ref_gtf %>% filter(type=="exon") %>%
  mutate(transcript_support_level=ifelse(is.na(transcript_support_level),"NA",transcript_support_level)) %>%
  group_by(transcript_id) %>% mutate(N_exons=n()) %>%
group_by(gene_type,gene_id,gene_name) %>% arrange(transcript_support_level) %>%
  summarise(min_level=min(level),
            min_transcript_support_level=transcript_support_level[1],
            max_N_exons=max(N_exons))

Ref_genes_summ <- Ref_genes %>% group_by(gene_type)%>%
  summarise(N_genes=n(),
            monoexonic=sum(max_N_exons==1),
            multiexonic=sum(max_N_exons>1),
            lev1=sum(min_level==1),
            lev2=sum(min_level==2),
            lev3=sum(min_level==3),
            TSL_1=sum(min_transcript_support_level==1),
            TSL_2=sum(min_transcript_support_level==2),
            TSL_3=sum(min_transcript_support_level==3),
            TSL_4=sum(min_transcript_support_level==4),
            TSL_5=sum(min_transcript_support_level==5),
            TSL_NA=sum(min_transcript_support_level=="NA"))

# write.table(Ref_genes,"data/references/gencode_vM36/Ref_genes_per_biotype_evidence_level.tsv",
#             quote = F,row.names = F,sep = "\t")
# write.table(Ref_genes_summ,"data/references/gencode_vM36/Ref_genes_per_biotype_evidence_level_summary.tsv",
#             quote = F,row.names = F,sep = "\t")

# write.table(Ref_genes$gene_name[Ref_genes$gene_type=="lncRNA"&Ref_genes$max_N_exons==1],
#             "data/references/gencode_vM36/GENCODE_vM36_monoexonic_lncRNAs.txt",quote = F,col.names = F,row.names = F)

tracking <- tracking %>% group_by(V1) %>%
  mutate(transfrag_biotype=ifelse(overlapRef,ref_biotype,paste(c(V4,na.omit(ref_biotype)),collapse = "-")))

table(tracking$transfrag_biotype)
biots=c("protein_coding","lncRNA","pseudogene","TEC")

tracking_biots <- tracking %>% filter(ifelse(!overlapRef,T,ifelse(ref_biotype%in%biots,T,F)))

most_freq_tranfrag_biots <- names(sort(table(tracking_biots$transfrag_biotype),decreasing = T))[1:10]

tracking_biots <- tracking_biots %>%
  mutate(transfrag_biotype=ifelse(transfrag_biotype%in%most_freq_tranfrag_biots,
                                  transfrag_biotype,
                                  "other"))
rm(ref_gtf)
gc()
nejm_pal=ggsci::pal_nejm()(8)

g=ggplot(tracking_biots,
       aes(fct_infreq(transfrag_biotype),
           fill = as.factor(N_libraries))) + geom_bar() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_fill_manual(values=colorRampPalette(colors = nejm_pal)(9)) +
  ylab("# transfrags")

svg("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/pre_filter_combined_transcriptome/combined_transcriptome_transfrag_distribution_N_libraries.svg")
print(g)
dev.off()

g=ggplot(tracking_biots %>% filter(N_libraries>1),
         aes(fct_infreq(transfrag_biotype),
             fill = as.factor(N_libraries))) + geom_bar() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_fill_manual(values=colorRampPalette(colors = nejm_pal)(9)[2:9]) +
  ylab("# transfrags")

g

svg("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/pre_filter_combined_transcriptome/combined_transcriptome_transfrag_distribution_N_libraries_atleast2.svg")
print(g)
dev.off()

# summarise at gene level ----

tracking_biots_2samps <- tracking_biots %>% filter(N_libraries>1)

class_order <- names(sort(table(tracking_biots_2samps$transfrag_biotype),
                          decreasing = T))

# N_libsGL <- get_n_samples_per_gene_from_tracking(tracking_biots_2samps %>%
#                                                    mutate(V2=gene_id),
#                                                  cols = 6:14)

cols = 6:14
sample_occurrence_per_gene = lapply(cols, function(sa) {
  stats::aggregate(tracking_biots_2samps[, sa],
                   by = list(gene_id = tracking_biots_2samps$gene_id),
                   function(x) any(x != "-"))
})

gene_id = sample_occurrence_per_gene[[1]]$gene_id
sample_occurrence_per_gene <- as.data.frame(sapply(sample_occurrence_per_gene,
                                                   function(x) x[,2]))
sample_occurrence_per_gene$gene_id = gene_id
Nsamps = rowSums(sample_occurrence_per_gene[, 1:length(cols)])
names(Nsamps) = sample_occurrence_per_gene$gene_id

tracking_biots_2samps_GL <- tracking_biots_2samps %>%
  group_by(gene_id) %>% arrange(match(transfrag_biotype,class_order)) %>%
  summarise(exonic_type=ifelse(any(N_exons>1),"multiexonic","monoexonic"),
          gene_biotype=transfrag_biotype[1])
tracking_biots_2samps_GL$N_libraries <-
  Nsamps[match(tracking_biots_2samps_GL$gene_id,names(Nsamps))]

g=ggplot(tracking_biots_2samps_GL,aes(x = fct_infreq(gene_biotype), fill = as.factor(N_libraries))) +
  geom_bar() +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_fill_manual(values=colorRampPalette(colors = nejm_pal)(9)[2:9]) +
  ylab("# Loci")
svg("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/pre_filter_combined_transcriptome/combined_transcriptome_transfrag_distribution_N_libraries_atleast2.GENE_LEVEL.svg")
print(g)
dev.off()

summarise_at_gene_level <- function(ftracking,cols=6:14){

  sample_occurrence_per_gene = lapply(cols, function(sa) {
    stats::aggregate(ftracking[, sa], by = list(gene_id = ftracking$gene_id),
                     function(x) any(x != "-"))
  })

  gene_id = sample_occurrence_per_gene[[1]]$gene_id
  sample_occurrence_per_gene <- as.data.frame(sapply(sample_occurrence_per_gene,
                                                     function(x) x[,2]))
  sample_occurrence_per_gene$gene_id = gene_id
  Nsamps = rowSums(sample_occurrence_per_gene[, 1:length(cols)])
  names(Nsamps) = sample_occurrence_per_gene$gene_id

  class_order <- names(sort(table(ftracking$transfrag_biotype),decreasing = T))

  tracking_GL <- ftracking %>%
    group_by(gene_id) %>% arrange(match(transfrag_biotype,class_order)) %>%
    summarise(exonic_type=ifelse(any(N_exons>1),"multiexonic","monoexonic"),
              gene_biotype=transfrag_biotype[1])

  best_cc <- ftracking %>% group_by(gene_id) %>%
    arrange(match(V4,ordered_classcodes)) %>% summarise(best_cc=V4[1])

  tracking_GL$N_libraries <-
    Nsamps[match(tracking_GL$gene_id,names(Nsamps))]

   tracking_GL <- left_join(tracking_GL, best_cc)

  return(tracking_GL)

}

# remove spurious classcodes ----
spurious_ccs <- c("e","p","o","s",".")
filtered_tracking <- tracking %>% filter(!V4%in%spurious_ccs)

# remove strand "*" ----
filtered_tracking <- filtered_tracking%>%filter(strand%in%c("+","-"))
length(unique(filtered_tracking$V2[!filtered_tracking$overlapRef]))

table(filtered_tracking$N_exons>1&filtered_tracking$N_libraries>1)

View(filtered_tracking%>%filter(!overlapRef,N_exons>1))
filtered_tracking_GL <- summarise_at_gene_level(filtered_tracking)
g=ggplot(filtered_tracking_GL,aes(x = fct_infreq(gene_biotype), fill = as.factor(N_libraries))) +
  geom_bar() +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_fill_manual(values=colorRampPalette(colors = nejm_pal)(9)) +
  ylab("# Loci")
g
svg("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/pre_filter_combined_transcriptome/combined_transcriptome_transfrag_distribution_N_libraries_atleast2.GENE_LEVEL.filter_spurious_ccs.svg")
print(g)
dev.off()

# remove transcripts with mean TPM≤0.2 ----
filtered_tracking <- filtered_tracking %>% filter(mean_tpm>=0.2)
table(filtered_tracking$N_exons,filtered_tracking$N_libraries)

length(unique(filtered_tracking$V2[!filtered_tracking$overlapRef]))

tr_support2=filtered_tracking%>%filter(N_libraries>1)

filtered_tracking_GL <- summarise_at_gene_level(filtered_tracking)
g=ggplot(filtered_tracking_GL,aes(x = fct_infreq(gene_biotype), fill = as.factor(N_libraries))) +
  geom_bar() +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_fill_manual(values=colorRampPalette(colors = nejm_pal)(9)) +
  ylab("# Loci")
g
svg("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/pre_filter_combined_transcriptome/combined_transcriptome_transfrag_distribution_N_libraries_atleast2.GENE_LEVEL.filter_spurious_ccs.TPM0.2.svg")
print(g)
dev.off()
# remove intronic to PCGs ----
intronic_sense <- filtered_tracking %>%
  filter(V4=="i"&ref_biotype=="protein_coding"&ref_strand==strand)
filtered_tracking <- filtered_tracking %>% filter(!V1%in%intronic_sense$V1)

filtered_tracking_GL <- summarise_at_gene_level(filtered_tracking)
g=ggplot(filtered_tracking_GL,aes(x = fct_infreq(gene_biotype), fill = as.factor(N_libraries))) +
  geom_bar() +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_fill_manual(values=colorRampPalette(colors = nejm_pal)(9)) +
  ylab("# Loci")
g
svg("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/pre_filter_combined_transcriptome/combined_transcriptome_transfrag_distribution_N_libraries_atleast2.GENE_LEVEL.filter_spurious_ccs.TPM0.2.notintronicPCG.svg")
print(g)
dev.off()
# remove all genes that are only found in one library ----
genes2keep <- filtered_tracking_GL$gene_id[filtered_tracking_GL$N_libraries>1]
filtered_tracking <- filtered_tracking[filtered_tracking$gene_id%in%genes2keep,]

filtered_tracking_GL <- summarise_at_gene_level(filtered_tracking)
g=ggplot(filtered_tracking_GL,aes(x = fct_infreq(gene_biotype), fill = as.factor(N_libraries))) +
  geom_bar() +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_fill_manual(values=colorRampPalette(colors = nejm_pal)(9)) +
  ylab("# Loci")
g
svg("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/pre_filter_combined_transcriptome/combined_transcriptome_transfrag_distribution_N_libraries_atleast2.GENE_LEVEL.filter_spurious_ccs.TPM0.2.notintronicPCG.2libsGL.svg")
print(g)
dev.off()
length(unique(filtered_tracking$V2[!filtered_tracking$overlapRef]))
length(unique(filtered_tracking$V2[filtered_tracking$overlapRef]))

filtered_tracking_GL %>% group_by(gene_biotype) %>% summarise(count=n())

# check how many annotated lncRNAs are assembled ----
lncRNAs_assembled <- filtered_tracking_GL$gene_id[filtered_tracking_GL$gene_biotype=="lncRNA"]
Ref_genes_lncRNAs_assembled <- Ref_genes%>%filter(gene_id%in%lncRNAs_assembled)
table(Ref_genes_lncRNAs_assembled$max_N_exons)
table(Ref_genes_lncRNAs_assembled$max_N_exons>1)
table(Ref_genes_lncRNAs_assembled$min_level)
table(Ref_genes_lncRNAs_assembled$min_transcript_support_level)
table(Ref_genes_lncRNAs_assembled$min_transcript_support_level,
      Ref_genes_lncRNAs_assembled$max_N_exons>1)

filtered_tracking_GL <- left_join(filtered_tracking_GL,Ref_genes)
write.table(filtered_tracking_GL,"outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation/filtered_assembled_loci_annotated_and_novel.tsv",
            quote = F,row.names = F,sep = "\t")
# keep only transcripts non-overlapping reference ----
filtered_tracking <- filtered_tracking %>% filter(!overlapRef)
filtered_tracking_GL <- summarise_at_gene_level(filtered_tracking)
g=ggplot(filtered_tracking_GL,aes(x = fct_infreq(gene_biotype), fill = as.factor(N_libraries))) +
  geom_bar() +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_fill_manual(values=colorRampPalette(colors = nejm_pal)(9)) +
  ylab("# Loci")
g
svg("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/pre_filter_combined_transcriptome/combined_transcriptome_transfrag_distribution_N_libraries_atleast2.GENE_LEVEL.filter_spurious_ccs.TPM0.2.notintronicPCG.2libsGL.noOlapRef.svg")
print(g)
dev.off()
# assign library occurrence ----
samples=c(paste0("LSK_",1:3),
          paste0("macro_",1:3),
          paste0("T.cell_",1:3))
filtered_tracking_libraries <- apply(filtered_tracking[6:14],1, function(x){paste(samples[x!="-"],collapse = ",")})
filtered_tracking$libraries <- filtered_tracking_libraries

filtered_tracking_bed <- filtered_tracking %>%
  mutate(start=start-1,exonic_type=ifelse(N_exons==1,"mono","multi")) %>% select(seqid,start,end,V2,libraries,strand,N_libraries,exonic_type,V1) %>%
  arrange(seqid,start)

# process mono and multiexonic clusters separately ----

# monoexonic clusters are easy to collapse - for multi, check later
# merge monoexons and filter out lowly supported ones
# assign unique transcript_id and gene name later (only one per gene)


# check if monoexons overlap exons from multiexons, if so annotate them with the same
# gene id, later we can further filter those.

# check if, after filtering, genes are splitted in two

# Process new filtered monoexons ----
monoexons_bed <- filtered_tracking_bed %>% filter(exonic_type=="mono")

# merge monoexons ----
merged_monoexons <- bt.merge(i = monoexons_bed,
                                             s = T,c = "4,5,6,7,8",
                                             o="distinct,distinct,distinct,max,distinct")

# count the number of distinct libraries after merging:
Nlibs <- sapply(merged_monoexons$V5,
                function(x){
                  splits=strsplit(x,split = ",")[[1]]
                  return(length(unique(splits)))
                })

merged_monoexons$N_libraries <- Nlibs

# keep only monoexons supported by at least 3 libraries ----
filtered_monoexons <- merged_monoexons %>% filter(N_libraries>2)

# Create unique gene ids
table(duplicated(filtered_monoexons$V4))
# 72

filtered_monoexons$new_gene_id <- ave(filtered_monoexons$V4, filtered_monoexons$V4, FUN = function(x) {
  if (length(x) > 1) {
    paste0(x, ".", seq_along(x))
  } else {
    x
  }
})
length(unique(filtered_monoexons$new_gene_id))

# process multiexons ----
multiexons_bed <- filtered_tracking_bed %>% filter(exonic_type=="multi")
# In this case I want to merge by gene id only, as bedtools merge would merge
# together intronic genes with the host gene

# Convert to GenomicRanges object
gr <- GRanges(seqnames = multiexons_bed$seqid,
              ranges = IRanges(start = multiexons_bed$start, end = multiexons_bed$end),
              gene_id = multiexons_bed$V2, strand = multiexons_bed$strand,
              transcript_id=multiexons_bed$V1)

# gene_ids <- unique(gr$gene_id)
# gene_overlaps <- sapply(gene_ids, function(gid) {
#   merge <- reduce(gr[gr$gene_id == gid])
#   length(merge)  # Number of merged regions
# })
#
# # Find genes with multiple loci
# split_genes <- names(gene_overlaps[gene_overlaps > 1])
# print(split_genes)  # List of genes that may need to be split
#
# length(split_genes)

# **Reduce separately per gene_id**
merged_gr_list <- gr %>%
  split(., .$gene_id) %>%  # Split by gene_id
  lapply(reduce) %>%       # Reduce within each gene group
  GRangesList() %>%
  unlist()

# merged_genes <- gene_overlaps <- sapply(gene_ids, function(gid) {
#   merge <- reduce(gr[gr$gene_id == gid])
#   return(merge)  # Number of merged regions
# })


# Preserve gene_id in the merged clusters
merged_gr_list$gene_id <- names(merged_gr_list)

# **Step 2: Assign Transcripts to Merged Clusters**
overlaps <- findOverlaps(gr, merged_gr_list, type = "any",ignore.strand=F)

# Create a dataframe linking transcripts to merged clusters
assignment_df <- data.frame(
  transcript_id = gr$transcript_id[queryHits(overlaps)],
  cluster_start = start(merged_gr_list)[subjectHits(overlaps)],
  cluster_end = end(merged_gr_list)[subjectHits(overlaps)],
  gene_cluster_id = merged_gr_list$gene_id[subjectHits(overlaps)]
)

# assign the original gene id to each transcript
assignment_df$gene_id <- multiexons_bed$V2[match(assignment_df$transcript_id,
                                                 multiexons_bed$V1)]

# keep only matches within the same gene id
assignment_df <- assignment_df %>%filter(gene_id==gene_cluster_id)
assignment_df$cluster_id <- paste(assignment_df$cluster_start,
                                  assignment_df$cluster_end,
                                  sep = "_")


#
gene_cluster <- assignment_df[,5:6]
# keep unique gene-cluster pairs:
gene_cluster <- gene_cluster[!duplicated(gene_cluster),]
dups=unique(gene_cluster$gene_id[duplicated(gene_cluster$gene_id)])
length(dups)
#length(split_genes)

#table(split_genes%in%dups)
# assign new gene ids to genes with more than 1 cluster
gene_cluster$new_gene_id <- ave(gene_cluster$gene_id, gene_cluster$gene_id,
                                 FUN = function(x) {
  if (length(x) > 1) {
    paste0(x, ".", seq_along(x))
  } else {
    x
  }
})

nrow(gene_cluster)==length(unique(gene_cluster$cluster_id))
assignment_df$new_gene_id <- gene_cluster$new_gene_id[match(assignment_df$cluster_id,
                                                            gene_cluster$cluster_id)]

# assign the newly defined gene_ids to the multiexons bed
multiexons_bed$new_gene_id <- assignment_df$new_gene_id[match(multiexons_bed$V1,
                                                              assignment_df$transcript_id)]

# summarise at gene level using the new gene cluster definition:
merged_multiexons <- multiexons_bed %>% group_by(new_gene_id) %>%
  summarise(seqid=seqid[1],
            start=min(start),
            end=max(end),
            libraries=paste(unique(libraries),collapse = ","),
            strand=strand[1],
            max_libraries=max(N_libraries),
            exonic_type=exonic_type[1])

Nlibs <- sapply(merged_multiexons$libraries,
                function(x){
                  splits=strsplit(x,split = ",")[[1]]
                  return(length(unique(splits)))
                })

merged_multiexons$N_libraries <- Nlibs

filtered_multiexons <- merged_multiexons %>%filter(N_libraries>1)
filtered_multiexons_bed <- multiexons_bed%>%filter(new_gene_id%in%filtered_multiexons$new_gene_id)

# merged_multiexons2 <- bt.merge(i = multiexons_bed,
#                              s = T,c = "5,6,7,9,10",
#                              o="distinct,distinct,max,distinct,distinct")
#N_libraries max_libraries

# check if multiexon blocks that share the same gene id are splitted
# If so, use findoverlaps to map the isoforms to new clusters,
# apply filtering



# for those, re-check the sample support
# change name if necessary
# multiexonic clusters should be supported by splicing from at least 2 libraries,
# while monoexonic clusters shoud be supported by 3 libraries or more,
# after this filtering, if a monoexon isoform overlaps an exon of a multiexons
# it gets that gene_id

filtered_multiexon_transcripts <- filtered_multiexons_bed$V1

multiexons_exon_bed <- gtf%>%
  filter(type=="exon",transcript_id%in%filtered_multiexon_transcripts) %>%
  mutate(start=start-1) %>%
  select(seqid,start,end,transcript_id,gene_id,strand)

multiexons_exon_bed$gene_id <- multiexons_bed$new_gene_id[match(multiexons_exon_bed$transcript_id,
                                                                    multiexons_bed$V1 )]
filtered_monoexons_bed <- filtered_monoexons %>%dplyr::select(V1,V2,V3,new_gene_id,N_libraries,V6)

# perform overlap between monoexons and exons from multiexons
multi_intersect_mono <- bt.intersect(a=filtered_monoexons_bed,multiexons_exon_bed,wo = T,s = T)

#summarise by new gene_id
summarised_multi_intersect_mono <- multi_intersect_mono %>% group_by(V4) %>%
  summarise(new_gene_id=paste0(unique(V11),collapse = ","),
            N_genes=length(unique(V11)))

# remove monoexons overlapping more than one multiexon gene if any:
monoexons2remove <- summarised_multi_intersect_mono$V4[summarised_multi_intersect_mono$N_genes>1]
filtered_monoexons <- filtered_monoexons %>%filter(!new_gene_id%in%monoexons2remove)

# if overlapping multiexon, assign gene id of multiexon:
filtered_monoexons$new_gene_id[filtered_monoexons$new_gene_id%in%summarised_multi_intersect_mono$V4] <-
  summarised_multi_intersect_mono$new_gene_id[match(filtered_monoexons$new_gene_id[filtered_monoexons$new_gene_id%in%summarised_multi_intersect_mono$V4],
                                                    summarised_multi_intersect_mono$V4)]

filtered_monoexons$width <- filtered_monoexons$V3 - filtered_monoexons$V2

# remove too short monoexons ----
filtered_monoexons <- filtered_monoexons %>% filter(width>=200)

# generate GTF ----
gtf_multi <- gtf%>% filter(transcript_id%in%filtered_multiexon_transcripts)
gtf_multi$gene_id <- filtered_multiexons_bed$new_gene_id[match(gtf_multi$transcript_id,
                                                           filtered_multiexons_bed$V1)]

gtf_multi$source="StemLinc"

gtf_multi <- gtf_multi[,1:10]
# modify gene_ids with ","
filtered_monoexons$new_gene_id[grep(",",filtered_monoexons$new_gene_id)] <-
  sapply(strsplit(filtered_monoexons$new_gene_id[grep(",",filtered_monoexons$new_gene_id)],split = ","),
         function(x)x[[1]][1])
filtered_monoexons$new_gene_id

dup_genes <- filtered_monoexons$new_gene_id[duplicated(filtered_monoexons$new_gene_id)]

filtered_monoexons$tr_id <- ave(filtered_monoexons$new_gene_id, filtered_monoexons$new_gene_id, FUN = function(x) {
  if (length(x) > 1) {
    paste0(x, "_", seq_along(x))
  } else {
    x
  }
})

filtered_monoexons$tr_id <- gsub("XLOC","MONO",filtered_monoexons$tr_id)
gtf_mono <- filtered_monoexons %>%
  mutate(source="StemLinc",type="transcript",score=NA,phase=NA) %>%
  dplyr::select(V1,source,type,V2,V3,score,V6,phase,tr_id,new_gene_id)

gtf_mono_exons <- gtf_mono %>%mutate(type="exon")

gtf_mono <- rbind(gtf_mono,gtf_mono_exons)

colnames(gtf_mono) <- colnames(gtf_multi)

gtf_new <- rbind(gtf_mono,gtf_multi)

chr_order=unique(ref_gtf$seqid)

gtf_StemLinc <- gtf_new %>% arrange(match(seqid,chr_order),start)

# export novel transcriptome ----
export(gtf_StemLinc,"outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/novel_StemLinc.gtf",format = "gtf")

# Merge with GENCODE VM36 ----
reference_gtf <- import(ref_gtf_path)

# Load novel transcriptome
new_gtf <- import("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/novel_StemLinc.gtf")

merged_transcriptome <- c(reference_gtf, new_gtf)

# export merged transcriptome ----
export(merged_transcriptome, "outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/novel_StemLinc.mergedVM36.gtf", format = "gtf")

# Write intermediate tables ----
write.table(filtered_multiexons_bed,"outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation/filtered_multiexons_bed.tsv",
            quote = F,row.names = F,sep = "\t")
write.table(filtered_monoexons,"outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation/filtered_monoexons.tsv",
            quote = F,row.names = F,sep = "\t")
# Write transcript and gene annotation files ----
# chr,start,end,gene_name,exonic_length,strand
old_raw_gene_annot=read.delim("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation_1/all_merged_genes_annotation.tsv")
gtf <- readGFF("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/novel_StemLinc.mergedVM36.gtf")
colnames(raw_gene_annot)
exonic_length <- trastools::compute_gene_non_redundant_exonic_length(gtf)
exons_per_tr <- gtf %>% filter(type=="exon") %>% group_by(gene_id,transcript_id) %>%
  summarise(N_exons=n())
exons_per_gene <- gtf %>% filter(type=="exon") %>%
  mutate(exon_id=paste(seqid,":",start,"-",end,":",strand)) %>% group_by(gene_id) %>%
  summarise(unique_exons=length(unique(exon_id)))
max_num_exons <- exons_per_tr %>% group_by(gene_id) %>% summarise(max_num_exons=max(N_exons))
raw_gene_annot <- gtf%>% filter(type=="exon") %>% group_by(gene_id) %>% summarise(chr=seqid[1],
                                                         start=min(start),
                                                         end=max(end),
                                                         strand=strand[1],
                                                         gencode_biotype=gene_type[1],
                                                         gencode_gene_name=gene_name[1])


colnames(old_raw_gene_annot)
table(is.na(raw_gene_annot$gencode_biotype))
# simplify ref biotypes ----
ref_biotypes=unique(raw_gene_annot$gencode_biotype)
#simplify biotypes
simpl_biotypes <- list("protein_coding" = c("protein_coding",
                                            grep("_gene|_V_",
                                                 grep("pseudogene",
                                                      ref_biotypes,
                                                      invert = T,value = T),value = T)),
                       "lncRNA" = "lncRNA",
                       "pseudogene" = grep("pseudogene",
                                           ref_biotypes,
                                           value = T),
                       "TEC" = grep("TEC",
                                    ref_biotypes,
                                    value = T),
                       "miRNA" =grep("miRNA",
                                     ref_biotypes,
                                     value = T),
                       "sno-snRNA" = c("snoRNA"  ,
                                       "snRNA",
                                       "ncRNA;snoRNA"),
                       "rRNA" = grep("rRNA",
                                     ref_biotypes,
                                     value = T),
                       "tRNA" = grep("tRNA",
                                     ref_biotypes,
                                     value = T))
#"other" = )

table(ref_biotypes%in%unlist(simpl_biotypes))

ref_biotypes[!ref_biotypes%in%unlist(simpl_biotypes)]

simpl_biotypes <- c(simpl_biotypes,
                    list(other = unique(ref_biotypes)[!unique(ref_biotypes)%in%
                                                        unlist(simpl_biotypes)&!is.na(unique(ref_biotypes))]))
simpl_biotypes_table <- data.frame(orig.biot=unlist(simpl_biotypes),
                                   final.biot=rep(names(simpl_biotypes),sapply(simpl_biotypes,length)))

write.table(simpl_biotypes_table,"data/references/gencode_vM36/simplified_biotypes_table.tsv",
            quote = F,row.names = F,sep = "\t")
raw_gene_annot$biotype=simpl_biotypes_table$final.biot[match(raw_gene_annot$gencode_biotype,
                                                           simpl_biotypes_table$orig.biot)]

raw_gene_annot$biotype[is.na(raw_gene_annot$biotype)] <- "potNovel"
table(raw_gene_annot$biotype)

raw_gene_annot <- left_join(raw_gene_annot,exons_per_gene)
raw_gene_annot <- left_join(raw_gene_annot,max_num_exons)
raw_gene_annot <- left_join(raw_gene_annot,exonic_length)


# assign unique gene names
raw_gene_annot <- raw_gene_annot %>% mutate(gene_name=ifelse(is.na(gencode_gene_name),gene_id,gencode_gene_name))
duplicated_gene_names <- raw_gene_annot$gene_name[duplicated(raw_gene_annot$gene_name)]
raw_gene_annot <- raw_gene_annot %>%mutate(gene_name=ifelse(gene_name%in%duplicated_gene_names,
                                                            paste(gene_name,gene_id,sep = "_"),
                                                            gene_name))

biots=c("protein_coding","lncRNA","potNovel","pseudogene","TEC")
write.table(raw_gene_annot,"outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation/all_merged_genes_annotation.tsv",
            quote = F,row.names = F,sep = "\t")

write.table(raw_gene_annot %>% filter(biotype%in%biots),"outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation/all_merged_genes_annotation.sel_biotypes.tsv",
            quote = F,row.names = F,sep = "\t")

old_raw_transcript_annot <- read.delim("outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation_1/all_merged_transcripts_annotation.tsv")
colnames(old_raw_transcript_annot)
raw_transcript_annot <- gtf %>%filter(type=="transcript") %>% select(transcript_id, seqid,start,end,strand,
                                                                     gene_id)
raw_transcript_annot <- left_join(raw_transcript_annot, exons_per_tr)
raw_transcript_annot <- left_join(raw_transcript_annot,raw_gene_annot %>%select(gene_id,gene_name,biotype))
exonic_type <- exons_per_tr %>%group_by(gene_id) %>% summarise(exonic_type=ifelse(any(N_exons>1),
                                                                                  "multiexonic",
                                                                                  "monoexonic"))

raw_gene_annot <- left_join(raw_gene_annot,exonic_type)

write.table(raw_transcript_annot,"outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation/all_merged_transcripts_annotation.tsv",
            quote = F,row.names = F,sep = "\t")

write.table(raw_transcript_annot %>% filter(biotype%in%biots),"outputs/transcriptome_characterization/filtered_LSK_T-cell_macrophage.combined/feb25/merged_transcriptome_annotation/all_merged_transcripts_annotation.sel_biotypes.tsv",
            quote = F,row.names = F,sep = "\t")
