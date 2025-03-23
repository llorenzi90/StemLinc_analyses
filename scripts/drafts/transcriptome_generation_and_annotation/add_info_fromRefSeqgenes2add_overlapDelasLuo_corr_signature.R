## Head -------------------------------------
##
##
## Purpose of script: collect info for known and potential novel lncRNAs 
## expressed in LSK StemLinc samples
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-04-14
##
## Email: lucialorenzi90@gmail.com
##
## Notes ---------------------------
##
## Things to check: 
##                  - Number of exons (for known ones, the one in annotation)
##                  - coding potential
##                  - distance to closest PCG and PCG id
##                  - distance to closest lncRNA and potnew lncRNA
##                  - correlation with closest PCG
##                  - overlap with CAGE, TFs, Enhancer Atlas
##                  - classcode
##                  - GO enrichment info
##                  - evidence within StemLinc dataset: transcript level/gene level
##                  - evidence from other datasets:
##                          - N samples with overlap (gene level)
##                          - N samples with classcode match (tr level)
##                          - Exonic type in public data
##                  - For gencode lncRNAs:
##                          - N StemLinc samples overlap (0-3)
##                          - N StemLinc samples exact match (0-3)
##                          - N public samples assembled (0-84)
##                          - Exonic type in StemLinc assemblies
##                          - level of evidence in GENCODE
##                  - expression across LSK samples
##                  - highest expression across other samples
##                  - highest public sample
##                  - orientation with respect to closest PCG (intergenic, 
##                    divergent, convergent, antisense, antisense to lncRNA or other)
##
## Setup ---------------------------

options(scipen = 999) 
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
require(rtracklayer)
## Load data---------------------------
setwd("~/work_local/APRIL/")
source("~/work_local/scripts/get_genomic_regions_function.R")
# write data function
wdat <- function(gene_level_info){
  write.table(gene_level_info,"~/work_local/APRIL/gene_level_info_APRIL.txt",
            sep = "\t",quote = F,row.names = F)
  }

# read GTF merged known and novel ----
GTF=readGFF("RefSeq_finaltranscripts2add.gtf")

# read info on lncRNAs and PCGs
#simpl_gencode_genes_info <- read.table('',header = T)


N_exons=GTF%>% filter(type=="exon") %>% group_by(transcript_id) %>% 
  summarise(N_exons=n())
N_exons$gene_name=GTF$gene_name[match(N_exons$transcript_id,
                                  GTF$transcript_id)]
gene_level_info=N_exons%>% group_by(gene_name) %>% 
  summarise(N_exons_tr=paste(N_exons,collapse = ","),
            exonic_type=ifelse(any(N_exons!="1"),"multiexonic","monoexonic"),
            N_transcripts=n())

gene_bed=get_genomic_regions(GTF,by = "gene_name")
gene_bed <- left_join(gene_bed, GTF%>% filter(!duplicated(gene_id))%>% dplyr::select(gene_name,gene_id,biotype))
# add  coords, biotype, gene_name ----
gene_level_info <- left_join(gene_level_info,gene_bed)

colnames(gene_level_info)[c(1,10)]=c("gene_id","gene_name")
gene_level_info$gene_name=gene_level_info$gene_id

# Add exonic_length ----
source("/home/llorenzi/work_local/scripts/compute_gene_non_redundant_exonic_length.function.R")
exonic_length <- compute_gene_non_redundant_exonic_length(in_gtf = GTF)
exonic_length$gene_id=GTF$gene_name[match(exonic_length$gene_id,
                                          GTF$gene_id)]
gene_level_info$width <- exonic_length$exonic_length[match(gene_level_info$gene_id,
                                                           exonic_length$gene_id)]

colnames(gene_level_info)[8] <- "exonic_length"
#wdat(gene_level_info)

#gene_info=read.table("~/work_local/potnovel_lncRNAs_info/gene_level_info_gencodeANDpotNovelStemLinc_lncRNAs.txt",sep = "\t",header = T)
# coding potential CPAT ----
CPAT_gene_level=read.table("CPAT/bestORF_gene_strand_RefSeq_finaltranscripts2add.gene_level.txt",
                           header = T)

# a few gene_ids are duplicated in this table, this is because some 
# genes have transcripts with the same ORFs
# collapse those
dupgeneids=CPAT_gene_level$gene_id[duplicated(CPAT_gene_level$gene_id)]
View(CPAT_gene_level%>%filter(gene_id%in%dupgeneids))
CPAT_gene_level_summ <- CPAT_gene_level %>% group_by(gene_id) %>% 
  summarise(ID_collapsed=paste(ID, collapse = ","),
            transcript_id_collapsed=paste(transcript_id,collapse = ","))
CPAT_gene_level <- left_join(CPAT_gene_level,CPAT_gene_level_summ)
CPAT_gene_level <- CPAT_gene_level %>% filter(!duplicated(gene_id))

gene_level_info <- left_join(gene_level_info , CPAT_gene_level%>% dplyr::select(gene_id,Coding_prob))

# distance to closest genes ----
closest_genes=read.table("closest_genes/knownANDpotNovel_lncRNAs_pseudo_closest.txt",header = T)
table(closest_genes$t_biotype)
#modify biotypes 
closest_genes$t_biotype[closest_genes$t_biotype%in%c("lncRNA","TEC")]="lncRNA"
closest_genes$t_biotype[grep("pseudo",closest_genes$t_biotype)]="pseudogene"

closest_genes$t_biotype[!closest_genes$t_biotype%in%c("lncRNA","potnew_lncRNA","pseudogene")]="PCG"


closest_genes_summary <- closest_genes %>% group_by(q_gene_id,t_biotype) %>% 
  summarise(closest_gene=paste(unique(t_gene_name[distance==min(distance)]),collapse = ","),
            distance=min(distance))

closest_genes_summary_clgene=closest_genes_summary%>%
  dplyr::select(q_gene_id,t_biotype,closest_gene,distance) %>% 
  pivot_wider(names_from = t_biotype,values_from = c(closest_gene,distance))


# Merge new info with previous info ----
# At this point I want to add the info on previous lncRNAs
# For now, I don't need to include all info 
# that I collected before because I want to focus only
# on the lncRNAs that are expressed/relevant for primary samples
# I want to clean up a bit before continuing with more analyses

# Read info on previous lncRNAs
prev_gene_info=read.table("/home/llorenzi/work_local/potnovel_lncRNAs_info/gene_level_info_gencodeANDpotNovelStemLinc_lncRNAs.txt",sep = "\t",header = T)
prev_gene_info_12=prev_gene_info[,1:12]
dim(prev_gene_info_12)

colnames(prev_gene_info_12)%in%colnames(gene_level_info)
colnames(prev_gene_info_12)[11]="biotype"

table(gene_level_info$gene_id%in%prev_gene_info$gene_name)

gene_level_info=rbind(gene_level_info,prev_gene_info_12)

gene_level_info <- left_join(gene_level_info,
                             closest_genes_summary_clgene%>%
                               dplyr::select(q_gene_id,
                                      closest_gene_PCG,
                                      distance_PCG,
                                      closest_gene_lncRNA,
                                      distance_lncRNA,
                                      closest_gene_potnew_lncRNA,
                                      distance_potnew_lncRNA,
                                      closest_gene_pseudogene,
                                      distance_pseudogene),
                             by=c(gene_id="q_gene_id"))


nrow(gene_level_info)-nrow(prev_gene_info_12)
prevrows=9160:nrow(gene_level_info)
# Conservation scores ----

phastCons35way <- read.table("phastCons35way/RefSeq_finalgenes2add.phastCons35way.txt",header = T)

phastCons35way$gene_id=phastCons35way$gene_name

phastCons35way=phastCons35way%>%group_by(gene_id)%>%
  summarise(meanPhastCons35way=mean(mean_score,na.rm=T))

gene_level_info <- left_join(gene_level_info,phastCons35way %>%dplyr::select(gene_id,meanPhastCons35way))


gene_level_info$meanPhastCons35way[is.na(gene_level_info$meanPhastCons35way)]=prev_gene_info$meanPhastCons35way[match(gene_level_info$gene_id[is.na(gene_level_info$meanPhastCons35way)],
                                                                                                                      prev_gene_info$gene_id)]
#wdat(gene_level_info)


# overlap with CAGE, TFs, Enhancer Atlas ----
ol_marks=read.table("/home/llorenzi/work_local/overlap_marks/RefSeq2add_overlapmarks.gene_level.txt",header = T)
colnames(ol_marks)[4]="gene_id"

gene_level_info <- left_join(gene_level_info, 
                             ol_marks%>% dplyr::select(gene_id,closest_CAGE,TF_rel.lev_A,TF_rel.lev_B,TF_rel.lev_C,TF_rel.lev_D,dist_closest_enhancer,closest_polyA))

colnames(gene_level_info)
gene_level_info[prevrows,22:28] <- prev_gene_info[,25:31]
wdat(gene_level_info = gene_level_info)


# RefSeq/gencode overlap with StemLinc ----
ol_RefSeq <- read.table('/home/llorenzi/work_local/APRIL/RefSeq_ol_StemLinc/RefSeq2add_overlap_StemLinc.gene_level.txt',header = T)
colnames(ol_RefSeq)[4:8] <- c("NStemLinc_exact_match",
                               "NStemLinc_other_match",
                               "NStemLinc_any_match",
                               "NStemLinc_antisense",
                               "max_Nexons_StemLinc")
colnames(ol_RefSeq)[1]="gene_id"

ol_refgenes_vs_StemLinc=rbind(ol_RefSeq %>% dplyr::select(gene_id,NStemLinc_exact_match,NStemLinc_other_match,NStemLinc_any_match,NStemLinc_antisense,max_Nexons_StemLinc),
                              prev_gene_info%>%dplyr::select(gene_id,NStemLinc_exact_match,NStemLinc_other_match,NStemLinc_any_match,NStemLinc_antisense,max_Nexons_StemLinc))

gene_level_info=read.table("~/work_local/APRIL/gene_level_info_APRIL.txt",header = T,sep = "\t")
gene_level_info <- left_join(gene_level_info,ol_refgenes_vs_StemLinc)
                             
wdat(gene_level_info = gene_level_info )                             

# Add overlap with Delas and Luo ----
### Delas ----
Delas_ol=read.table("/home/llorenzi/work_local/APRIL/Delas_ol_classcodes_relative_to_RefStemLinc_genes.txt",header = T)

# add annotation on Delas genes
classes=c("AML_enriched","HSC_enriched","lymphoid_enriched","progenitor_vs_differentiated")

Delas_classes=lapply(1:4,function(x)readxl::read_xls("/home/llorenzi/work_local/public_datasets/Delas_et_al/elife-25607-supp2-v2.xls",sheet = x,col_names = F))
names(Delas_classes)=classes
Delas_classes=lapply(names(Delas_classes),function(x){
  dc=Delas_classes[[x]]
  colnames(dc)="Delas_gene"
  dc=unlist(strsplit(dc$Delas_gene,split = ","))
  return(data.frame(Delas_gene=dc,
                    class=x))
})

Delas_classes=do.call("rbind",Delas_classes)
length(unique(Delas_classes$Delas_gene))
Delas_classes=Delas_classes%>%group_by(Delas_gene)%>%
  summarise(class=paste(class,collapse = ","))

# add annotation to overlap
Delas_ol=left_join(Delas_ol,Delas_classes)


Delas_5_selected=read.csv("/home/llorenzi/work_local/public_datasets/Delas_et_al/Delas_5_manual_selected.csv")

# For these 5 lncRNAs add extra info

Delas_ol$Delas_annot=Delas_5_selected$Delas_annot[match(Delas_ol$Delas_gene,
                                                        Delas_5_selected$Delas.2017.Annotation)]

colnames(Delas_ol)[5]="Delas_enrichment"

# separate in overlap and antisense
table(Delas_ol$best.cc.Delas2Ref)
overlapping_ccs=c("=","c","j","k","m","n","o","p")
table(Delas_ol$best.cc.Delas2Ref[Delas_ol$best.cc.Delas2Ref%in%overlapping_ccs])
table(Delas_ol$best.cc.Delas2Ref[!Delas_ol$best.cc.Delas2Ref%in%overlapping_ccs])

Delas_non_ol=Delas_ol[!Delas_ol$best.cc.Delas2Ref%in%overlapping_ccs,]
Delas_ol=Delas_ol[Delas_ol$best.cc.Delas2Ref%in%overlapping_ccs,]

Delas_ol <- Delas_ol%>% group_by(Ref_gene)%>%
  summarise(Delas_overlap=paste0(Delas_gene,collapse = ","),
            classcodeDelas=paste0(best.cc.Delas2Ref,collapse = ","),
            Delas_enrichment=paste(Delas_enrichment,collapse = ","),
            Delas_annot=paste(Delas_annot,collapse = ","))

Delas_non_ol <- Delas_non_ol%>% group_by(Ref_gene)%>%
  summarise(Delas_other=paste0(Delas_gene,collapse = ","),
            otherclasscodeDelas=paste0(best.cc.Delas2Ref,collapse = ","),
            otherDelas_enrichment=paste(Delas_enrichment,collapse = ","),
            otherDelas_annot=paste(Delas_annot,collapse = ","))


gene_level_info=left_join(gene_level_info,Delas_ol,by=c(gene_id="Ref_gene"))
gene_level_info=left_join(gene_level_info,Delas_non_ol,by=c(gene_id="Ref_gene"))
wdat(gene_level_info)

### Luo ----
luomm39=read.table("/home/llorenzi/work_local/public_datasets/Goodell/annotations/Goodell_159_LncHSCs.liftOver.mm39.bed")
luomm39_sorted=luomm39%>% arrange(V1,V2)

gene_level_info_bed=gene_level_info%>%dplyr::select(seqnames,
                                                    start,
                                                    end,
                                                    gene_id,
                                                    exonic_length,
                                                    strand)

gene_level_info_bed=gene_level_info_bed%>%arrange(seqnames,start)
luo_ol=bedtoolsr::bt.intersect(a=luomm39,b=gene_level_info_bed,wo=T, 
                              s=T)
luo_ol=luo_ol%>%group_by(V10)%>%summarise(Luo_gene=paste(V4,collapse = ","))
colnames(luo_ol)[1]="gene_id"
gene_level_info=left_join(gene_level_info,luo_ol)


delasvsluo=as.data.frame(table(gene_level_info$Delas_overlap,gene_level_info$Luo_gene))
delasvsluo=delasvsluo%>%filter(Freq!=0)
wdat(gene_level_info = gene_level_info )                             

gene_level_info=read.delim("/home/llorenzi/work_local/APRIL/gene_level_info_APRIL.txt")

# Classify lncRNAs ----

gene_info=read.table("../potnovel_lncRNAs_info/gene_level_info_gencodeANDpotNovelStemLinc_lncRNAs.txt",header = T,sep = "\t")
gene_level_info=read.table("gene_level_info_APRIL.txt",header = T)
lncRNA_classes=read.table("gene_classif/lncRNA_classification.txt",header = T)
colnames(lncRNA_classes)[4]="ol_gene"
gene_level_info=left_join(gene_level_info,lncRNA_classes%>%dplyr::select(lncRNA_id,class,ol_gene,best_class),
                          by=c(gene_id="lncRNA_id"))
wdat(gene_level_info = gene_level_info )

# add info on network analysis ----

intramodular_connectivity <- read.delim("~/work_local/APRIL/WGCNA/intramodular_connectivity.txt")
gene_level_info=left_join(gene_level_info,
                          intramodular_connectivity %>%
                            dplyr::select(gene,is_primary_corr,blue_module_membership,
                                          assigned_module_membership,kWithin,kTotal,assigned_module),
                          by=c(gene_id="gene"))

wdat(gene_level_info)

# Repeats ----
library(rtracklayer)

#repeats=read.table("/home/llorenzi/work_local/UCSC_RepeatMasker_mm39.bed")
repeats=read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/RepeatMasker/UCSC_RepeatMasker_mm39.bed")
gencodelncRNA_gtf=readGFF("/home/llorenzi/references/gencode.vM31.long_noncoding_RNAs.gtf")
RefSeq_gtf=readGFF("/home/llorenzi/work_local/APRIL/RefSeq_finaltranscripts2add.gtf")
StemLincPotNov_gtf=readGFF("/home/llorenzi/work_local/LSK_StemLinc.combined.potNovel_3reps.gtf")
table(gencodelncRNA_gtf$gene_id%in%gene_level_info$gene_id)
#source("/home/llorenzi/work_local/scripts/new_scripts/compute_gene_non_redundant_exons_and_exonic_length.functions.R")
source("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/scripts/StemLinc_code/new_scripts/compute_gene_non_redundant_exons_and_exonic_length.functions.R")
get_overlap_with_repeats_from_GTF=function(gtf){
  require(bedtoolsr)
  require(tidyverse)
  pgme=get_per_gene_merged_exons(gtf)
  pgme=pgme[,c(3:5,2,6,7)]
  pgme=pgme[order(pgme$seqnames,pgme$start),]
  exonic_length=pgme%>%group_by(group_name)%>%
    summarise(exonic_length=sum(width))
  ol=bt.intersect(a = pgme,b=repeats,wo = T)
  ol=ol%>% rowwise()%>% mutate(st=max(c(V2,V8)),en=min(c(V3,V9)))
  olsm=ol%>%group_by(V4)%>%
    summarise(repeats=paste(unique(V10[order(-V13)]),
                            collapse = ","),
              repeats_classes=paste(unique(V11[order(-V13)]),
                                    collapse = ","),
              t.overlap_length=sum(IRanges::width(IRanges::reduce(IRanges(
                start = st,end = en)))))
  olsm=left_join(olsm,exonic_length,by=c(V4="group_name"))
  olsm$repeat.fraction=olsm$t.overlap_length/olsm$exonic_length
  return(olsm)
}


#repeats_chr=unique(repeats$V1)
repeats_chr=unique(repeats$V6)

# repeats=repeats[order(match(repeats$V1,repeats_chr),
#                       repeats$V2),]

repeats=repeats[order(match(repeats$V6,repeats_chr),
                      repeats$V7),]

repeats=repeats[,c(6,7,8,11,12,10)]
gencodelncRNA_ol=get_overlap_with_repeats_from_GTF(gencodelncRNA_gtf)
StemLincPotNov_ol=get_overlap_with_repeats_from_GTF(StemLincPotNov_gtf)
RefSeq_gtf$gene_id=RefSeq_gtf$gene_name[match(RefSeq_gtf$gene_id,
                                              RefSeq_gtf$gene_id)]
RefSeq_ol=get_overlap_with_repeats_from_GTF(RefSeq_gtf)

ol_repeats=rbind(StemLincPotNov_ol,
                 gencodelncRNA_ol,
                 RefSeq_ol)
names(ol_repeats)[1]="gene_id"
write.table(ol_repeats,"overlap_repeats_all_noncoding_genes.txt",sep = "\t",row.names = F)
ol_repeats=read.table("/home/llorenzi/work_local/APRIL/overlap_repeats_all_noncoding_genes.txt",sep = "\t",header = T)

gene_level_info=read.table("/home/llorenzi/work_local/APRIL/gene_level_info_APRIL.txt",header = T)

gene_level_info=left_join(gene_level_info,
                          ol_repeats%>%dplyr::select(gene_id,repeats,repeat.fraction))
gene_level_info$repeat.fraction[is.na(gene_level_info$repeat.fraction)]=0

wdat(gene_level_info)

# Correlation with HSC signature ----
vst=read.table("/home/llorenzi/work_local/expression_files/featureCounts/PCGslncRNAsPseudogenes.110424.vst.txt")
gene_biot=read.table("/home/llorenzi/work_local/expression_files/featureCounts/PCGslncRNAsPseudogenes.110424.annot.txt")
table(gene_biot$simpl_biotype[!gene_biot$Geneid%in%gene_level_info$gene_id])
gene_translation_file <- "/home/llorenzi/references/conversion_table_ENSEMBL_NCBI_ids.csv"
gene_translation <- fread(gene_translation_file)
gene_biot$gene_name=gene_translation$gene_name[match(gene_biot$Geneid,
                                                     gene_translation$gene_id)]
vst_lncRNAs=vst[rownames(vst)%in%gene_biot$Geneid[gene_biot$simpl_biotype!="PCG"],]

gene_signature=read.table("/home/llorenzi/work_local/HSC_signature.txt")
#gene_signature=read.table("HSC_signature_Sergi.txt")
gene_sign_ids=gene_biot$Geneid[match(gene_signature$V1,gene_biot$gene_name)]
vst_signature=vst[rownames(vst)%in%gene_sign_ids,]

corr=cor(t(vst_lncRNAs),t(vst_signature))

corr=as.data.frame(corr)
corr$lncRNA=rownames(corr)
corr=pivot_longer(corr,cols=1:(ncol(corr)-1),
                  values_to = "corr",names_to = "PCG_id")
corr$PCG_name=gene_translation$gene_name[match(corr$PCG_id,gene_translation$gene_id_vM31)]
corr$lncRNA_name=gene_level_info$gene_name[match(corr$lncRNA,gene_level_info$gene_id)]

corr_summ=corr %>% arrange(-corr)%>%group_by(lncRNA)%>% 
  filter(any(corr>=0.8)) %>%
  summarise(N_0.8=sum(corr>=0.8),
            N_0.65=sum(corr>=0.65),
            PCGs=paste(PCG_name,collapse = ","),
            corr=paste(round(corr,2),collapse = ","))

corr_summ=left_join(corr_summ,gene_level_info%>%dplyr::select(gene_id,closest_gene_PCG,distance_PCG),
                    by=c(lncRNA="gene_id"))

corr_summ$lncRNA_name=gene_level_info$gene_name[match(corr_summ$lncRNA,
                                                      gene_level_info$gene_id)]

corr_summ=corr_summ%>%arrange(-N_0.8,-N_0.65)


corr_summ=left_join(corr_summ,intramodular_connectivity %>%
                      dplyr::select(gene,is_primary_corr,blue_module_membership,
                                    assigned_module_membership,kWithin,kTotal,assigned_module),
                    by=c(lncRNA="gene"))

corr_summ_extended=left_join(corr_summ,gene_level_info%>%
                      dplyr::select(gene_id,exonic_type,
                                    Coding_prob,
                                    meanPhastCons35way,
                                    closest_CAGE,
                                    dist_closest_enhancer,
                                    Delas_overlap,
                                    Delas_enrichment,
                                    Luo_gene),by = c(lncRNA="gene_id"))

counts=read.table("/home/llorenzi/work_local/expression_files/featureCounts/PCGslncRNAsPseudogenes.110424.counts.txt")
mean_counts_LSK_StL=apply(counts[,grep("LSK_StL",colnames(counts))],1,mean)

corr_summ_extended$mean_counts_LSK_StL=mean_counts_LSK_StL[match(corr_summ_extended$lncRNA,
                                                                 names(mean_counts_LSK_StL))]

mean_corr_signature=corr%>%group_by(lncRNA)%>%summarise(mean_corr=mean(corr))
mean_corr_signature=mean_corr_signature%>%arrange(-mean_corr)

View(mean_corr_signature%>%filter(mean_corr>0.645))
corr_summ_extended$mean_corr=mean_corr_signature$mean_corr[match(corr_summ_extended$lncRNA,
                                                                 mean_corr_signature$lncRNA)]
gene_level_info=left_join(gene_level_info,mean_corr_signature,by=c(gene_id="lncRNA"))
colnames(gene_level_info)[colnames(gene_level_info)=="mean_corr"]="mean_corr_signature"
wdat(gene_level_info)
#gene_info=read.table("potnovel_lncRNAs_info/gene_level_info_gencodeANDpotNovelStemLinc_lncRNAs.txt",header = T,sep = "\t")

# transcript level support
prev_gene_info=read.table("/home/llorenzi/work_local/potnovel_lncRNAs_info/gene_level_info_gencodeANDpotNovelStemLinc_lncRNAs.txt",header = T,sep = "\t")

gene_level_info=left_join(gene_level_info,prev_gene_info%>%dplyr::select(gene_id,tr_level_support))
wdat(gene_level_info)

length(unique(corr$lncRNA[corr$corr>=0.8]))
table(grepl("XLOC",unique(corr$lncRNA[corr$corr>=0.8])))
potnov=grep("XLOC",unique(corr$lncRNA[corr$corr>=0.8]),value = T)

gene_info18=read.table("potnovel_lncRNAs_info/gene_level_info_18potNovel.txt",sep = "\t",header = T)
table(potnov%in%gene_info18$lncRNA)

potnov0.8=corr%>%filter(lncRNA%in%potnov)
table(potnov0.8$lncRNA,potnov0.8$corr>=0.8)
potnov[potnov%in%gene_info18$lncRNA]

for (cutoff in c(0.65,0.7,0.8)) {
  # how many 
  
  print("######################")
  print(paste0("cutoff: ",cutoff))
  print("")
  print(paste0("lncRNAs with at least 1 correlated signature PCG >= ", cutoff, ":"))
  print(length(unique(corr$lncRNA[corr$corr>=cutoff])))
  print("Of which potential novel:")
  print(table(grepl("XLOC",unique(corr$lncRNA[corr$corr>=cutoff]))))
  print(paste0("potential novel lncRNAs with N PCgenes with correlation >= ",cutoff))
  print(table(table(grep("XLOC",corr$lncRNA[corr$corr>=cutoff],value = T))))
  print(paste0("Annotated lncRNAs with N PCgenes with correlation >= ",cutoff))
  print(table(table(grep("XLOC",corr$lncRNA[corr$corr>=cutoff],value = T,invert = T))))
  print("")
}

# add intramodular connectivity 

# Check expression values ----
# remove genes that are virtually non-expressed
# classify genes in terms of position relative to
# PCGs
# Check correlation with signature PCGs, compare known and novel

# perform DGEA and check if extra interesting genes arise
# play around with WGCNA

TPMS=read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/featureCounts/PCGslncRNAsANDpotNovelLSK.TPM.txt")

TPMS <- TPMS[rownames(TPMS)%in%gene_level_info$gene_id,]
meanStemLinc=apply(TPMS[,1:3],1,mean)

max_StemLinc=apply(TPMS[,1:3],1,max)

mean_othersamps=apply(TPMS[,4:ncol(TPMS)],1,mean)
max_othersamps=apply(TPMS[,4:ncol(TPMS)],1,max)
max_samp_othersamps=apply(TPMS[,4:ncol(TPMS)],1,function(r){
  colnames(TPMS)[4:ncol(TPMS)][which.max(r)]})

summary_TPMS=cbind(gene_id=names(meanStemLinc),
                   meanStemLincTPM=round(meanStemLinc,2),
                   maxStemLincTPM=round(max_StemLinc,2),
                   meanOthersampsTPM=round(mean_othersamps,2),
                   maxOthersampsTPM=round(max_othersamps,2),
                   maxOthersamp=max_samp_othersamps)

summary_TPMS=as.data.frame(summary_TPMS)
gene_level_info <- left_join(gene_level_info,summary_TPMS)

wdat(gene_level_info = gene_level_info)

summary(as.numeric(gene_level_info$meanStemLincTPM))

View(gene_level_info %>% arrange(-as.numeric(meanStemLincTPM)))


# class code and transcript level support for novel ones ----
tracking=read.table('/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare/LSK_StemLinc.tracking')
tracking=tracking %>% filter(V1%in%GTF$transcript_id)
tracking=separate(tracking,col = V3,sep = "\\|",into = c("ref_gene_id","ref_transcript_id"))
ref_id_biot=read.table('/home/llorenzi/shares/INVESTIGACIO/Cytuartero Group/CUARTERO GROUP/references/mouse/annotation/merged_refs_annotation/merged_refs_gene_biotypes.txt',header = T)
tracking <- left_join(tracking,ref_id_biot,by=c(ref_gene_id="merged_gene_id"))

tracking$Nsamples=apply(tracking[,6:8],1,function(x)sum(x!="-"))

tracking_gl_cc=tracking %>% group_by(V2) %>% 
  summarise(classcodes=paste(V4,collapse = ","),
            ref_genes=paste(gene_name,collapse = ","),
            ref_biotype=paste(gene_biotype,collapse = ","),
            tr_level_support=max(Nsamples))
colnames(tracking_gl_cc)[1]="gene_id"
gene_level_info <- left_join(gene_level_info,tracking_gl_cc)

# evidence from other datasets ----
ol_pub_datasets <- read.table('/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/overlap_with_public_data/StemLin_potNovel_vs_potNovelpublicdataassemblies.gene_level.txt',header = T)
colnames(ol_pub_datasets)[c(3:6)] <- c("type_ol_pubdatasets",
                                       "pubdatasets_genes",
                                       "maxNexons_pubdatasets",
                                       "maxNsamps_pubdatasets")
gene_level_info <- left_join(gene_level_info,ol_pub_datasets%>%dplyr::select(gene_id,type_ol_pubdatasets,maxNexons_pubdatasets,maxNsamps_pubdatasets))






# correlation with PCGs ----
#corr_0.8=read.table("Experimento4_LSK_190224/analyses/correlation_PCG_lncRNA/corr_lncRNA_PCG.0.8.txt")
corr=fread("Experimento4_LSK_190224/analyses/correlation_PCG_lncRNA/corr_lncRNA_PCG.txt")

# correlation closest PCG ----
# get gene_id for closest PCG for each lncRNA
closest_genes_gene_id <- closest_genes %>% filter(t_biotype=="PCG") %>%group_by(q_gene_id) %>%
  summarise(closest_PCG=paste(t_gene_id[distance==min(distance)],collapse = ","))

closest_genes_gene_id <- closest_genes_gene_id %>% filter(q_gene_id %in% colnames(corr))

closest_genes_gene_id$Ngenes=sapply(strsplit(closest_genes_gene_id$closest_PCG,split = ","),length)
closest_genes_gene_id=as.data.frame(closest_genes_gene_id)
corr=as.data.frame(corr)

table(unique(unlist(strsplit(closest_genes_gene_id$closest_PCG,split = ",")))%in%corr$V1)

get_corr_val=function(mytrip, corr){
  
  lnc=as.character(mytrip[1])
  if(as.numeric(mytrip[3])>1){
    pcg=unlist(strsplit(as.character(mytrip[2]),split = ","))
    
  } else {
    pcg=as.character(mytrip[2])
    
  }
  if(any(pcg%in%corr$V1)){
    return(paste(round(corr[match(pcg,corr$V1,nomatch = 0),colnames(corr)==lnc],2),collapse = ","))
  } else return(NA)
}

closest_genes_corr=c()
for(i in 1:nrow(closest_genes_gene_id)){
  closest_genes_corr=c(closest_genes_corr, get_corr_val(closest_genes_gene_id[i,],corr))
}

closest_genes_gene_id$corr_vals=closest_genes_corr
colnames(closest_genes_gene_id) <- c("gene_id",
                                     "closest_PCG",
                                     "N_closestPCG",
                                     "corr_closestPCG")
gene_level_info <- left_join(gene_level_info,closest_genes_gene_id)

# N corr >0.8 ----
cor_thres=0.8
corr <- pivot_longer(corr,2:ncol(corr),names_to = "lncRNA",values_to = "corr")
colnames(corr)[1]="PCG"
corr$PCG_name=simpl_gencode_genes_info$gene_name[match(corr$PCG,
                                                       simpl_gencode_genes_info$gene_id)]

corr=corr %>% arrange(-abs(corr))

N_highly_correlated_PCG=corr %>%
  group_by(lncRNA) %>% 
  summarise(N_cor0.8_PCG=sum(abs(corr)>cor_thres),
            top_corr_PCGs=paste(na.omit(PCG_name[1:5]),collapse = ","),
            top_corr_vals=paste(round(na.omit(corr[1:5]),2),collapse = ",")) 

gene_level_info <- left_join(gene_level_info,N_highly_correlated_PCG,
                             by=c(gene_id="lncRNA"))



# GO enrichment info ----
all_GO_enrich=fread("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/correlation_PCG_lncRNA/GO_analyses/GOfuncR_all_matches.txt")
GO_enrich=read.table('/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/correlation_PCG_lncRNA/GO_analyses/GOfuncR_all_matches.summary.txt',sep = "\t",header = T)

colnames(GO_enrich)[c(1,2,4,5,8,12)] <- c("gene_id",
                                          "N_GOmatchs",
                                          "N_GOsignif_p.val",
                                          "N_GOsignif_corrected.p.val",
                                          "total_GOgenes",
                                          "all_GOgenes")
gene_level_info <- left_join(gene_level_info,GO_enrich%>%select(gene_id,N_GOmatchs,N_GOsignif_p.val,N_GOsignif_corrected.p.val,most_signif_term_Ngenes,total_GOgenes,most_signif_term,genes_in_most_signif_term,most_signif_term_corr_vals,all_GOgenes,NOTE))




# modules from WGCNA ----
gene_level_info=read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/gene_level_info_gencodeANDpotNovelStemLinc_lncRNAs.txt",sep = "\t",header = T)

modules_WGCNA=read.table("Experimento4_LSK_190224/analyses/WGCNA/module.gene.mapping.txt",header = T)

gene_level_info <- left_join(gene_level_info, modules_WGCNA%>% dplyr::select(gene_id,bwnet.colors))

# calculate GO enrichment for genes excluded from prev analyses ----
library(GOfuncR)
library(data.table)
gene_translation_file <- "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/mouse/conversion_table_ENSEMBL_NCBI_ids.csv"
gene_translation <- fread(gene_translation_file)

gene_universe=rownames(TPMS)
gene_universe <- gene_translation$gene_name[match(gene_universe,
                                                  gene_translation$gene_id)]
gene_universe <- na.omit(gene_universe)

get_go_res=function(mod_genes){
  input_hyper_mouse = data.frame(gene_id=gene_universe, 
                                 is_candidate=ifelse(gene_universe%in%mod_genes,1,0))
  
  
  res_hyper_mouse = go_enrich(input_hyper_mouse, orgDb = "org.Mm.eg.db")
  go_res=res_hyper_mouse$results
  
  
  anno_genes = get_anno_genes(go_ids=go_res$node_id, 
                              genes=mod_genes,database = "org.Mm.eg.db")
  
  anno_genes$GO_name=go_res$node_name[match(anno_genes$go_id,
                                            go_res$node_id)]
  
  anno_genes <- anno_genes %>% filter(!GO_name%in%c("molecular_function",
                                                    "cellular_component",
                                                    "biological_process"))
  
  go_id_gene_ids=anno_genes %>% group_by(go_id) %>%
    summarise(genes=paste0(gene,collapse = ","),
              Ngenes=length(gene))
  
  
  go_res <- left_join(go_res,go_id_gene_ids,by=c("node_id" = "go_id"))
  
  return(go_res)
  
}

go_res=get_go_res(mod_genes =  unlist(strsplit(gene_level_info$top_corr_PCGs[gene_level_info$gene_id=="XLOC_002566"],split = ",")))


wdat(gene_level_info)


# Add info on correlation with primary samples  ----
isprimary_assoc=read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/WGCNA/isprimaryCorrANDbrownmembership.txt",header = T)

gene_level_info <- left_join(gene_level_info,
                             isprimary_assoc %>% dplyr::select(-c(assigned_module,biotype)),
                             by=c(gene_id="gene"))
wdat(gene_level_info)


# Intramodular connectivity ----
all_network_info=read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/WGCNA/info_on_intramodular_connectivity.txt",header = T)
gene_level_info=left_join(gene_level_info,all_network_info%>%dplyr::select(gene_id,kTotal,kWithin,kOut,kDiff))
wdat(gene_level_info)