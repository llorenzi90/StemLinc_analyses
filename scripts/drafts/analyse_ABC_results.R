## Purpose of script: analize results from ABC pipeline
##
## Author: Lucia Lorenzi
##
## Date Created: 2024-10-30
##
## Email: lucialorenzi90@gmail.com
##
#  Notes ---------------------------
##
##
##
#  Setup ---------------------------

options(scipen = 999)
# Load packages---------------------------
# load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
library(rtracklayer)
library(bedtoolsr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

# Define file paths---------------------------
#ABC_results_path <- "~/ABC-Enhancer-Gene-Prediction/results/LT.HSC_all_data/"
ABC_results_path <- "~/ABC-Enhancer-Gene-Prediction/results/LT.HSC_with_HiC.R1/"
res_files <- list.files(paste0(ABC_results_path,"Predictions"),full.names = T)
res_file <- grep("Full",
                 grep("tsv",res_files,
                      value = T),
                 value = T)
hematopoietic_gene_list_path <- "data/gene_lists/HSC_genes.match_filtered_genes_440.txt"

# Define the color palette
nejm_pal <- c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF")

# Load data---------------------------

# Load hematopoietic gene list
hematopoietic_genes <- read.table(hematopoietic_gene_list_path, header = FALSE, stringsAsFactors = FALSE)
colnames(hematopoietic_genes) <- "TargetGene"

# Load ABC results, promoters, and gene boundaries
full_predictions <- read.table(res_file,
                               sep = "\t",
                               header = T)
promoters <- read.table("outputs/ABC_files/LSK_CollapsedGeneBounds.mm39.TSS500bp.bed",header = F,sep = "\t")
genes <- read.table("outputs/ABC_files/LSK_CollapsedGeneBounds.mm39.bed",header = F,sep = "\t")

# Adjust column names
colnames(promoters) <- c("chr", "start", "end", "gene_name", "score", "strand", "gene_id","biotype")
colnames(genes) <- c("chr", "start", "end", "gene_name", "score", "strand", "gene_id","biotype")

colnames(full_predictions)[colnames(full_predictions)%in%colnames(genes)] <-
  paste0("Enh_",colnames(full_predictions)[colnames(full_predictions)%in%colnames(genes)])

# Perform bedtools intersections to identify overlapping regions
ol_prom_pred <- bt.intersect(a = promoters, b = full_predictions, wo = TRUE)
ol_gene_pred <- bt.intersect(a = genes, b = full_predictions, wo = TRUE)
colnames(ol_gene_pred) <- c(colnames(genes),
                                  colnames(full_predictions),
                            "overlap_length")

# Mark enhancer - gene self relationships
ol_gene_pred_self <- ol_gene_pred %>%filter(gene_name==TargetGene) %>%
  mutate(Gene_Enhc_id=paste(gene_name,name,sep = "_"))

# Add ID to full predictions
full_predictions <- full_predictions %>% mutate(Gene_Enhc_id=paste(TargetGene,name,sep = "_"))


# Add GeneEnhancer
Enhancer_annot <- ol_gene_pred %>% dplyr::select(gene_name,name, biotype)
Enhancer_annot <- Enhancer_annot[!duplicated(Enhancer_annot),]
Enhancer_annot <- Enhancer_annot %>% group_by(name) %>%
  summarise(EnhancerGene=paste0(sort(unique(gene_name)),collapse = ","),
            EnhancerBiotype=paste0(unique(sort(biotype)),collapse = ","))

# Add GeneEnhancer to full predictions
full_predictions <- left_join(full_predictions, Enhancer_annot)

# Remove self-enhancers (self-promoter or overlap with its own gene) and non-expressed targets
non_self_enhancers <- full_predictions%>%
  filter(isSelfPromoter == "False" & !Gene_Enhc_id%in%ol_gene_pred_self$Gene_Enhc_id
         &TargetGeneIsExpressed=="True")



length(unique(non_self_enhancers$TargetGene))

# Number of enhancer targets per biotype ----
non_self_enhancers <- left_join(non_self_enhancers,
                                genes %>%dplyr::select(gene_name,biotype),
                                by=c(TargetGene="gene_name"))

# Fraction of each biotype that is target of enhancer ----
# Load gene info ----
gene_level_info <- read.csv("outputs/gene_level_info_all_evidence_oct24.csv")

gene_level_info <- gene_level_info %>% mutate(Enhancer_target =
                                                gene_name%in%non_self_enhancers$TargetGene)

Gene_as_targets <-  non_self_enhancers %>% group_by(gene_name = TargetGene) %>%
  arrange(-ABC.Score) %>%
  summarise(N_enh_asTarget = n(),
            maxABCscore_asTarget = ABC.Score[1],
            maxEnh_asTarget = name[1],
            maxEnh_biotype_asTarget = EnhancerBiotype[1])

gene_level_info <- left_join(gene_level_info, Gene_as_targets)

Fraction_targets <- gene_level_info %>% group_by(biotype) %>%
  summarise(Fraction_targets=sum(Enhancer_target)/n())

ggplot(Fraction_targets, aes(x=biotype,y=Fraction_targets)) + geom_col()+
  labs(title = "Fraction of targets of enhancer per biotype",
       x = "biotype", y = "Fraction of enhancer targets")


N_enhancers <- non_self_enhancers %>% group_by(biotype,TargetGene) %>%
  summarise(N_enhancers=n()) %>%
  ungroup()



N_enhancers %>% group_by(biotype) %>% summarise(median=median(N_enhancers),
                                                mean=mean(N_enhancers))

# Histogram Number of enhancers per biotype ----
ggplot(N_enhancers, aes(x = N_enhancers)) +
  geom_histogram(binwidth = 1, fill = nejm_pal[1], color = "black") +
  facet_wrap(~biotype)+
  labs(title = "Number of Enhancers per Gene biotype",
       x = "Number of Enhancers", y = "Frequency")

# ABC scores per biotype ----
ggplot(non_self_enhancers, aes(x = biotype, y=ABC.Score)) +
  geom_boxplot() +

  labs(title = "ABC score per Target Gene biotype",
       x = "biotype", y = "ABC score") + scale_y_log10()


# Enhancer location per biotype ----

# Plot: Enhancer Location Distribution
ggplot(non_self_enhancers, aes(x = class)) +
  geom_bar(fill = nejm_pal[2]
           ) + facet_wrap(~biotype)+
  labs(title = "Location of Enhancers Relative to Genes per target biotype",
       x = "Location", y = "Interactions Count")

non_self_enhancers <- non_self_enhancers%>%
  mutate(EnhancerBiotype=ifelse(class=="intergenic",
                                "intergenic",
                                ifelse(grepl("protein_coding",EnhancerBiotype),
                                "protein_coding",
                                ifelse(grepl("lncRNA", EnhancerBiotype),
                                             "lncRNA",
                                       ifelse(grepl("potNovel",EnhancerBiotype),
                                              "potNovel",
                                              ifelse(grepl("TEC",EnhancerBiotype),
                                                     "TEC",EnhancerBiotype))))))

ggplot(non_self_enhancers, aes(x = EnhancerBiotype, fill = EnhancerBiotype)) +
  geom_bar() + facet_wrap(~biotype)+
  labs(title = "Location of Enhancers Relative to Genes per target biotype",
       x = "Location", y = "Interactions Count") +
  scale_fill_manual(values = nejm_pal) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggplot(non_self_enhancers,aes(x=EnhancerBiotype,fill=EnhancerBiotype,
                              y=ABC.Score))+
  geom_boxplot() + facet_wrap(~biotype)+
  labs(title = "ABC score per target biotype and enhancer location",
       x = "Location", y = "ABC score") +
  scale_fill_manual(values = nejm_pal) + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

# Number of targets per enhancer ----
Targets_per_enhancer <- non_self_enhancers%>% group_by(name,EnhancerBiotype,EnhancerGene) %>%
  summarise(N_targets=length(unique(TargetGene)),
            N_PCG_targets=length(unique(TargetGene[biotype=="protein_coding"])),
          maxABC_score=max(ABC.Score))

ggplot(Targets_per_enhancer, aes(x = N_targets)) +
  geom_histogram(binwidth = 1, fill = nejm_pal[1], color = "black") +
  facet_wrap(~EnhancerBiotype)+
  labs(title = "Number of Targets per Enhancer biotype",
       x = "Number of Targets", y = "Frequency")

ggplot(Targets_per_enhancer%>%filter(!EnhancerBiotype%in%c("intergenic","protein_coding")), aes(x = N_targets)) +
  geom_histogram(binwidth = 1, fill = nejm_pal[1], color = "black") +
  facet_wrap(~EnhancerBiotype)+
  labs(title = "Number of Targets per Enhancer biotype",
       x = "Number of Targets", y = "Frequency")

Targets_per_enhancer%>% group_by(EnhancerBiotype) %>% summarise(mean(N_targets),
                                                                median(N_targets),
                                                                n(),
                                                                sum(N_targets>1),
                                                                sum(N_PCG_targets>1))


ggplot(Targets_per_enhancer,aes(x=N_targets,y=maxABC_score)) + geom_point()

ggplot(Targets_per_enhancer,aes(x=N_targets,y=maxABC_score)) + geom_point() +
  facet_wrap(~EnhancerBiotype)

ggplot(Targets_per_enhancer,aes(x=N_PCG_targets,y=maxABC_score)) + geom_point() +
  facet_wrap(~EnhancerBiotype)

# Select those potential novel harboring enhancers that ----
# have a score higher than 0.1 and regulate 5 or more PCGs

potNovel_confident_enhancers <- Targets_per_enhancer%>%filter(EnhancerBiotype=="potNovel",
                                                              N_PCG_targets>4,
                                                              maxABC_score>0.1)
table(gene_level_info$exonic_type[gene_level_info$gene_name%in%potNovel_confident_enhancers$EnhancerGene])

genes_regulated_by_potNovel_confident_enhancers <- unique(non_self_enhancers$TargetGene[non_self_enhancers$name%in%potNovel_confident_enhancers$name])

go_enrichment <- enrichGO(gene = genes_regulated_by_potNovel_confident_enhancers,
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)

# View the top GO terms
head(go_enrichment)

# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) +
  ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")

# Targets per EnhancerGene ----
non_self_enhancers_long <- separate_longer_delim(non_self_enhancers,
                                                 col = "EnhancerGene",
                                                 delim = ",")
Targets_per_enhancerGene <- non_self_enhancers_long %>% filter(!is.na(EnhancerGene)) %>%
  group_by(EnhancerGene) %>% summarise(N_targets_asEnh=length(unique(TargetGene)),
                                       N_PCGtargets_asEnh=length(unique(TargetGene[biotype=="protein_coding"])),
                                       maxABC_score_asEnh=max(ABC.Score))

# modify gene_level_info ----
gene_level_info <- gene_level_info %>% mutate(harborsEnhancer =
                                                gene_name%in%Targets_per_enhancerGene$EnhancerGene)

gene_level_info <- left_join(gene_level_info, Targets_per_enhancerGene, by=c(gene_name="EnhancerGene"))


# Write enhancer data to database ----
library(DBI)
library(RSQLite)
db <- dbConnect(RSQLite::SQLite(), "outputs/dbs/StemLinc.db")

gene_ABC_predictions_enhancers_and_or_targets <- gene_level_info[,c(1,57:65)]
gene_ABC_predictions_enhancers_and_or_targets <- gene_ABC_predictions_enhancers_and_or_targets %>%
  filter(Enhancer_target|harborsEnhancer)
head(gene_ABC_predictions_enhancers_and_or_targets)

dbExecute(db, "
    CREATE TABLE gene_ABC_predictions_enhancers_and_or_targets (
        gene_name TEXT PRIMARY KEY,
        Enhancer_target BOOLEAN,
        N_enh_asTarget INTEGER,
        maxABCscore_asTarget REAL,
        maxEnh_asTarget TEXT,
        maxEnh_biotype_asTarget TEXT,
        harborsEnhancer BOOLEAN,
        N_targets_asEnh INTEGER,
        N_PCGtargets_asEnh INTEGER,
        maxABC_score_asEnh REAL,
        FOREIGN KEY (gene_name) REFERENCES gene_core_data(gene_name)
    );
")

dbWriteTable(
  conn = db,
  name = "gene_ABC_predictions_enhancers_and_or_targets",
  value = gene_ABC_predictions_enhancers_and_or_targets,
  append = TRUE,
  row.names = FALSE
)

metadata <- data.frame(
  table_name = "gene_ABC_predictions_enhancers_and_or_targets",
  script_name = "analyse_ABC_results.R",
  level = "gene",
  filter_status = " 3 replicates gene level, >= 0.2 TPM, only biotypes of interest, only genes that harbor ABC-predicted enhancers and/or are targets of them",
  analysis_version = "v1.0",
  timestamp = Sys.time(),
  source_data = "ABC-Enhancer-Gene-Prediction/results/LT.HSC_with_HiC.R1/Predictions/parameters.predict.txt",
  orig_scripts = "scripts on ABC-Enhancer-Gene-Prediction were run on the IJC hpc",  # Replace with your actual script
  description = "Summary of genes that harbor ABC-predicted enhancers and/or are targets of them, including number of enhancers and targets and max ABC scores."
)

# Append the metadata to the database
dbWriteTable(db, "metadata", metadata, append = TRUE, row.names = FALSE)

# lncRNAs that are enhancers of PCGs
non_self_enhancers %>% filter(biotype=="protein_coding") %>%
  group_by(EnhancerBiotype) %>% summarise(N_genes=length(unique(EnhancerGene)),
                                          N_enhancers=n())

# All enhancers in common ----
N_targets_Enhancers <- non_self_enhancers %>%
  group_by(name,biotype) %>%
  summarise(count=n())

N_targets_Enhancers <- pivot_wider(N_targets_Enhancers,
                                   names_from = biotype,
                                   values_from = count,
                                   id_cols = name)

N_targets_Enhancers[is.na(N_targets_Enhancers)] <- 0

N_targets_Enhancers$TotCount <- rowSums(N_targets_Enhancers[,2:ncol(N_targets_Enhancers)])

N_targets_Enhancers$ncCount <- rowSums(N_targets_Enhancers[,c(2,5,6)])


N_targets_Enhancers <- N_targets_Enhancers %>%
  mutate(Enhancer_targets=
           ifelse(TotCount==protein_coding,
                  "PCG_only",
                  ifelse(protein_coding!=0&ncCount!=0,
                         "PCG_and_nc",
                         ifelse(protein_coding!=0&pseudogene!=0,
                                "PCG_and_pseudo","nc_only" ))))

ggplot(N_targets_Enhancers,aes(x=Enhancer_targets)) + geom_bar() +
  labs(title="All enhancers", y="Enhancer count",
       x="Enhancer targets")

# Genes regulated by all PCG_and_nc enhancers ----
PCG_and_nc_enhancers <- N_targets_Enhancers$name[N_targets_Enhancers$Enhancer_targets=="PCG_and_nc"]

genes_regulated_Enh_common_PCG_and_nc <- non_self_enhancers %>% filter(name%in%PCG_and_nc_enhancers) %>%
  group_by(biotype) %>% summarise(N_targets=length(unique(TargetGene)))


# GO analysis of PCG that have enhancers in common with ncRNAs ----
PCGs_with_enhancers_common2lncRNAs <-
  non_self_enhancers$TargetGene[non_self_enhancers$name%in%PCG_and_nc_enhancers&
                                  non_self_enhancers$biotype=="protein_coding"]

PCGs_with_enhancers_common2lncRNAs=unique(PCGs_with_enhancers_common2lncRNAs)


go_enrichment <- enrichGO(gene = PCGs_with_enhancers_common2lncRNAs,
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)
# View the top GO terms
head(go_enrichment)

# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) +
  ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")


# Intergenic enhancers in common ----
N_targets_intergenic_Enhancers <- non_self_enhancers %>%
  filter(class=="intergenic") %>% group_by(name,biotype) %>%
  summarise(count=n())

N_targets_intergenic_Enhancers <- pivot_wider(N_targets_intergenic_Enhancers,
                                              names_from = biotype,
                                              values_from = count,
                                              id_cols = name)

N_targets_intergenic_Enhancers[is.na(N_targets_intergenic_Enhancers)] <- 0

N_targets_intergenic_Enhancers$TotCount <- rowSums(N_targets_intergenic_Enhancers[,2:ncol(N_targets_intergenic_Enhancers)])

N_targets_intergenic_Enhancers$ncCount <- rowSums(N_targets_intergenic_Enhancers[,c(2,5,6)])


N_targets_intergenic_Enhancers <- N_targets_intergenic_Enhancers %>%
  mutate(Enhancer_targets=ifelse(TotCount==protein_coding,
                                 "PCG_only",
                                 ifelse(protein_coding!=0&ncCount!=0,
                                        "PCG_and_nc",
                                        ifelse(protein_coding!=0&pseudogene!=0,
                                               "PCG_and_pseudo","nc_only" ))))

ggplot(N_targets_intergenic_Enhancers,aes(x=Enhancer_targets)) + geom_bar() +
  labs(title="Intergenic Enhancers", ylab="Enhancer count",
       xlab="Enhancer targets")

# Genes regulated by intergenic PCG_and_nc enhancers
PCG_and_nc_intergenic_enhancers <- N_targets_intergenic_Enhancers$name[N_targets_intergenic_Enhancers$Enhancer_targets=="PCG_and_nc"]

genes_regulated_IntEnh_common_PCG_and_nc <- non_self_enhancers %>% filter(name%in%PCG_and_nc_intergenic_enhancers) %>%
  group_by(biotype) %>% summarise(N_targets=length(unique(TargetGene)))

# GO analysis of PCG that have intergenic enhancers in common with lncRNAs ----

PCGs_with_intergenic_enhancers_common2lncRNAs <-
  non_self_enhancers$TargetGene[non_self_enhancers$name%in%PCG_and_nc_intergenic_enhancers&
                                  non_self_enhancers$biotype=="protein_coding"]

PCGs_with_intergenic_enhancers_common2lncRNAs=unique(PCGs_with_intergenic_enhancers_common2lncRNAs)



go_enrichment <- enrichGO(gene = PCGs_with_intergenic_enhancers_common2lncRNAs,
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)

# View the top GO terms
head(go_enrichment)

# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) +
  ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")


# Genic enhancers in common ----
N_targets_genic_Enhancers <- non_self_enhancers %>%
  filter(class=="genic") %>% group_by(name,biotype) %>%
  summarise(count=n())

N_targets_genic_Enhancers <- pivot_wider(N_targets_genic_Enhancers,
                                         names_from = biotype,
                                         values_from = count,
                                         id_cols = name)

N_targets_genic_Enhancers[is.na(N_targets_genic_Enhancers)] <- 0

N_targets_genic_Enhancers$TotCount <- rowSums(N_targets_genic_Enhancers[,2:ncol(N_targets_genic_Enhancers)])

N_targets_genic_Enhancers$ncCount <- rowSums(N_targets_genic_Enhancers[,c(2,5,6)])


N_targets_genic_Enhancers <- N_targets_genic_Enhancers %>% mutate(Enhancer_targets=ifelse(TotCount==protein_coding,
                                                                                          "PCG_only",
                                                                                          ifelse(protein_coding!=0&ncCount!=0,
                                                                                                 "PCG_and_nc",
                                                                                                 ifelse(protein_coding!=0&pseudogene!=0,
                                                                                                        "PCG_and_pseudo","nc_only" ))))

ggplot(N_targets_genic_Enhancers,aes(x=Enhancer_targets)) + geom_bar() +
  labs(title="genic Enhancers", ylab="Enhancer count",
       xlab="Enhancer targets")

# Genes regulated by genic PCG_and_nc enhancers
PCG_and_nc_genic_enhancers <- N_targets_genic_Enhancers$name[N_targets_genic_Enhancers$Enhancer_targets=="PCG_and_nc"]

genes_regulated_GenEnh_common_PCG_and_nc <- non_self_enhancers %>% filter(name%in%PCG_and_nc_genic_enhancers) %>%
  group_by(biotype) %>% summarise(N_targets=length(unique(TargetGene)))




# GO analysis of PCG that have genic enhancers in common with lncRNAs ----
PCGs_with_genic_enhancers_common2lncRNAs <-
  non_self_enhancers$TargetGene[non_self_enhancers$name%in%PCG_and_nc_genic_enhancers&
                                  non_self_enhancers$biotype=="protein_coding"]

PCGs_with_genic_enhancers_common2lncRNAs=unique(PCGs_with_genic_enhancers_common2lncRNAs)


go_enrichment <- enrichGO(gene = PCGs_with_genic_enhancers_common2lncRNAs,
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)

# View the top GO terms
head(go_enrichment)

# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) +
  ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")

# GO analysis of PCGs only regulated by enhancers exclusive to PCGs ----
PCG_only_enhancers <- N_targets_Enhancers$name[N_targets_Enhancers$Enhancer_targets=="PCG_only"]

PCGs_with_PCG_only_enhancers <-
  non_self_enhancers$TargetGene[non_self_enhancers$name%in%PCG_only_enhancers&
                                  non_self_enhancers$biotype=="protein_coding"]

PCGs_with_PCG_only_enhancers=unique(PCGs_with_PCG_only_enhancers)

PCGs_with_PCG_only_enhancers <- PCGs_with_PCG_only_enhancers[!PCGs_with_PCG_only_enhancers%in%PCGs_with_enhancers_common2lncRNAs]

go_enrichment <- enrichGO(gene = PCGs_with_PCG_only_enhancers,
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)



# View the top GO terms
head(go_enrichment)

# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) +
  ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")

# GO analysis of PCG harboring enhancers regulating potNovel ----
PCG_with_enhancers_reg_potNovel <-
  non_self_enhancers$EnhancerGene[non_self_enhancers$biotype=="potNovel"&
                                                                     non_self_enhancers$EnhancerBiotype=="protein_coding"]

PCG_with_enhancers_reg_potNovel=unique(unlist(strsplit(PCG_with_enhancers_reg_potNovel,split = ",")))

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

go_enrichment <- enrichGO(gene = PCG_with_enhancers_reg_potNovel,
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)

# View the top GO terms
head(go_enrichment)

# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")


# GO analyses of PCGs regulated by enhancers in potNovel ----
PCGs_targets_of_potNovel_enhancers=unique(non_self_enhancers$TargetGene[non_self_enhancers$biotype=="protein_coding"&
                                                                   non_self_enhancers$EnhancerBiotype=="potNovel"])

go_enrichment <- enrichGO(gene = PCGs_targets_of_potNovel_enhancers,
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)

# View the top GO terms
head(go_enrichment)

# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")


# ABC score for potNovel genes ----

## as enhancers ----
potNovel_enhancers <- non_self_enhancers %>% filter(EnhancerBiotype=="potNovel")

potNovel_enhancers <- separate_longer_delim(potNovel_enhancers,cols = "EnhancerGene",delim = ",")
maxABC_score_potNovel_enhancers <- potNovel_enhancers %>% group_by(EnhancerGene) %>%
  summarise(maxABC_score=max(ABC.Score),
            N_targets=length(unique(TargetGene)))

maxABC_score_potNovel_enhancers <- maxABC_score_potNovel_enhancers%>%filter(grepl("XLOC",maxABC_score_potNovel_enhancers$EnhancerGene))

maxABC_score_potNovel_enhancers <- left_join(maxABC_score_potNovel_enhancers,
                                             gene_level_info%>%dplyr::select(gene_name,
                                                                             max_corr_signature,
                                                                             mean_corr_signature, exonic_type,
                                                                             best_classif_to_PCG,best_closest_PCG),
                                             by=c(EnhancerGene="gene_name"))


### correlation corr signature ABC score ----
ggplot(maxABC_score_potNovel_enhancers,aes(maxABC_score,max_corr_signature)) + geom_point()
cor(maxABC_score_potNovel_enhancers$maxABC_score,maxABC_score_potNovel_enhancers$max_corr_signature)

ggplot(maxABC_score_potNovel_enhancers,aes(maxABC_score,mean_corr_signature)) + geom_point()
cor(maxABC_score_potNovel_enhancers$maxABC_score,maxABC_score_potNovel_enhancers$mean_corr_signature)

mean(maxABC_score_potNovel_enhancers$maxABC_score)

## as targets ----
potNovel_targets <- non_self_enhancers %>% filter(biotype=="potNovel")

maxABC_score_potNovel_targets <- potNovel_targets %>% group_by(TargetGene) %>%
  summarise(maxABC_score=max(ABC.Score))

maxABC_score_potNovel_targets <- left_join(maxABC_score_potNovel_targets,
                                             gene_level_info%>%dplyr::select(gene_name,
                                                                             max_corr_signature,
                                                                             mean_corr_signature),
                                             by=c(TargetGene="gene_name"))


ggplot(maxABC_score_potNovel_targets,aes(maxABC_score,max_corr_signature)) + geom_point()
cor(maxABC_score_potNovel_targets$maxABC_score,maxABC_score_potNovel_targets$max_corr_signature)

ggplot(maxABC_score_potNovel_targets,aes(maxABC_score,mean_corr_signature)) + geom_point()
cor(maxABC_score_potNovel_targets$maxABC_score,maxABC_score_potNovel_targets$mean_corr_signature)

mean(maxABC_score_potNovel_targets$maxABC_score)

# How are enhancers of potNovel genes ----
potNovel_targets %>% group_by(EnhancerBiotype) %>% summarise(mean_dist=mean(distance),
                                                             median_dist=median(distance),
                                                             mean_ABC_score=mean(ABC.Score),
                                                             median_ABC_score=median(ABC.Score))
ggplot(potNovel_targets,aes(distance,ABC.Score)) + geom_point() +
  scale_x_log10()
cor(potNovel_targets$distance, potNovel_targets$ABC.Score)

# GO analyses of PCGs regulated by enhancers in annotated lncRNAS ----
PCGs_targets_of_lncRNA_enhancers=unique(non_self_enhancers$TargetGene[non_self_enhancers$biotype=="protein_coding"&
                                                                          non_self_enhancers$EnhancerBiotype=="lncRNA"])

go_enrichment <- enrichGO(gene = PCGs_targets_of_lncRNA_enhancers,
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)

# View the top GO terms
head(go_enrichment)

# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")


# Add TAD info ----
TADs_Chen <- read.table("data/public_data/HiC_data/GSE119347_BMHSC_TADs.mm39.bed")

TADs_Chen$TAD_ID=paste0(TADs_Chen$V1,":",TADs_Chen$V2,"-",TADs_Chen$V3)

colnames(TADs_Chen) <- c("TAD_chr","TAD_start","TAD_end","TAD_ID")

non_self_enhancers_coords <- non_self_enhancers[,1:4]
non_self_enhancers_coords <- non_self_enhancers_coords[!duplicated(non_self_enhancers_coords$name),]

enhancers_in_TADs <- bt.intersect(a = non_self_enhancers_coords,
                                  b=TADs_Chen,wo = T)

colnames(enhancers_in_TADs) <- c(colnames(non_self_enhancers_coords),
                                 colnames(TADs_Chen),"ol_length")

genes_in_TADs <- bt.intersect(a = genes,
                              b = TADs_Chen,
                              wo = T)

colnames(genes_in_TADs) <- c(colnames(genes),colnames(TADs_Chen),"ol_length")

enhancers_in_TADs <- enhancers_in_TADs %>% group_by(name) %>%
  summarise(TADs=paste0(TAD_ID,collapse = ","))

genes_in_TADs <- genes_in_TADs %>% group_by(gene_name) %>%
  summarise(TADs=paste0(TAD_ID,collapse = ","))

non_self_enhancers$Enhancer_TAD <- enhancers_in_TADs$TADs[match(non_self_enhancers$name,
                                                                enhancers_in_TADs$name)]

non_self_enhancers$Gene_TAD <- genes_in_TADs$TADs[match(non_self_enhancers$TargetGene,
                                                                genes_in_TADs$gene_name)]


# Filter for hematopoietic genes---------------------------
# Only keep predictions for genes in the hematopoietic list
non_self_enhancers_hemat <- non_self_enhancers %>%
  filter(TargetGene %in% hematopoietic_genes$TargetGene)


# Number of Enhancers per Hematopoiesis Gene---------------------------
# Calculate number of unique enhancers per gene
enhancer_counts <-  non_self_enhancers_hemat %>%
  group_by(TargetGene) %>%
  summarise(n_enhancers = n()) %>%
  ungroup()

# Plot: Number of Enhancers per Gene
ggplot(enhancer_counts, aes(x = n_enhancers)) +
  geom_histogram(binwidth = 1, fill = nejm_pal[1], color = "black") +
  labs(title = "Number of Enhancers per HSC-related Gene",
       x = "Number of Enhancers", y = "Frequency")

# Enhancer Location Relative to Genes---------------------------

# Plot: Enhancer Location Distribution
ggplot(non_self_enhancers_hemat, aes(x = class)) +
  geom_bar(fill = nejm_pal[2]) +
  labs(title = "Location of Enhancers Relative to Genes",
       x = "Location", y = "Count")

non_self_enhancers_hemat_genic <- non_self_enhancers_hemat %>% filter(class=="genic")
non_self_enhancers_hemat_genic <- non_self_enhancers_hemat_genic %>%
  mutate(EnhancerBiotype=ifelse(grepl("protein_coding",EnhancerBiotype),
                                "protein_coding",ifelse(grepl("lncRNA",EnhancerBiotype),
                                                        "lncRNA",EnhancerBiotype)))

ggplot(non_self_enhancers_hemat_genic, aes(x = EnhancerBiotype)) +
  geom_bar(fill = nejm_pal[2]) +
  labs(title = "Biotypes of genic enhancers",
       x = "Biotype", y = "Count")

non_self_enhancers_hemat_genic %>% group_by(EnhancerBiotype) %>%
  summarise(N_enhancers=length(unique(name)))

# ABC Score Distribution for Hematopoiesis Genes---------------------------
# Plot the distribution of ABC scores
ggplot(non_self_enhancers_hemat, aes(x = ABC.Score)) +
  geom_density(fill = nejm_pal[4]) +
  labs(title = "ABC Score Distribution for HSC-related Genes",
       x = "ABC Score", y = "Density")


# enhancers of 79903

xloc="XLOC_079903"

xloc_enhancers <- non_self_enhancers%>%filter(TargetGene==xloc) %>% pull(name)

genes_coregulated_xloc <- non_self_enhancers%>%filter(name%in%xloc_enhancers) %>% pull(TargetGene)
gene_level_info[gene_level_info$gene_name%in%genes_coregulated_xloc,c("gene_name","mean_corr_signature")]

# # write non-self enhancers bed ----
# non_self_enhancers_bed <- non_self_enhancers %>%filter(!duplicated(enh_idx)) %>% dplyr::select(Enh_chr,Enh_start,Enh_end,enh_idx)
# # write.table(non_self_enhancers_bed,
# #             "outputs/bed_files/non_self_enhancers_HSC_all_data_powerlaw.bed",quote = F,col.names = F,row.names = F,sep = "\t")
# write.table(non_self_enhancers_bed,
#             "outputs/bed_files/non_self_enhancers_HSC_all_data_HiCR1.bed",quote = F,col.names = F,row.names = F,sep = "\t")


# Enhancers of hematopoietic genes ----
# select enhancers that regulate the previously defined 440 hematopoiesis related genes ----
enhancers_of_hematPCGs <- unique(
  non_self_enhancers %>%
    filter(TargetGene%in%hematopoietic_genes$TargetGene) %>%
    pull(name))

length(enhancers_of_hematPCGs)

length(unique(non_self_enhancers$name))
length(enhancers_of_hematPCGs)/length(unique(non_self_enhancers$name))*100

# check biotype of these enhancers ----
Enhancer_annot <- non_self_enhancers %>% dplyr::select(name, EnhancerBiotype,EnhancerGene)
Enhancer_annot <- Enhancer_annot[!duplicated(Enhancer_annot),]

table(Enhancer_annot$EnhancerBiotype[Enhancer_annot$name%in%enhancers_of_hematPCGs])

# how many of these enhancers also regulate other genes? ----
genes_regulated_by_hemat_enhancers <- non_self_enhancers%>%
  filter(name%in%enhancers_of_hematPCGs)

genes_regulated_by_hemat_enhancers_N_targets <- genes_regulated_by_hemat_enhancers %>%
  group_by(name,EnhancerBiotype) %>% summarise(N_targets=length(unique(TargetGene)))

ggplot(genes_regulated_by_hemat_enhancers_N_targets,aes(x=N_targets)) + geom_histogram() +
  facet_wrap(~EnhancerBiotype)

other_genes_regulated_by_hemat_enhancers <- genes_regulated_by_hemat_enhancers %>% filter(!TargetGene%in%hematopoietic_genes$TargetGene)
length(unique(other_genes_regulated_by_hemat_enhancers$TargetGene))
# these enhancers also regulate other 7975 genes
# how are these genes
table(other_genes_regulated_by_hemat_enhancers$biotype[!duplicated(other_genes_regulated_by_hemat_enhancers$TargetGene)])

other_genes_regulated_by_hemat_enhancers_info <- gene_level_info %>% filter(gene_name%in%other_genes_regulated_by_hemat_enhancers$TargetGene)

table(other_genes_regulated_by_hemat_enhancers_info$biotype,
      other_genes_regulated_by_hemat_enhancers_info$exonic_type)
# Do all enhancers that regulate genes related to hematopoiesis regulate also other genes?
table(enhancers_of_hematPCGs%in%other_genes_regulated_by_hemat_enhancers$name)

# for each gene generate a list of hemat genes co-regulated by enhancers ----

# Step 2: Join with Original Data to Identify Common Enhancers
data <- genes_regulated_by_hemat_enhancers%>%dplyr::select(name,TargetGene)
hemat_related <- data %>%
  filter(TargetGene %in% hematopoietic_genes$TargetGene)

Genes_with_coregulated_hemat_genes <- data %>%
  inner_join(hemat_related, by = "name", suffix = c("", "_hemat")) %>%
  filter(TargetGene != TargetGene_hemat) %>%
  group_by(TargetGene) %>%
  summarise(
    Hemat_Targets_shared_enhancer = paste(unique(TargetGene_hemat), collapse = ", "),
    N_Hemat_Targets_shared_enhancer = length(unique(TargetGene_hemat)),
    N_Hemat_enhancers = length(unique(name))
  )

Genes_with_coregulated_hemat_genes <- left_join(Genes_with_coregulated_hemat_genes,
                                                gene_level_info %>% dplyr::select(gene_name,
                                                                                  biotype,
                                                                                  best_classif_to_PCG,
                                                                                  best_closest_PCG,
                                                                                  dist_PCG,
                                                                                  width,
                                                                                  CAGE_within_100bp,
                                                                                  polyAsite_within_100bp,
                                                                                  H3K4me3_at_promoter,
                                                                                  H3K4me1_at_promoter,
                                                                                  H3K27ac_at_promoter,
                                                                                  H3K36me3_at_geneBody,
                                                                                  Enhancer_Atlas_at_promoter,
                                                                                  FANTOM_enhancers_at_promoter,
                                                                                  H3K27me3_at_promoter,
                                                                                  N_TFs,
                                                                                  TF_rel.lev_A,

                                                                                  harborsEnhancer,
                                                                                  N_enhTargets,
                                                                                  enh_maxABC_score
                                                                                  ), by=c(TargetGene="gene_name"))


#write.csv(Genes_with_coregulated_hemat_genes,"outputs/ABC_files/Genes_with_coregulated_hemat_genes.txt",row.names = F)

# select all the hemat genes coregulated potNovel and perform GO analysis ----
hemat_genes_common_potNovel <- Genes_with_coregulated_hemat_genes %>%
  filter(biotype=="potNovel") %>% pull(Hemat_Targets_shared_enhancer)
hemat_genes_common_potNovel=unlist(strsplit(hemat_genes_common_potNovel,split = ", "))
length(hemat_genes_common_potNovel)

# Count occurrences to use as weights
gene_weights <- table(hemat_genes_common_potNovel)
gene_weights <- as.numeric(gene_weights)

# Create a named vector of weights
names(gene_weights) <- names(table(hemat_genes_common_potNovel))

go_enrichment <- enrichGO(gene = hemat_genes_common_potNovel,
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",  # Use "ENSEMBL" if your IDs are Ensembl
                          ont = "BP",          # "BP" for Biological Process
                          pAdjustMethod = "BH", # Benjamini-Hochberg adjustment
                          pvalueCutoff = 0.05,  # p-value cutoff
                          qvalueCutoff = 0.05)
# View the top GO terms
head(go_enrichment)

# Plot results with a bar plot
barplot(go_enrichment, showCategory = 10) +
  ggtitle("Top 10 Enriched GO Terms")

# Alternatively, dot plot
dotplot(go_enrichment, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")

View(as.data.frame(go_enrichment))

