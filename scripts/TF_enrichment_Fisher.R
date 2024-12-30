#load data ----
source("scripts/load_gene_level_data_biotypes_of_interest.R")




separate_data <- function(data,col,split=","){
  data <- as.data.frame(data)
  col2sep <- data[,col]
  seplist <- strsplit(col2sep,split = split)
  lens <- sapply(seplist,length)
  long_data <- data[rep(1:nrow(data),lens),]
  if(is.numeric(col)) col=colnames(data)[col]
  new_col=paste0(col,"_split")
  long_data[,new_col] <- unlist(seplist)
  return(long_data)
}
#test_data=TF_data[1:10,1:2]
#test_out=separate_data(test_data,"TF_rel.lev_A")
# Function to perform Fisher's Exact Test for a given TF and biotype

perform_fisher_test <- function(tf, biot, gene_biotypes, tf_promoters) {

  # Total number of genes
  total_genes <- nrow(gene_biotypes)

  # Genes in the biotype of interest
  biotype_genes <- gene_biotypes %>% dplyr::filter(biotype == biot) %>% dplyr::pull(gene_name)
  num_biotype <- length(biotype_genes)

  # Genes not in the biotype of interest
  other_genes <- setdiff(gene_biotypes$gene_name, biotype_genes)
  num_other <- total_genes - num_biotype

  # Genes in biotype with the TF
  a <- tf_promoters %>%
    filter(TF == tf, gene_name %in% biotype_genes) %>%
    distinct(gene_name) %>%
    nrow()

  # Genes in biotype without the TF
  b <- num_biotype - a

  # Genes not in biotype with the TF
  c <- tf_promoters %>%
    filter(TF == tf, gene_name %in% other_genes) %>%
    distinct(gene_name) %>%
    nrow()

  # Genes not in biotype without the TF
  d <- num_other - c

  # Create contingency table
  contingency_table <- matrix(c(a, b, c, d), nrow = 2,
                              dimnames = list(
                                Biotype = c(biot, "Other"),
                                TF_Present = c("Yes", "No")
                              ))

  # Perform Fisher's Exact Test
  test <- fisher.test(contingency_table, alternative = "greater")  # 'greater' tests for enrichment

  # Return results
  return(data.frame(
    TF = tf,
    Biotype = biot,
    a = a,
    b = b,
    c = c,
    d = d,
    p_value = test$p.value,
    odds_ratio = test$estimate
  ))
}


# preprocess data ----
gene_biotypes=gene_level_info%>%select(gene_name,biotype)
table(gene_biotypes$biotype)

TF_data=read.table("outputs/overlap_marks/LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif.TFs_rel_level_500_250.tsv",header = T)
table(TF_data$gene_name%in%gene_biotypes$gene_name)

TF_data <- TF_data%>%filter(gene_name%in%gene_biotypes$gene_name)

# For now use only rel level A:
TF_data_long <- TF_data %>% select(gene_name,TF_rel.lev_A)
TF_data_long <- separate_data(TF_data_long,2)
TF_data_long <- TF_data_long[,c(1,3)]
colnames(TF_data_long)[2]="TF"
TF_data_long <- TF_data_long %>% filter(!duplicated(TF_data_long))
TF_data_long <- TF_data_long %>% filter(!is.na(TF))

write.table(TF_data_long,"outputs/TF_enrichment/TF_data_long.txt",quote = F,row.names = F,sep = "\t")
# run ----
library(purrr)

## all biotypes ----
# Get list of unique TFs and biotypes
unique_tfs <- unique(TF_data_long$TF)
unique_biotypes <- unique(gene_biotypes$biotype)

# Create all combinations of TFs and biotypes
combinations <- expand.grid(tf = unique_tfs, biotype = unique_biotypes, stringsAsFactors = FALSE)

# Apply the function to each combination
results <- combinations %>%
  pmap_dfr(~ perform_fisher_test(..1, ..2, gene_biotypes, TF_data_long))

# Adjust p-values using Benjamini-Hochberg method
results <- results %>%
  group_by(Biotype) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
  ungroup()


# Define significance threshold
significance_threshold <- 0.05

# Filter for significant results
significant_results <- results %>%
  filter(adj_p_value < significance_threshold)


write.table(results,"outputs/TF_enrichment/results_TF_rel.levelA.fisher.txt",
            quote = F,row.names = F,sep = "\t")
write.table(significant_results,"outputs/TF_enrichment/results_TF_rel.levelA.fisher.significant.txt",
            quote = F,row.names = F,sep = "\t")

# exonic type ----
## all genes ----
gene_exonic_type <- gene_level_info %>% select(gene_name,exonic_type)

colnames(gene_exonic_type)[2] <- "biotype"


# Get list of unique TFs and biotypes
unique_tfs <- unique(TF_data_long$TF)
unique_biotypes <- unique(gene_exonic_type$biotype)

# Create all combinations of TFs and biotypes
combinations <- expand.grid(tf = unique_tfs, biotype = unique_biotypes, stringsAsFactors = FALSE)

# Apply the function to each combination
results <- combinations %>%
  pmap_dfr(~ perform_fisher_test(..1, ..2, gene_exonic_type, TF_data_long))

# Adjust p-values using Benjamini-Hochberg method
results <- results %>%
  group_by(Biotype) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
  ungroup()


# Define significance threshold
significance_threshold <- 0.05

# Filter for significant results
significant_results <- results %>%
  filter(adj_p_value < significance_threshold)


write.table(results,"outputs/TF_enrichment/results_TF_rel.levelA.fisher.exonic_type_allGenes.txt",
            quote = F,row.names = F,sep = "\t")
write.table(significant_results,"outputs/TF_enrichment/results_TF_rel.levelA.fisher.significant.exonic_type_allGenes.txt",
            quote = F,row.names = F,sep = "\t")


## per biotype ----

for (biot in unique(gene_level_info$biotype)) {

  gene_table <- gene_level_info %>% filter(biotype==biot) %>% select(gene_name, exonic_type)
  colnames(gene_table)[2] <- "biotype"


  # Get list of unique TFs and biotypes
  unique_tfs <- unique(TF_data_long$TF)
  unique_biotypes <- unique(gene_table$biotype)

  # Create all combinations of TFs and biotypes
  combinations <- expand.grid(tf = unique_tfs, biotype = unique_biotypes, stringsAsFactors = FALSE)

  # Apply the function to each combination
  results <- combinations %>%
    pmap_dfr(~ perform_fisher_test(..1, ..2, gene_table, TF_data_long))

  # Adjust p-values using Benjamini-Hochberg method
  results <- results %>%
    group_by(Biotype) %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
    ungroup()


  # Define significance threshold
  significance_threshold <- 0.05

  # Filter for significant results
  significant_results <- results %>%
    filter(adj_p_value < significance_threshold)


  write.table(results,paste0("outputs/TF_enrichment/results_TF_rel.levelA.fisher.exonic_type_",biot,".txt"),
              quote = F,row.names = F,sep = "\t")
  write.table(significant_results,paste0("outputs/TF_enrichment/results_TF_rel.levelA.fisher.significant.exonic_type_",biot,".txt"),
              quote = F,row.names = F,sep = "\t")


}

# High corr signature vs low corr ----
mean_corr_co=0.65
gene_level_info <- read.csv("outputs/gene_level_info_all_evidence_oct24.csv")
TF_data_long <- left_join(TF_data_long,gene_level_info%>%dplyr::select(gene_name,
                                                                       biotype,
                                                                       mean_corr_signature))
TF_data_nc <- TF_data_long %>% filter(biotype!="protein_coding")
TF_data_nc <- TF_data_nc %>% mutate(corr_class=ifelse(mean_corr_signature>mean_corr_co,"highCorr",
                                    ifelse(mean_corr_signature<0,"lowCorr","medianCorr")))



TF_data_nc$biotype=TF_data_nc$corr_class

TF_data_nc <- TF_data_nc %>% filter(biotype!="medianCorr")
TF_data <- TF_data_nc %>% select(gene_name,TF)

gene_data <- TF_data_nc %>% select(gene_name,biotype)

# Get list of unique TFs and biotypes
unique_tfs <- unique(TF_data$TF)
unique_biotypes <- unique(gene_data$biotype)

# Create all combinations of TFs and biotypes
combinations <- expand.grid(tf = unique_tfs,
                            biotype = unique_biotypes,
                            stringsAsFactors = FALSE)

# Apply the function to each combination
results <- combinations %>%
  pmap_dfr(~ perform_fisher_test(..1, ..2, gene_data, TF_data))

# Adjust p-values using Benjamini-Hochberg method
results <- results %>%
  group_by(Biotype) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
  ungroup()


# Define significance threshold
significance_threshold <- 0.05

# Filter for significant results
significant_results <- results %>%
  filter(adj_p_value < significance_threshold)


write.table(results,"outputs/TF_enrichment/results_TF_rel.levelA.fisher.corrClass_ncGenes.txt",
            quote = F,row.names = F,sep = "\t")
write.table(significant_results,"outputs/TF_enrichment/results_TF_rel.levelA.fisher.significant.corrClass_ncGenes.txt",
            quote = F,row.names = F,sep = "\t")

sort(table(TF_data_nc$TF[TF_data_nc$biotype=="highCorr"]))
