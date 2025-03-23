# Prepare the data
log_TPM_data <- log2(TPM_data + 0.001)
TPM_data_df <- as.data.frame(log_TPM_data)
TPM_data_df$gene_id <- rownames(TPM_data)

# Combine with biotypes
table(TPM_data_df$gene_id==gene_level_info$gene_name)
biotype_df <- data.frame(gene_id = rownames(TPM_data),
                         biotype = gene_level_info$biotype)
merged_data <- merge(TPM_data_df, biotype_df, by = "gene_id")

set.seed(123)  # Set seed for reproducibility
sampled_genes <- merged_data %>%
  group_by(biotype) %>%
  sample_n(1000, replace = FALSE)
# Reshape the data to long format
long_data <- melt(merged_data, id.vars = c("gene_id", "biotype"))
# Reshape the data to long format
long_data <- melt(sampled_genes, id.vars = c("gene_id", "biotype"))

# Create heatmap using ggplot
g=ggplot(long_data, aes(x = variable, y = gene_id, fill = value)) +
  geom_tile() +
  facet_wrap(~ biotype, scales = "free_y") +  # Separate heatmaps for each biotype
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = median(long_data$value, na.rm = TRUE)) +
  theme_minimal() +
  labs(x = "Samples", y = "Genes", fill = "TPM") +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5 ,
                                   size =2)) +
  theme(axis.text.y = element_blank(),  # Remove y-axis text
        axis.ticks.y = element_blank())
g
outdir="outputs/plots"

save_tiff_svg(g,outdir = outdir,filename = "test_heatmap.ggplot")


# summarize by class of sample 3 groups, make heatmap and boxplot
