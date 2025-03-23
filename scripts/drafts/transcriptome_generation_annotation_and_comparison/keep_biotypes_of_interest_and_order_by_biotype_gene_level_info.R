source("scripts/load_functions_and_data_for_overlap_marks.R")
biots_of_interest
gene_level_info <- gene_level_info %>% filter(biotype%in%biots_of_interest)

gene_level_info <- gene_level_info %>% arrange(match(biotype,biots_of_interest))

write.table(gene_level_info,
            gsub(".tsv",
                 ".biotypes_of_interest.tsv",
                 gene_level_data_path),
            quote = F,
            row.names = F,
            sep = "\t")



