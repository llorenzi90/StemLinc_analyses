## plot density ----
outdir="outputs/plots/"
gene_level_info <- read.csv("outputs/gene_level_info_all_evidence_oct24.csv")
gene_level_info_nc <- gene_level_info %>%filter(biotype!="protein_coding")
g=ggplot(gene_level_info_nc,
         aes(x=mean_corr_signature)) + geom_density(linewidth = 2) +theme_minimal()

g
save_tiff_svg(g,outdir = outdir,filename = "density_mean_corr_signature",h = 5,w=8)

g=g + geom_vline(xintercept = 0.65,linetype = "dashed")

g
save_tiff_svg(g,outdir = outdir,filename = "density_mean_corr_signature_dashed0.65",h = 5,w=8)

