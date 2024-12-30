source("scripts/load_gene_level_data_biotypes_of_interest.R")


datadir <- "outputs/overlap_marks/histone_marks/"
sample_pref <- "LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif."


outdir="outputs/plots/histone_marks"
dir.create(outdir)

# functions ----
ecdf_plot <- function(plot_data, save_plot=T) {

  g <- ggplot(plot_data, aes(x=mark, col=biotype))+
    stat_ecdf(linewidth=2) + theme_minimal() +
    scale_color_manual(values = nejm_pal) +
    ggtitle(paste0("Fraction of ", datatype, " covered by ",mark," signal"))

  if(save_plot)  save_tiff_svg(g,outdir = outdir,filename = paste0("ecdf_",mark,"_",datatype))
  return(g)
}

fraction_plot_cov <- function(plot_data,
                          frac_co=0.5,
                          save_plot=T,
                          return_data=T,
                          colors=nejm_pal){

  plot_data <- plot_data %>% group_by(biotype) %>%
    summarise(nmark=sum(mark>=frac_co),
              fraction=nmark/n())

  frac_plot <- plot_data %>%
    ggplot(aes(x=biotype, y=fraction, fill=biotype)) +
    geom_col() +
    scale_fill_manual(values = colors) +
    ggtitle(paste("Fraction of genes with",
                  frac_co,
                  mark,"coverage at",datatype)) +
    theme(axis.text.x = element_text(angle = 45 , hjust = 1))


 if(save_plot){
   save_tiff_svg(frac_plot,outdir = outdir,
                 filename =
                   paste0("Fraction",mark,"_at_least_",frac_co,"at",datatype ))
 }

  return_list=list(frac_plot)

  if (return_data) {
    return_list=c(return_list, plot_data)
  }

  return(return_list)
}

fraction_plot <- function(plot_data,

                              save_plot=T,
                              return_data=T,
                              colors=nejm_pal){

  plot_data <- plot_data %>% group_by(biotype) %>%
    summarise(nmark=sum(mark),
              fraction=nmark/n())

  frac_plot <- plot_data %>%
    ggplot(aes(x=biotype, y=fraction, fill=biotype)) +
    geom_col() +
    scale_fill_manual(values = colors) +
    ggtitle(paste("Fraction of genes with",
                  mark,"at",datatype)) +
    theme(axis.text.x = element_text(angle = 45 , hjust = 1))


  if(save_plot){
    save_tiff_svg(frac_plot,outdir = outdir,
                  filename =
                    paste0("Fraction_",mark,"_at_",datatype ))
  }

  return_list=list(frac_plot)
  if (return_data) {
    return_list=c(return_list, plot_data)
  }

  return(return_list)
}

# plots histone marks ----

mark="H3K4me3"
for(mark in c("H3K4me3","H3K36me3","H3K27me3","H3K27ac","H3K4me1")){
  # TSS ----
  datatype="TSS"
  data_path<- paste0(datadir,sample_pref,mark,"at",datatype,".tsv")
  da <- read.table(data_path, header = T)
  gene_level_info <- left_join(gene_level_info, da)
  cols=c(colnames(da),"biotype")
  plot_da <- gene_level_info[,cols]
  colnames(plot_da)[2]="mark"
  fraction_plot(plot_da)

  # Promoters ----
  datatype="Promoters"
  data_path<- paste0(datadir,sample_pref,mark,"at",datatype,".tsv")
  da <- read.table(data_path, header = T)
  gene_level_info <- left_join(gene_level_info, da)
  cols=c(colnames(da),"biotype")
  plot_da <- gene_level_info[,cols]
  colnames(plot_da)[2]="mark"
  fraction_plot_cov(plot_da)

  # GeneBody ----
  datatype="GeneBody"
  data_path<- paste0(datadir,sample_pref,mark,"at",datatype,".tsv")
  da <- read.table(data_path, header = T)
  colnames(da)
  # [1] "gene_name"                      "H3K4me3_total_cov"
  # [3] "H3K4me3_total_len"              "H3K4me3_total_fraction_covered"
  # [5] "H3K4me3_max_fraction_covered"   "H3K4me3_max_tr_covered"
  # [7] "H3K4me3_max_tr_len"

  # I tipically choose column 4 total fraction covered
  gene_level_info <- left_join(gene_level_info, da)
  cols=c(colnames(da),"biotype")
  plot_da <- gene_level_info[,cols]
  colnames(plot_da)[4]="mark"
  ecdf_plot(plot_da)
  fraction_plot_cov(plot_da)

}

# write extended gene level info ----
write_info_table(gene_level_info,data_path = gene_level_data_path,modify_path = T,
                 suffix = "histone_marks")
