source("scripts/load_gene_level_data_biotypes_of_interest.R")


datadir <- "outputs/overlap_marks/Enhancers/"
sample_pref <- "LSK_StemLinc.combined_annotated_tracking.gene_level_info.20240930_173216.filtered.gene_classif."


outdir="outputs/plots/Enhancers"
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

fraction_plot_within <- function(plot_data,
                          dist_co=500,
                          save_plot=T,
                          return_data=F,
                          colors=nejm_pal){

  plot_data <- plot_data %>% group_by(biotype) %>%
    summarise(nmark=sum(abs(mark)<=dist_co),
              fraction=nmark/n())

  frac_plot <- plot_data %>%
    ggplot(aes(x=biotype, y=fraction, fill=biotype)) +
    geom_col() +
    scale_fill_manual(values = colors) +
    ggtitle(paste("Fraction of genes with",
                  mark,"peak within",dist_co,"of TSS")) +
    theme(axis.text.x = element_text(angle = 45 , hjust = 1))


  if(save_plot){
    save_tiff_svg(frac_plot,outdir = outdir,
                  filename =
                    paste0("Fraction_",mark,"_within_",dist_co,"_TSS" ))
  }

  return_list=list(frac_plot)
  if (return_data) {
    return_list=c(return_list, plot_data)
  }

  return(return_list)
}

plot_closest_mark <- function(data=gene_level_info,
                              mark,
                              outdir="outputs/plots/",
                              outfile="plot",
                              filter_var= biotype,
                              col_var = {{filter_var}},
                              filter_list=biotypes2select,
                              xlims=c(-500,500),
                              color_pal=nejm_pal,
                              h=5,
                              w=8,
                              save_plot=T){
  for (biots in filter_list) {
    tmp_outfile=paste0(outfile,".",paste0(biots,collapse = "_"),"_",paste0(xlims,collapse = "_"))
    g=ggplot(data%>%filter({{filter_var}}%in%biots)
             ,
             aes({{mark}},col={{col_var}})) + geom_density(linewidth = 1) +
      xlim(xlims) +
      theme_minimal() +scale_color_manual(values = color_pal)
    print(g)
    if(save_plot) save_tiff_svg(g,outdir = outdir,filename = tmp_outfile, h=h,w=w)

  }
}
# plots enhancers ----

for(mark in c("Enhancer_Atlas","FANTOM_enhancers")){
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

  # Distance to closest enhancer ----
  datatype="closest"
  data_path=paste0(datadir,sample_pref,datatype,"_",mark,".tsv")
  da=read.table(data_path,header = T)

  gene_level_info <- left_join(gene_level_info, da)
  feat_name=colnames(gene_level_info)[ncol(gene_level_info)]

  colnames(gene_level_info)[ncol(gene_level_info)]="closest_enhancer"
  ### distance closest Enhancer ----
  gene_level_info %>%
    plot_closest_mark(mark = closest_enhancer,
                      outfile = paste0("closest_",mark),
                      outdir = outdir,
                      save_plot = T)

  ### by gene classif ----

  gene_level_info %>% filter(biotype%in%ncbiots) %>%
    plot_closest_mark(mark = closest_enhancer,
                      outdir = outdir,
                      outfile = paste0("closest_",mark,"_ncBiotypes_classif"),
                      filter_var = classif,
                      filter_list = list(unique(gene_level_info$classif)),
                      save_plot = T)

  # fraction of TSS with enhancer peaks within 500 bp ----
  colnames(gene_level_info)[ncol(gene_level_info)]="mark"
  fraction_plot_within(gene_level_info,
                       save_plot = T)

  ## return feature name to original ----
  colnames(gene_level_info)[ncol(gene_level_info)]=feat_name
}




# Enhancers defined by marks ----

# write extended gene level info ----
write_info_table(gene_level_info,data_path = gene_level_data_path,modify_path = T,
                 suffix = "enhancers")
