save_tiff_svg <- function(plo,outdir=".",filename="plot",
                       h=10,w=10){
  grDevices::tiff(res = 300, height = h ,width = w, units = "in",
       paste0(outdir,"/",filename,".tiff"))
  print(plo)
  dev.off()

  grDevices::svg( height = h ,width = w,
       paste0(outdir,"/",filename,".svg"))
  print(plo)
  dev.off()

}

nejm_pal=ggsci::pal_nejm()(8)

# functions useful for ordering boxplots within facets
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}



scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}



scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

#source: https://github.com/dgrtwo/drlib/blob/master/R/reorder_within.R
