datadir="data/raw"

gtf_files=list.files(datadir,pattern = "gtf$",full.names = T)


list_render_params <- lapply(gtf_files, function(gtf_path){

  return(list(gtf_path = gtf_path,
              tracking_path = gsub(".combined.gtf",".tracking",gtf_path)))})

sampnames=basename(gsub(".combined", "", gsub(".gtf", "", gtf_files)))

names(list_render_params) <- sampnames

list_render_params <- list_render_params[names(list_render_params)!="LSK_StemLinc.sorted"]

script_path="Transcriptome_characterization_report.Rmd"

for (pref in names(list_render_params)) {

  tmp_list_params=list_render_params[[pref]]
  rmarkdown::render(input = script_path,params = tmp_list_params,
                    output_file =paste0("outputs/transcriptome_characterization/",
                                        pref,".transcriptome_characteristics_report.html" ))


}


# T-cell and macrophages
datadir="data/raw"

gtf_files=list.files(datadir,pattern = "gtf$",full.names = T)

gtf_files <- grep("T_cell|Macro",gtf_files,value = T)

list_render_params <- lapply(gtf_files, function(gtf_path){

  return(list(gtf_path = gtf_path,
              tracking_path = gsub(".combined.gtf",".tracking",gtf_path)))})

sampnames=basename(gsub(".combined", "", gsub(".gtf", "", gtf_files)))

names(list_render_params) <- sampnames


script_path="Transcriptome_characterization_report.Rmd"

for (pref in names(list_render_params)) {

  tmp_list_params=list_render_params[[pref]]
  rmarkdown::render(input = script_path,params = tmp_list_params,
                    output_file =paste0("outputs/transcriptome_characterization/",
                                        pref,".transcriptome_characteristics_report.html" ))


}


# Combined transcriptome from filtered LSK, T-cells and macrophages

datadir="data/hpc_data/"

gtf_files=list.files(datadir,pattern = "gtf$",full.names = T)

gtf_files <- grep("filtered_LSK_T-cell_macrophage",gtf_files,value = T)

list_render_params <- lapply(gtf_files, function(gtf_path){

  return(list(gtf_path = gtf_path,
              tracking_path = gsub(".combined.gtf",".tracking",gtf_path)))})

sampnames=basename(gsub(".combined", "", gsub(".gtf", "", gtf_files)))

names(list_render_params) <- sampnames


script_path="Transcriptome_characterization_report.Rmd"

for (pref in names(list_render_params)) {

  tmp_list_params=list_render_params[[pref]]
  rmarkdown::render(input = script_path,params = tmp_list_params,
                    output_file =paste0("outputs/transcriptome_characterization/",
                                        pref,".transcriptome_characteristics_report.html" ))


}
