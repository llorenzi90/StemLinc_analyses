datadir="data/raw"

gtf_files=list.files(datadir,pattern = "gtf",full.names = T)


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
