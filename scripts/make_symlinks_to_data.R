list_render_params=list(Kli_MarkedDups_guided =
                          list(tracking_path = '/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/analyses/gffcompare_Kli_MarkedDups_guided/Kli_MarkedDups_guided.tracking',
                               gtf_path= '/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/analyses/gffcompare_Kli_MarkedDups_guided/Kli_MarkedDups_guided.combined.gtf'),

                        LSK_StemLinc.guided =
                          list(tracking_path= '/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare/LSK_StemLinc.tracking' ,
                               gtf_path ='/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare/LSK_StemLinc.combined.gtf' ),

                        LSK_StemLinc.blind =
                          list(tracking_path="/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare_blind/LSK_StemLinc.blind.tracking",
                               gtf_path="/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare_blind/LSK_StemLinc.blind.combined.gtf"))
list_render_params

linkdir="data/raw/"
lapply(list_render_params,function(fi)sapply(paste0("ln -s '",fi,"' ",linkdir),function(cm)system(cm)))
