# list_render_params=list(Kli_MarkedDups_guided =
#                           list(tracking_path = '/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/analyses/gffcompare_Kli_MarkedDups_guided/Kli_MarkedDups_guided.tracking',
#                                gtf_path= '/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/analyses/gffcompare_Kli_MarkedDups_guided/Kli_MarkedDups_guided.combined.gtf'),
#
#                         LSK_StemLinc.guided =
#                           list(tracking_path= '/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare/LSK_StemLinc.tracking' ,
#                                gtf_path ='/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare/LSK_StemLinc.combined.gtf' ),
#
#                         LSK_StemLinc.blind =
#                           list(tracking_path="/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare_blind/LSK_StemLinc.blind.tracking",
#                                gtf_path="/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare_blind/LSK_StemLinc.blind.combined.gtf"))
# list_render_params
#
# linkdir="data/raw/"
# lapply(list_render_params,function(fi)sapply(paste0("ln -s '",fi,"' ",linkdir),function(cm)system(cm)))

linkdir="data/references/merged_refs_annotation/"
#ref_gtf='/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/references/mouse/annotation/merged_refs_annotation/merged_refs.combined.annot_trnames.gtf'
ref_gtf='/run/user/1608803857/gvfs/smb-share:server=10.110.20.7,share=bdcuartero/references/mouse/annotation/merged_refs_annotation/merged_refs.combined.annot_trnames.gtf'
system(paste0("ln -s '",ref_gtf,"' ",linkdir))
annot='/run/user/1608803857/gvfs/smb-share:server=10.110.20.7,share=bdcuartero/references/mouse/annotation/merged_refs_annotation/annotated_tracking_file.updated_gene_names.txt'
system(paste0("ln -s '",annot,"' ",linkdir))
ref_gtf_with_gene_name='/run/user/1608803857/gvfs/smb-share:server=10.110.20.7,share=bdcuartero/references/mouse/annotation/merged_refs_annotation/merged_refs.combined.with_gene_name.gtf'
system(paste0("ln -s '",ref_gtf_with_gene_name,"' ",linkdir))

linkdir="data/references/genomes/"
mm39_genome='/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/references/mouse/genomes/mm39.fa'
system(paste0("ln -s '",mm39_genome,"' ",linkdir))


linkdir="data/public_data/CAGE/"
CAGE_data='/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/FANTOM5_CAGE/TSS_mouse.liftOver.mm39.bed'
system(paste0("ln -s '",CAGE_data,"' ",linkdir))

linkdir="data/public_data/cistromes/"
dir.create(linkdir)
cistromes_data='/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/cistromes/ALL_TFs_mm39_cistrome.bed'
system(paste0("ln -s '",cistromes_data,"' ",linkdir))

linkdir="data/public_data/Enhancer_data/"
dir.create(linkdir)
source_data='/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/Enhancer_data/Enhancer_atlas_LSK.mm39.bed'
system(paste0("ln -s '",source_data,"' ",linkdir))



linkdir="data/public_data/Enhancer_data/"
dir.create(linkdir)
source_data='/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/Enhancer_data/FANTOM5/mouse_permissive_enhancers_phase_1_and_2.mm39.bed'
system(paste0("ln -s '",source_data,"' ",linkdir))


linkdir="data/public_data/ChIP_seq_data/"
dir.create(linkdir,recursive = T)
source_data='/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/ChIP_seq_data/me3_marks_GSE47765_HSCs_Luo/mm39_beds'
for(sd in list.files(source_data,full.names = T)){
  system(paste0("ln -s '",sd,"' ",linkdir))

}


# H3K27ac
linkdir="data/public_data/ChIP_seq_data/"
dir.create(linkdir,recursive = T)
source_data='/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/ChIP_seq_data/K27ac_K4me1_GSE63276_RAW/mm39_beds/GSM1544999_m24_H3K27ac.narrowPeak.mm39.bed'
system(paste0("ln -s '",source_data,"' ",linkdir))


# H3K4me1
linkdir="data/public_data/ChIP_seq_data/"
dir.create(linkdir,recursive = T)
source_data='/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/ChIP_seq_data/K27ac_K4me1_GSE63276_RAW/mm39_beds/GSM1545000_m24_H3K4me1.narrowPeak.mm39.bed'
system(paste0("ln -s '",source_data,"' ",linkdir))

# ployAsite
linkdir="data/public_data/polyAsite2.0/"
dir.create(linkdir,recursive = T)
source_data='/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/polyASite/atlas.clusters.2.0.chr.liftOver.mm39.bed'
system(paste0("ln -s '",source_data,"' ",linkdir))


# Delas data
linkdir="data/public_data/lncRNA_studies/Delas_et_al"
dir.create(linkdir,recursive = T)
source_data <- '/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/metadata_RNA_seq_studies/Delas_et_al/Delas_all_lncRNAs_2017_from_website.txt'
system(paste0("ln -s '",source_data,"' ",linkdir))

source_data <- '/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/metadata_RNA_seq_studies/Delas_et_al/elife-25607-supp1-v2_Delas_lncRNA_catalog.liftovermm39.gtf'
system(paste0("ln -s '",source_data,"' ",linkdir))

linkdir="data/public_data/lncRNA_studies/Delas_et_al"
source_data <- '/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/strtie_transcriptome_assembly/gffcompare/StemLincpotNovLSK_vs_Delas/StemLincpotNovLSK_vs_Delas.LSK_StemLinc.combined.gtf.tmap'
system(paste0("ln -s '",source_data,"' ",linkdir))

source_data <- '/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/metadata_RNA_seq_studies/Delas_et_al/elife-25607-supp2-v2.xls'
system(paste0("ln -s '",source_data,"' ",linkdir))

linkdir="data/public_data/lncRNA_studies/Luo_et_al"
dir.create(linkdir)
source_data <- '/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/metadata_RNA_seq_studies/Goodell/annotations/Goodell_159_LncHSCs.liftOver.mm39.bed'
system(paste0("ln -s '",source_data,"' ",linkdir))

source_data <- '/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/metadata_RNA_seq_studies/Goodell/503_novel_lncRNAs.liftOver.mm39.bed'
system(paste0("ln -s '",source_data,"' ",linkdir))


linkdir="data/references/conservation_bigWigs/"
dir.create(linkdir)
source_data <- '/home/llorenzi/Documentos/references/mm39.phastCons35way.bw'
system(paste0("ln -s '",source_data,"' ",linkdir))

source_data <- '/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/public_datasets/RepeatMasker/UCSC_RepeatMasker_mm39.bed'
linkdir="data/references/RepeatMasker/"
dir.create(linkdir)
system(paste0("ln -s '",source_data,"' ",linkdir))


list_render_params=list(T_cell.guided =
                          list(tracking_path = '/home/llorenzi/Documentos/T-cells_macrophages/gffcompare/T_cell.tracking',
                               gtf_path= '/home/llorenzi/Documentos/T-cells_macrophages/gffcompare/T_cell.combined.gtf'),

                        Macro.guided =
                          list(tracking_path = '/home/llorenzi/Documentos/T-cells_macrophages/gffcompare/Macro.tracking',
                               gtf_path= '/home/llorenzi/Documentos/T-cells_macrophages/gffcompare/Macro.combined.gtf'),

                        T_cell.blind =
                          list(tracking_path = '/home/llorenzi/Documentos/T-cells_macrophages/gffcompare_blind/T_cell.blind.tracking',
                               gtf_path= '/home/llorenzi/Documentos/T-cells_macrophages/gffcompare_blind/T_cell.blind.combined.gtf'),

                        Macro.blind =
                          list(tracking_path = '/home/llorenzi/Documentos/T-cells_macrophages/gffcompare_blind/Macro.blind.tracking',
                               gtf_path= '/home/llorenzi/Documentos/T-cells_macrophages/gffcompare_blind/Macro.blind.combined.gtf')
)
list_render_params

linkdir="data/raw/"
lapply(list_render_params,function(fi)sapply(paste0("ln -s '",fi,"' ",linkdir),function(cm)system(cm)))
