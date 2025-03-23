# load data ----
source("scripts/load_gene_level_data_biotypes_of_interest.R")
#exp_path="outputs/expression_data/LSK_StemLinc.combined.20240930_173216.filtered.gene_name.all_blood_cells.primary.txt"
exp_path="outputs/expression_data/featureCounts_jan25/LSK_StemLinc.combined.filtered.20250127_112131.all_blood_cells.txt"
exp_dir="outputs/expression_data/featureCounts_jan25/"
featC <- read.table(exp_path,header = T)
head(featC)
#modify sample names ----
inputbams=read.table("outputs/expression_data/featureCounts_jan25/all_public_and_SL_bams_jan2025.txt")
cnames=gsub("Aligned.sortedByCoord.out.MarkedDups.minq1.primarychrs.bam","",basename(inputbams$V1))
cnames[84:87] <- c("LT-HSC_StL_1","MPP3-4_StL_1","MPP3-4_StL_3","MPP3-4_StL_2")
cnames[1:3] <- c("LSK_StL_3","LSK_StL_2","LSK_StL_1")
cnames[88:93] <- c("macrophage_StL_1","macrophage_StL_2","macrophage_StL_3",
                   "T-cell_StL_1","T-cell_StL_2","T-cell_StL_3")
colnames(featC)[7:ncol(featC)]=cnames

# #
# rownames(featC)=featC$Geneid
#
# # separate only count values
# counts_data <- featC[,7:ncol(featC)]
#
# colnames(counts_data) <- gsub("-",".",colnames(counts_data))
#
# ## create sample metadata ----
# sample_type=sapply(strsplit(colnames(counts_data),split = "_"),function(x)x[[1]])
# study=sapply(strsplit(colnames(counts_data),split = "_"),function(x)x[[2]])
# rep=sapply(strsplit(colnames(counts_data),split = "_"),function(x)x[[3]])
# coldata=data.frame(sample=colnames(counts_data),
#                    sample_type=sample_type,
#                    study=study,
#                    rep=rep)
# rownames(coldata)=colnames(counts_data)
#
# ## sort samples first by sample type ----
# sample_type_order=c("LT.HSC","LSK","MPP1","MPP2",
#                     "MPP3","MPP3.4","MPP4","CMP",
#                     "CLP","MEP","GMP","megakaryocyte",
#                     "erythrocyte","granulocyte","monocyte",
#                     "macrophage","dendritic","NK","PreB",
#                     "ProB","B.cell","T.cell")
#
# table(sample_type%in%sample_type_order)
# sample_type[!sample_type%in%sample_type_order]#run this
# #and modify until sample_type is all in sample_type_order
# coldata=coldata[order(match(sample_type,sample_type_order),
#                       study,rep),]
#
# counts_data <- counts_data[,match(coldata$sample,colnames(counts_data))]
#
# ## write featureCounts output with all genes and cleaned samples ----
# featC <- cbind(featC[,1:6],
#                counts_data)
# out_path=gsub("txt","ordered_samples.txt",exp_path)
#
# write.table(featC,out_path,quote = F,row.names = F,
#             sep = "\t")
# # keep only biotypes of interest ----
#
# featC <- featC %>% filter(Geneid%in%gene_level_info$gene_name)
# featC <- featC %>% arrange(match(Geneid,
#                                  gene_level_info$gene_name
# ))
#
# ## write featureCounts output with filtered biotypes genes and cleaned samples ----
#
# out_path=gsub("txt","ordered_samples.selected_biotypes.txt",exp_path)
#
# write.table(featC,out_path,quote = F,row.names = F,
#             sep = "\t")
#
# ## write counts with filtered biotypes genes and cleaned samples ----
#
# counts_data <- featC[,7:ncol(featC)]
# counts_data2write <- cbind(featC%>%select(Geneid),
#                            counts_data)
# out_path=paste0(exp_dir,"featureCounts_",get_timestamp(exp_path),
#                 "ordered_samples.selected_biotypes.counts.txt")
#
# write.table(counts_data2write,
#             out_path,
#             sep = "\t",quote = F,row.names = F)
# # calculate TPMs ----
# # 1. Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
# RPK=apply(counts_data,2,function(x)x/featC$Length/1000)
# # 2. Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
# perMscaling=apply(RPK, 2, function(x)sum(x)/1000000)
# # 3. Divide the RPK values by the “per million” scaling factor. This gives you TPM.
# TPM_data=t(apply(RPK, 1, function(x)x/perMscaling))
#
# TPM_data=as.data.frame(TPM_data)
# ## write TPM data ----
# TPM_data2write <- cbind(featC%>%select(Geneid),
#                         TPM_data)
# out_path=paste0(exp_dir,"featureCounts_",get_timestamp(exp_path),
#                 "ordered_samples.selected_biotypes.TPM.txt")
#
# write.table(TPM_data2write,
#             out_path,
#             sep = "\t",quote = F,row.names = F)
#
#
# # vst ----
# library(DESeq2)
# vst_data=varianceStabilizingTransformation(as.matrix(counts_data))
#
# ## write vst data ----
# vst_data2write <- cbind(featC%>%select(Geneid),
#                         vst_data)
# out_path=paste0(exp_dir,"featureCounts_",get_timestamp(exp_path),
#                 "ordered_samples.selected_biotypes.vst.txt")
#
# write.table(vst_data2write,
#             out_path,
#             sep = "\t",quote = F,row.names = F)
#
#
# # add sample class to coldata ----
# coldata <- coldata %>% mutate(cell_class = ifelse(grepl("LT|LSK|MPP",sample),
#                                                   "HSPC",
#                                                   ifelse(grepl("CMP|CLP|MEP|GMP",sample),
#                                                          "progenitor",
#                                                          "differentiated")))
#
#
# # write coldata ----
# out_path=paste0(exp_dir,"coldata_",get_timestamp(exp_path),
#                 "ordered_samples.txt")
#
# write.table(coldata,
#             out_path,
#             sep = "\t",quote = F,row.names = F)
