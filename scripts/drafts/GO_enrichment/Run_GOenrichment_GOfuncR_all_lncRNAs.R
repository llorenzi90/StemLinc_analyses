## Head -------------------------------------
##
##
## Purpose of script: Given a list of lncRNAs and correlated 
##                    PCGs generate a GO table with enriched GO terms
##                    and/or associated terms for each PCG
## Author: Lucia Lorenzi
##
## Date Created: 2024-03-13
##
## Email: lucialorenzi90@gmail.com
##
## Notes ---------------------------
##  INPUTS:
##          - corr table (rows PCGs, cols lncRNAs)
##          - threshold
##            
##
## Setup ---------------------------

options(scipen = 999) 
## Load packages---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
require(GOfuncR)
library(org.Mm.eg.db)
## Load data---------------------------
default_thresh=0.6
args=commandArgs(trailingOnly = T)
if(length(args)==0){
  args=c(
    "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/correlation_PCG_lncRNA/corr_lncRNA_PCG.16brownPotNovel.txt",
    default_thresh)
}else{if(length(args)==1) args[2]=default_thresh}

gene_translation_file <- "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/featureCounts/PCGslncRNAs.mapIDfile.txt"
gene_translation <- fread(gene_translation_file)

# some setup data ----
# select all PCG gene_ids as universe
gene_universe=gene_translation$gene_name[!gene_translation$biotype%in%c("TEC","lncRNA")]
gene_universe <- na.omit(gene_universe)


# define functions ----

get_genes_corrvals <- function(ln){
  cn=which(colnames(cor_PCG_lncRNA)==ln)
  corix=which(abs(cor_PCG_lncRNA[,cn])>=cor_thresh)
  corr=as.numeric(cor_PCG_lncRNA[corix,cn])
  gene_ids=rownames(cor_PCG_lncRNA)[corix]
  gene_names=gene_translation$gene_name[match(gene_ids,
                                              gene_translation$gene_id_vM31)]
  return(data.frame(corr,gene_names))
}


GOanalysis <- function(gene_names,
                       corr_vals,
                       gene_universe=gene_universe,
                       gene_translation=gene_translation,
                       orgDb="org.Mm.eg.db"){
  
  entrez_id=na.omit(gene_translation$ncbi_GeneID[match(gene_names,
                                               gene_translation$gene_name)])
  
  go_res_annot=NULL
  
  if(sum(!is.na(entrez_id))==0){
    print("no genes with ENTREZID")
  }else{
    input_hyper_mouse = data.frame(gene_id=gene_universe, 
                                   is_candidate=ifelse(gene_universe%in%gene_names,1,0))
    
    
    # test the gene set for enrichment in GO categories
    res_hyper_mouse = go_enrich(input_hyper_mouse, orgDb = orgDb)
    
    go_res=res_hyper_mouse$results
    
    # map genes to GO terms
    anno_genes = get_anno_genes(go_ids=go_res$node_id, 
                                genes=gene_names,database = orgDb)
    
    anno_genes$GO_name=go_res$node_name[match(anno_genes$go_id,
                                              go_res$node_id)]
    
    anno_genes <- anno_genes %>% filter(!GO_name%in%c("molecular_function",
                                                      "cellular_component",
                                                      "biological_process"))
    
    anno_genes$corr=corr_vals[match(anno_genes$gene,
                                                 gene_names)]
    #add info on genes associated with each GO term
    
    anno_genes <- anno_genes %>% arrange(-abs(corr))
    
    
    summ_go_ids <- anno_genes %>% group_by(go_id) %>% 
      summarise(gene_names=paste0(gene,collapse = ","),
                N_genes=n(),
                corr_vals=paste0(round(corr,2),collapse = ","))
    
    
    #add everuthing in a single table
    go_res_annot=inner_join(go_res,summ_go_ids,by=c(node_id="go_id"))
    
    }
    
   return(go_res_annot)
}


# run code ----
fi=args[1]
#fi="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/correlation_PCG_lncRNA/corr_lncRNA_PCG.9extraPotNovel.txt"
fi="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/Lucia/Marie_Curie/Experimento4_LSK_190224/analyses/correlation_PCG_lncRNA/corr_lncRNA_PCG.txt"

cor_PCG_lncRNA=fread(fi,data.table = F)
cor_thresh=args[2]

all_corr_genes <- lapply(colnames(cor_PCG_lncRNA),get_genes_corrvals)
names(all_corr_genes) <- colnames(cor_PCG_lncRNA)

N_corr=apply(cor_PCG_lncRNA,2,function(x)sum(abs(x)>=cor_thresh))
allres=lapply(names(all_corr_genes),function(x){
  da=all_corr_genes[[x]]
  gores=GOanalysis(da$gene_names,da$corr,gene_universe = gene_universe,
                   gene_translation = gene_translation,orgDb = "org.Mm.eg.db")
  gores$lncRNA=x
  gores$totalNcorr=length(da$gene_names)
  gores <- gores[,c((ncol(gores)-1),ncol(gores),1:(ncol(gores)-2))]
  return(gores)
})
# gather and write results ----
allres=do.call("rbind",allres)

write.table(allres,paste0(tools::file_path_sans_ext(basename(fi)),".GOenrich_results.",cor_thresh,".txt"))
