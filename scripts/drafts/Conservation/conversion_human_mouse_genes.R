convert_mouse_to_human <- function(gene_list) {
  output = c()
  mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

  for(gene in gene_list) {
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "mouse, laboratory"))[['DB.Class.Key']]
    if( !identical(class_key, integer(0)) ) {
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes) {
        output = rbind(c(gene, human_gene), output)
      }
    }
  }
  return (output)
}

convert_human_to_mouse <- function(gene_list) {
  output = c()
  mouse_human_genes = read.csv("~/HOM_MouseHumanSequence.rpt",sep = "\t")

  for(gene in gene_list) {
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "human"))[['DB.Class.Key']]
    if( !identical(class_key, integer(0)) ) {
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(human_gene in human_genes) {
        output = rbind(c(gene, human_gene), output)
      }
    }
  }
  return (output)
}

abc_ub_genes_tr_mouse <- convert_human_to_mouse(abc_ub_genes$V1)

mouse_human_genes <- read.csv("~/HOM_MouseHumanSequence.rpt",sep = "\t")
conversion_table <- fread('/run/user/1608803857/gvfs/smb-share:server=10.110.20.13,share=investigacio/Cuartero Group/CUARTERO GROUP/references/conversion_tables_human_mouse/Human_Ensembl_genes_104_GRCh38.p13_NCBI_and_mouse_orthologs.tsv',data.table = F)

table(abc_ub_genes$V1%in%c(conversion_table$`Gene name.x`,conversion_table$`Gene name.y`))
table(abc_ub_genes$V1%in%mouse_human_genes$Symbol)
not_in_mine=abc_ub_genes$V1[!abc_ub_genes$V1%in%c(conversion_table$`Gene name.x`,conversion_table$`Gene name.y`)]
not_in_pub=abc_ub_genes$V1[!abc_ub_genes$V1%in%c(mouse_human_genes$Symbol)]
table(not_in_pub%in%not_in_mine)
table(abc_ub_genes_tr_mouse%in%hrt_ub_genes$Gene)

Human_Mouse_common_hrt <- read.csv('/home/llorenzi/Descargas/Human_Mouse_Common.csv',sep = ";")
table(abc_ub_genes$V1%in%Human_Mouse_common_hrt$Human)


hrt_human_ub_genes <- read.csv("~/Descargas/Housekeeping_GenesHuman.csv",sep = ";")
table(abc_ub_genes$V1%in%hrt_human_ub_genes$Gene.name)
abc_ub_genes$V1[!abc_ub_genes$V1%in%hrt_human_ub_genes$Gene.name]

