#load data ----
source("scripts/source_all_functions.R")
source("scripts/load_gene_level_data_biotypes_of_interest.R")

gene_level_info <- read.csv("outputs/gene_level_info_all_evidence_oct24.csv")
g=ggplot(gene_level_info,aes(x=biotype,fill=exonic_type,y=N_TFs)) + geom_boxplot() +
  scale_fill_manual(values = nejm_pal) + theme_minimal() + ylab("# TFs at promoter")
g

save_tiff_svg(g,outdir = outdir,"N_TF_all_genes")

gene_level_info %>% group_by(biotype,exonic_type) %>% summarise(mean(N_TFs))

# N TFs at promoters selected nc genes ----
gene_level_info_sel <- read.csv("outputs/gene_level_info_all_evidence_oct24.selected114.csv")


g=ggplot(gene_level_info_sel,aes(x=biotype,fill=exonic_type,y=N_TFs)) + geom_boxplot() +
  scale_fill_manual(values = nejm_pal) + theme_minimal() + ylab("# TFs at promoter")
g

save_tiff_svg(g,outdir = outdir,"N_TF_114selected_genes")
gene_level_info_sel %>% group_by(biotype,exonic_type) %>% summarise(mean(N_TFs))

g=gene_level_info_sel %>% ggplot(aes( biotype, fill=classif)) + geom_bar() +facet_wrap(~exonic_type)+
  scale_fill_manual(values = nejm_pal) +theme_minimal()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
g
save_tiff_svg(g,outdir,"class_biotype_distribution_114_highcorr")

# Venn 3 enhancer sources selected genes ----
logic_marks=read.delim("outputs/overlap_marks/gene_level_info.20240930_173216.logicmarks.tsv")

logic_marks <- left_join(gene_level_info_sel%>%select(gene_name,biotype),logic_marks)

ggplot(logic_marks,aes(x=biotype,fill=FANTOM_enhancers_at_promoter)) +geom_bar()
ggplot(logic_marks,aes(x=biotype,fill=Enhancer_Atlas_at_promoter)) +geom_bar()
ggplot(logic_marks,aes(x=biotype,fill=Enhancer_from_marks_at_promoter)) +geom_bar()

fgenes=logic_marks$gene_name[logic_marks$FANTOM_enhancers_at_promoter]
agenes=logic_marks$gene_name[logic_marks$Enhancer_Atlas_at_promoter]
mgenes=logic_marks$gene_name[logic_marks$Enhancer_from_marks_at_promoter]

library(VennDiagram)
venn.diagram(
  x = list(FANTOM = fgenes, Enhancer_Atlas = agenes, Enhancer_marks = mgenes),
  category.names = c("", "", ""),
  filename = "outputs/plots/venn_diagram_enhancers_selected.png",  # Save as PNG
  output = TRUE,
  fill = c("#BC3C29FF", "#0072B5FF",  "#20854EFF"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 1.5,
  cat.pos = c(-45, 45, 180),
  cat.dist = c(0.05, 0.05, 0.05)
)
#"#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF" "#FFDC91FF" "#EE4C97FF

# Bar plot any enhancer seleced genes ----
g=ggplot(gene_level_info_sel,aes(biotype,fill=Enhancer)) + geom_bar() +
  scale_fill_manual(values = c("#E18727FF", "#0072B5FF")) + theme_minimal()
g
save_tiff_svg(g,outdir = outdir, "bar_plot_any_enhancer_selected114")

# Bar plot evidence level facet exonic type ----
g=ggplot(gene_level_info_sel,aes(biotype,fill=evidence_level)) +geom_bar() +
  scale_fill_manual(values = c("#E18727FF", "#0072B5FF","grey")) +
  theme_minimal() +facet_wrap(~exonic_type) +theme(axis.text.x = element_text(angle = 45, hjust = 1))
g
save_tiff_svg(g,outdir,"level_of_evidence_candidatessel114")

