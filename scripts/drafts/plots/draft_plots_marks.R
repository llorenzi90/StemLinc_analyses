ggplot(gene_level_info,aes(closest_CAGE,col=biotype)) + geom_density()
biots_of_interest=c("protein_coding","lncRNA","TEC","potNovel","pseudogene")
ncbiots=c("lncRNA","TEC","potNovel","pseudogene")
lncRNA_potNovel_TEC=c("lncRNA","TEC","potNovel")
potNovel_TEC=c("TEC","potNovel")


biotypes2select=list(biots_of_interest,
                     ncbiots,
                     lncRNA_potNovel_TEC,
                     potNovel_TEC)

ggplot(gene_level_info%>%filter(biotype%in%biots_of_interest),
       aes(closest_CAGE,col=biotype)) + geom_density() + xlim(c(-10000,10000))
ggplot(gene_level_info%>%filter(biotype%in%biots_of_interest),
       aes(closest_CAGE,col=biotype)) + geom_density() + xlim(c(-500,500))+
  scale_color_manual(values = nejm_pal)

biots_of_interest

g=ggplot(gene_level_info%>%filter(biotype%in%biots_of_interest)%>%filter(biotype!="protein_coding"),
       aes(closest_CAGE,col=biotype)) + geom_density(linewidth = 1) + xlim(c(-10000,10000)) +
  theme_minimal()

g

save_tiff_svg(plo = g,outdir = "outputs/plots/",
              filename = "closest_CAGE_lncRNA_potNovel_pseudogene_TEC",h=5,w=5)


g=ggplot(gene_level_info%>%filter(biotype%in%biots_of_interest)%>%
           filter(!biotype%in%c("protein_coding","lncRNA","pseudogene")),
         aes(closest_CAGE,col=biotype)) + geom_density(linewidth = 1) + xlim(c(-10000,10000)) +
  theme_minimal()

g


ggplot(gene_level_info%>%filter(biotype%in%biots_of_interest), aes(closest_Enhancer_Atlas,col=biotype)) + geom_density(linewidth = 1) + xlim(c(-10000,10000)) +
  theme_minimal() +scale_color_manual(values = nejm_pal)

g=ggplot(gene_level_info%>%filter(biotype%in%biots_of_interest)%>%
           filter(!biotype%in%c("protein_coding","lncRNA","pseudogene")),
         aes(closest_Enhancer_Atlas,col=biotype)) + geom_density(linewidth = 1) + xlim(c(-10000,10000)) +
  theme_minimal() +scale_color_manual(values = nejm_pal)

g


gene_level_info <- left_join(gene_level_info, distance_to_closest_enhancer2)

for (biots in biotypes2select) {
  g=ggplot(gene_level_info%>%filter(biotype%in%biots)
             ,
           aes(closest_FANTOM_enhancer,col=biotype)) + geom_density(linewidth = 1) + xlim(c(-10000,10000)) +
    theme_minimal() +scale_color_manual(values = nejm_pal)
  print(g)

}
