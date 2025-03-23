library(UpSetR)
movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"),
                   header = T, sep = ";")

upset(movies, nsets = 6, number.angles = 30, point.size = 3.5, line.size = 2,
      mainbar.y.label = "Genre Intersections", sets.x.label = "Movies Per Genre",
      text.scale = c(1.3, 1.3, 1, 1, 2, 0.75))
# example of list input (list of named vectors)
listInput <- list(one = c(1, 2, 3, 5, 7, 8, 11, 12, 13), two = c(1, 2, 4, 5,
                                                                 10), three = c(1, 5, 6, 7, 8, 9, 10, 12, 13))

# example of expression input
expressionInput <- c(one = 2, two = 1, three = 2, `one&two` = 1, `one&three` = 4,
                     `two&three` = 1, `one&two&three` = 2)

upset(fromList(listInput), order.by = "freq")

# intronic genes vs number of samples,

intronic_genes <- tracking_GL %>% filter(combined_cc%in%c("i;i","i;NA","NA;i"))
StemLinc_intronic=lapply(1:3,function(i)intronic_genes%>%filter(maxNsamps1==i) %>% pull(gene_name))
Klimmeck_intronic=lapply(1:3,function(i)intronic_genes%>%filter(maxNsamps2==i) %>% pull(gene_name))
intronic_genes_list <- c(StemLinc_intronic,Klimmeck_intronic)
names(intronic_genes_list)=c(paste0("SL_intronic_maxSamps=",1:3),
                             paste0("Kli_intronic_maxSamps=",1:3))

upset(fromList(intronic_genes_list), order.by = "freq",nsets = 6)


intergenic_genes <- tracking_GL %>% filter(combined_cc%in%c("u;u","u;NA","NA;u"))
StemLinc_intergenic=lapply(1:3,function(i)intergenic_genes%>%filter(maxNsamps1==i) %>% pull(gene_name))
Klimmeck_intergenic=lapply(1:3,function(i)intergenic_genes%>%filter(maxNsamps2==i) %>% pull(gene_name))
intergenic_genes_list <- c(StemLinc_intergenic,Klimmeck_intergenic)
names(intergenic_genes_list)=c(paste0("SL_intergenic_maxSamps=",1:3),
                             paste0("Kli_intergenic_maxSamps=",1:3))

upset(fromList(intergenic_genes_list), order.by = "freq",nsets = 6)

# exonic type and overlap Ref

tracking_GL_potNovel <- tracking_GL%>%filter(is.na(biotype))
potNovel_exonic_type_list <- list(StemLinc_mono=tracking_GL_potNovel%>%filter(maxNexons1==1) %>% pull(gene_name),
                                  StemLinc_multi=tracking_GL_potNovel%>%filter(maxNexons1>1) %>% pull(gene_name),
                                  Klimmeck_mono=tracking_GL_potNovel%>%filter(maxNexons2==1) %>% pull(gene_name),
                                  Klimmeck_multi=tracking_GL_potNovel%>%filter(maxNexons2>1) %>% pull(gene_name))
potNovel_exonic_type <- fromList(potNovel_exonic_type_list)
up=upset(fromList(potNovel_exonic_type_list), order.by = "freq",nsets = 6, )
up
grid.text("Potential novel genes and exonic type",x = 0.65, y=0.95, gp=gpar(fontsize=12))

str(potNovel_exonic_type_list)
# df2 <- data.frame(gene=unique(unlist(potNovel_exonic_type_list)))
#
# head(df2)
#
# usp=upset(fromList(potNovel_exonic_type_list))
# usp$New_data
#
# df1 <- lapply(potNovel_exonic_type_list,function(x){
#   data.frame(gene = x)
# }) %>%
#   bind_rows(.id = "set")
#
# df_int <- lapply(df2$gene,function(x){
#   # pull the name of the intersections
#   intersection <- df1 %>%
#     dplyr::filter(gene==x) %>%
#     arrange(set) %>%
#     pull("set") %>%
#     paste0(collapse = "|")
#
#   # build the dataframe
#   data.frame(gene = x,int = intersection)
# }) %>%
#   bind_rows()
# unlist(potNovel_exonic_type_list, use.names = FALSE)

# add annotation to each set combination in upsetR
# plot expression levels
# automate the gene summarization and script in general

sets=fromList(potNovel_exonic_type_list)
sets$gene_name=unique(unlist(potNovel_exonic_type_list,use.names = F))
sets$maxNsampsStemLinc=tracking_GL_potNovel$maxNsamps1[match(sets$gene_name,
                                                             tracking_GL_potNovel$gene_name)]


upset(sets,
      query.legend = "bottom", nsets = 6, number.angles = 30, point.size = 3.5, line.size = 2,
      queries = list(
        list(
          query = elements,
          params = list("maxNsampsStemLinc"),
          color = "#Df5286",
          active = T,
          query.name = "Maximum replicates StemLinc"
        )
      )
)
sets$maxNsampsStemLinc2=sets$maxNsampsStemLinc==2
sets$maxNsampsStemLinc3=sets$maxNsampsStemLinc==3

upset(sets,
      query.legend = "bottom", nsets = 6, number.angles = 30, point.size = 3.5, line.size = 2,
      queries = list(
        list(
          query = elements,
          params = list("maxNsampsStemLinc1",T),
          color = "#Df5286",
          active = T,
          query.name = "1 Maximum replicate StemLinc"
        )

      )
)

upset(sets,
      query.legend = "bottom", nsets = 6, number.angles = 30, point.size = 3.5, line.size = 2,
      queries = list(
        list(
          query = elements,
          params = list("maxNsampsStemLinc2",T),
          color = "#Df5286",
          active = T,
          query.name = "2 Maximum replicates StemLinc"
        )


      )
)
#
# "#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF"
# [7] "#FFDC91FF" "#EE4C97FF"
#

upset(sets,
      query.legend = "bottom", nsets = 6, number.angles = 30, point.size = 3.5, line.size = 2,
      queries = list(
        list(
          query = elements,
          params = list("maxNsampsStemLinc1",T),
          color = "#EE4C97FF",
          active = T,
          query.name = "1 Maximum replicate StemLinc"
        ),

        list(
          query = elements,
          params = list("maxNsampsStemLinc2",T),
          color =  "#E18727FF",
          active = T,
          query.name = "2 Maximum replicates StemLinc"
        )

      )
)



upset(sets,
      query.legend = "bottom", nsets = 6, number.angles = 30, point.size = 3.5, line.size = 2,
      queries = list(



        list(
          query = elements,
          params = list("maxNsampsStemLinc3",T),
          color =  "#20854EFF",
          active = T,
          query.name = "3 Maximum replicates StemLinc"
        )
      )
)


#

exonic_type_list <- list(StemLinc_mono=tracking_GL%>%filter(maxNexons1==1) %>% pull(gene_name),
                                  StemLinc_multi=tracking_GL%>%filter(maxNexons1>1) %>% pull(gene_name),
                                  Klimmeck_mono=tracking_GL%>%filter(maxNexons2==1) %>% pull(gene_name),
                                  Klimmeck_multi=tracking_GL%>%filter(maxNexons2>1) %>% pull(gene_name))
sets <- fromList(exonic_type_list)
sets$gene_name=unique(unlist(exonic_type_list,use.names = F))
tracking_GL$overlapRef=!is.na(tracking_GL$biotype)
sets$overlapRef=tracking_GL$overlapRef[match(sets$gene_name,
                                             tracking_GL$gene_name)]
up=upset(sets, order.by = "freq",nsets = 6, )
up
grid.text("Potential novel genes and exonic type",x = 0.65, y=0.95, gp=gpar(fontsize=12))

upset(sets,
      query.legend = "bottom", nsets = 6, number.angles = 30, point.size = 3.5, line.size = 2,
      queries = list(

        list(
          query = elements,
          params = list("overlapRef",T),
          color =  "#20854EFF",
          active = T,
          query.name = "Overlaps Reference"
        )
      )
)
