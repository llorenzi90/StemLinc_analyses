get_max_per_group_by_arrange <- function(tracking,
                                         group_var,
                                         arr_var,
                                         summ_var={{arr_var}},
                                         descending = T,
                                         outname="max_val"){

  tracking <- tracking %>% group_by({{group_var}})

  if(descending){
    tracking <- tracking %>% arrange(desc({{arr_var}}))
  } else tracking <- tracking %>% arrange({{arr_var}})

  tracking <- tracking %>% summarise(!!outname:= {{summ_var}}[1])

  return(tracking)
}



summarize_GL <- function(tracking,
                         ol_cc=overlapping_class_codes,
                         sampnames=sample_names,
                         ordered_classcodes=c("=",
                                              "j",
                                              "k",
                                              "c",
                                              "m",
                                              "n",
                                              "e",
                                              "o",
                                              "p",
                                              "s",
                                              "x",
                                              "i",
                                              "y",
                                              "u",
                                              "r",
                                              ".")){

  samp1=sampnames[1]
  samp2=sampnames[2]

  tracking <- tracking %>% mutate(gene_name = ifelse(V4%in%ol_cc,
                                                     ifelse(!is.na(ref_gene_name_1),
                                                            ref_gene_name_1,
                                                            ref_gene_name_2),
                                                     V2))

  tracking_GL <- tracking %>% group_by(gene_name) %>% summarise(!!samp1:=any(q1_gene!="-"),
                                                                !!samp2:=any(q2_gene!="-"))

  # get maximum values per gene for variables of interest
  variables= c("mean_tpm_1","mean_tpm_2","Nsamps_1","Nsamps_2","Nexons_1","Nexons_2")

  max_vals <- map(variables, function(var) {
    onam=paste0("max_",var)
    var=sym(var)
    get_max_per_group_by_arrange(tracking,
                                 group_var = gene_name,
                                 arr_var = !!var,
                                 outname = onam)
  })

  tracking_GL <- left_join(tracking_GL,Reduce(merge,max_vals))

  # best classcode per sample

  c2 <- tracking %>% mutate(!!samp1:=V5!="-",!!samp2:=V6!="-") %>%
    pivot_longer(cols = c(!!samp1,!!samp2),names_to = "sample",values_to = "is") %>%
    filter(is) %>% group_by(gene_name,sample) %>% arrange(match(V4,ordered_classcodes)) %>%
    summarise(best_cc=V4[1])
  c3 <- c2 %>% pivot_wider(names_from = sample,values_from = best_cc)
  c3[is.na(c3)] = "NA" # try to preserve the sample order
  c3 <- c3 %>% select(gene_name,!!samp1,!!samp2)
  colnames(c3)[2:3]=paste0(colnames(c3)[2:3],"_cc")

  tracking_GL <- left_join(tracking_GL,c3)

  return(tracking_GL)
}


# The following are split by feature bc doing max(feat,na.rm=T) took a lot of time,
# I realized arrange method is better, need to generalize
# mean tpm



test=get_max_per_group_by_arrange(tracking,V2,mean_tpm_1,outname = "maxTPM1",summ_var = ref_biotype_1)
test=get_max_per_group_by_arrange(tracking,gene_name,mean_tpm_1,outname = "maxTPM1")

var="mean_tpm_1"
test=get_max_per_group_by_arrange(tracking,gene_name,mean_tpm_1,outname = "maxTPM1")
var=sym(var)
onam="maxTPM1"
test=get_max_per_group_by_arrange(tracking,gene_name,!!var,outname = onam)

variables= c("mean_tpm_1","mean_tpm_2","Nsamps_1","Nsamps_2","Nexons_1","Nexons_2")
# group by gene_name summarise each sample occurrence
results <- map(variables, function(var) {
  onam=paste0("max_",var)
  var=sym(var)
  get_max_per_group_by_arrange(tracking,
                               group_var = gene_name,
                               arr_var = !!var,
                               outname = onam)
})

Reduce(merge,results)

tragl1 <- tracking %>% group_by(gene_name) %>% arrange(-mean_tpm_1) %>% summarise(maxTPM1=mean_tpm_1[1])
tragl2 <- tracking %>% group_by(gene_name) %>% arrange(-mean_tpm_2) %>% summarise(maxTPM2=mean_tpm_2[1])
GL1=tragl1
GL2=tragl2

# max Nsamps
tragl1 <- tracking %>% group_by(gene_name) %>% arrange(-Nsamps_1) %>% summarise(maxNsamps1=Nsamps_1[1])
tragl2 <- tracking %>% group_by(gene_name) %>% arrange(-Nsamps_2) %>% summarise(maxNsamps2=Nsamps_2[1])
GL1 <- left_join(GL1,tragl1)
GL2 <- left_join(GL2,tragl2)

# max Nexons
tragl1 <- tracking %>% group_by(gene_name) %>% arrange(-Nexons_1) %>% summarise(maxNexons1=Nexons_1[1])
tragl2 <- tracking %>% group_by(gene_name) %>% arrange(-Nexons_2) %>% summarise(maxNexons2=Nexons_2[1])
GL1 <- left_join(GL1,tragl1)
GL2 <- left_join(GL2,tragl2)

# pull everything into a single table
tracking_GL <- left_join(tracking_GL,GL1)
tracking_GL <- left_join(tracking_GL,GL2)

# best classcode per sample
c2 <- tracking %>% mutate(StemLinc=V5!="-",Klimmeck=V6!="-") %>%
  pivot_longer(cols = c(StemLinc,Klimmeck),names_to = "sample",values_to = "is") %>%
  filter(is) %>% group_by(gene_name,sample) %>% arrange(match(V4,ordered_classcodes)) %>%
  summarise(best_cc=V4[1])
c3 <- c2 %>% pivot_wider(names_from = sample,values_from = best_cc)
c3[is.na(c3)] = "NA" # try to preserve the sample order
c3 <- c3 %>% select(gene_name,StemLinc,Klimmeck)
colnames(c3)[2:3]=paste0(colnames(c3)[2:3],"_cc")

tracking_GL <- left_join(tracking_GL,c3)
