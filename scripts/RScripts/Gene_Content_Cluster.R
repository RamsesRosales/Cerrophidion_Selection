#!/bin/env Rscript
library(argparser)

argu<-arg_parser(description = 'Script to produce a Gene Content Cluster for one species samples', name = 'gene_cluster')

argu<-add_argument(argu,arg = '--input', short = '-i', default = "Ctzot_consense_df.csv", help = 'Data Frame input that contains in the first 3 colums the gene.id, class, and toxin_family, and then in each colum the rsen results of one individual')

argu<-add_argument(argu,arg = '--number', short = '-n', default = 3, help = 'number of individual that where sampled in the table, excluding average')

argu<-add_argument(argu, arg = '--species_name', short = '-s', default = 'Ctzot', help = 'species name that will be added to the tree and the csv output')

arg<-parse_args(argu)

#setwd('~/Desktop/Cerrophidion/Ctzot_1/17rsem/')

a<-read.csv(as.character(arg$input), header = T)
#a<-read.csv('Ctzot_consense_df.csv', header = T)

a<-a[which(a$class == 'Toxin'), c(4:(3 + arg$number))]

#x<-c(c(colnames(a)),
#     c(length(a[which(a$Ctzot.CLP2364 > 0),1]),length(a[which(a$Ctzot.CLP2366 > 0),2]),length(a[which(a$Ctzot.CLP2383 > 0),2])),
#     c(length(a[which(a$Ctzot.CLP2364 == 0),1]),length(a[which(a$Ctzot.CLP2366 == 0),2]),length(a[which(a$Ctzot.CLP2383 == 0),2])))

#x<-matrix(x,ncol = 3,nrow = 3)
#x

a<- t(a)

a[a>0] <-1

#
A<-t(data.frame(a))
A

Acn<-colnames(A)
un<-c()
cer<-c()
for (i in 1:ncol(A)){
  unos<-sum(A[,i])
  ceros<-length(A[,1])-unos
  un<-c(un,unos)
  cer<-c(cer,ceros)
}

GN_df<-data.frame('Individuals' = Acn, 'Num_Genes' = un, 'Miss_Genes' = cer)


library(vegan)

d=vegdist(a,method = 'jaccard')
d


h=hclust(d, method = 'average')

library(ape)

p=as.phylo(h)
p

write.tree(p,paste0(as.character(arg$species_name),'_toxin_content_cluster.tree'))
write.csv(GN_df,paste0(as.character(arg$species_name),'_toxin_content.csv'),row.names = F)

