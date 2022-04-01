#!/bin/env Rscript
##### script to add the toxins families
library(argparser)

argu<-arg_parser(description = 'mix busco loci with sum of toxin families from several species', name = 'DifExp')

argu<-add_argument(parser = argu, arg = "--input", short = '-i',help = 'input file, should be a csv with the expected counts for each individual from rsem, from all the species',default = 'Cgodm_Ctzot_consense_ExpCount_df.csv')

argu<-add_argument(parser = argu, arg = "--Tox", short = '-t',help = 'list of files, should be a list of csv files with the expected counts of the toxins,',default = 'list')

argu<-add_argument(argu, arg = '--species', short = '-s', help = 'name of the species in the tabla, example Cgodm', default = 'Cgodm_Ctzot')

argmts<- parse_args(argu)

print('loading libraries')

# Read Libraries and Functions.
library(readr)
library(dplyr)
library(tidyverse)

print('libraries loaded')

print('import expected count data')

# Prepare Expected Counts
expcounts <- read_csv(paste(argmts$input))
x<-ncol(expcounts)
tmp_df<-expcounts

a<-colnames(tmp_df[,4:ncol(tmp_df)])
a
print(paste0("table ",argmts$input," imported,it contains the individuals:"))
print(a)

# Read tox data set
print("reading toxins data bases from the list")

toxlist<- scan(argmts$Tox,what = "character",sep = "\n")

for (i in 1:length(toxlist)){
  cmd = paste0("tox",i,"<- read_csv(toxlist[",i,"])")
  eval(parse(text=cmd))
}

print(paste0(length(toxlist)," files imported"))


for (i in 1:length(toxlist)){
  cmd = paste0("expcounts<- tox",i)
  eval(parse(text=cmd))
  expcountsTF<-as.data.frame(expcounts %>% group_by(toxin_family,class) %>%
                               summarize_if(is_numeric,sum))
  expcountsTF <- expcountsTF %>% filter(class=="Toxin")
  cmd = paste0("TF",i,"<- expcountsTF")
  eval(parse(text = cmd))


  
}
toxins<-c()
for (i in 1:length(toxlist)){
  cmd = paste0("toxins <- c(toxins, TF",i,"$toxin_family)")
  eval(parse(text = cmd))
}

toxins<-unique(toxins)
tmpTF<- data.frame(gene_id = toxins, class = rep ("Toxin",length(toxins)), toxin_family = toxins)

#####
print("mixing the data sets of the toxins into one data.frame")

tmplist<-c()
for (i in 1:length(toxlist)){
  cmd = paste0("tmplist<- TF",i)
  eval(parse(text = cmd))
  n<-ncol(tmpTF)
  for (j in 3:(ncol(tmplist)-1)){
    tmpTF<- data.frame(tmpTF,rep(0, length(toxins)))
    colnames(tmpTF)[n+j-2]<-colnames(tmplist)[j]
    for (k in 1:length(toxins)){
      if(length(which(tmplist$toxin_family == toxins[k])) != 0){
        tmpTF[which(tmpTF$toxin_family == toxins[k]),(n+j-2)]<- tmplist[which(tmplist$toxin_family == toxins[k]),j] 
      }
    }
  }
}

z<-c()

print("adding Average colum")

for (i in 1:nrow(tmpTF)){
  x <- tmpTF[i,4:ncol(tmpTF)]
  #  y <- levels(as.factor(x))
  z <- c(z,round(mean(as.numeric(x)),4))
}

tmpTF<-tmpTF %>% mutate(Average = z)

###
print("correcting posible incongruences in the data sets")

colnames(tmpTF)<-gsub(".CLP","-CLP",colnames(tmpTF))
colnames(tmp_df)<-gsub("Cgodm-CLP2903","Ctzot-CLP2903",colnames(tmp_df))

print("merging toxins and busco data sets")

output<-rbind(tmp_df,tmpTF)

print("writing output csv")

write.csv(output,paste0(argmts$species,'_consenseTF_df.csv'), row.names = F)

