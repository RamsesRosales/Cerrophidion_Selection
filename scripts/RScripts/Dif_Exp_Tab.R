#!/bin/env Rscript
##### program to do a DEanalysis
library(argparser)

argu<-arg_parser(description = 'Take the output of the program of Dif_Exp_Cerrophidion.R, it must be run inside the directory created by this program', name = 'DifExp')

argu<-add_argument(parser = argu, arg = "--input", short = '-i',help = 'input file, should be a csv with the expected counts for each individual from rsem',default = '../Cgodm_consense_ExpCount_df.csv')

argu<-add_argument(parser = argu, arg = "--Pop", short = '-p',help = 'Set as true if your comparing populations, ser false if you only need SVL',default = T)

argu<-add_argument(argu, arg = '--species', short = '-s', help = 'name of the species in the tabla, example Cgodm, will be used to name the output', default = 'Cgodm')

argu<-add_argument(argu, arg = '--ToxFam', short = '-TF', help = 'set as True if you want to run analysis of Toxin Families too', default = T)

argmts<- parse_args(argu)

print('load libraries')

library(tidyverse)

print('loaded libraries')

print('import expected count table')

# Prepare Expected Counts
expcounts <- read_csv(paste(argmts$input))
#expcounts <- expcounts %>% filter(class=="Toxin") we did not run this as we get different results if we use all the genes and only toxins.
x<-ncol(expcounts)-1
tmp_df<-as.matrix(expcounts[,4:x])
rownames(tmp_df)<-expcounts$gene_id

expcountsTF<-as.data.frame(expcounts %>% group_by(toxin_family,class) %>%
                             summarize_if(is_numeric,sum))
y<-ncol(expcountsTF)-1
expcountsTF <- expcountsTF %>% filter(class=="Toxin")
tmp_df_TF<-as.matrix(expcountsTF[,3:y])
rownames(tmp_df_TF)<-expcountsTF$toxin_family



x<-rownames(tmp_df_TF)
y<-rownames(tmp_df)


system("ls * | grep '.csv' > tmp")
z<-scan('tmp',what = character())
system("rm tmp")
ind<-grep('TF',z)
filesTF<-z[ind]
files<-z[-ind]
ind<-grep('table',files)
files<-files[ind]
files

print('importing files')
print(files)

for (i in 1:length(files)){
  cmd = paste0("DE_",i,"<- read_csv('",files[i],"',col_names = T)")
  eval(parse(text = cmd))
}

print('files imported')

print('function to compare tables')

compare_tab<- function(tab,y){
  w<-tab$X1
  n<-c()
  m<-c()
  m<-c()
  for (i in 1:length(y)){
    n<-c()
    for (j in 1:length(w)){
      if (y[i] == w[j]){
        n<-c(n,1)
      }
      if(y[i] != w[j]){
        n<-c(n,0)
      }
    }
    if(sum(n) == 1){
      m<-c(m,1)
    }
    if(sum(n) == 0)
      m<-c(m,0)
  }
  a<-as.data.frame(m,row.names = y)
}

print(files)

files1<- gsub('.csv','',files)
files1<- gsub('_tab','',files1)


print('start analysing the files:')
print(files1)



tab1<- data.frame(y)
for (i in 1:length(files1)) {
  tab<-get(eval(paste0('DE_',i)))
  colnames(tab)[1]<- 'X1'
  print(colnames(tab))
  cmd = paste0(files1[i],'<-compare_tab(tab,y)')
  print(cmd)
  eval(parse(text = cmd))
  cmd1 = paste0("tab1<- data.frame(tab1,",files1[i],")")
  print(cmd1)
  eval(parse(text = cmd1))
}

print('analisis finalized')


tab1<-tab1[,-1]
colnames(tab1)<-files1

ind<- grep('EDGER',colnames(tab1))
EdgeRn<-colnames(tab1)[ind]

ind<- grep('DESeq',colnames(tab1))
DESeqn<-colnames(tab1)[ind]

print('ading colums for DESeq and EdgeR aggrement')

ind2<- grep('SVL',colnames(tab1))
b<-c()
for (i in 1:nrow(tab1)){
  if (sum(tab1[i,ind2]) == 2){
    b<- c(b,1)
  }
  if (sum(tab1[i,ind2]) <= 1){
    b<-c(b,0)
  }
}

tab1<- mutate(tab1,SVL = b)

if (argmts$Pop){
  ind1<- grep('Pop',colnames(tab1))
  a<-c()
  for (i in 1:nrow(tab1)){
    if (sum(tab1[i,ind1]) == 2){
      a<- c(a,1)
    }
    if (sum(tab1[i,ind1]) <= 1){
      a<-c(a,0)
    }
  }
  tab1<- mutate(tab1,Pop = a)
}

print('write the csv file as:')
print(paste0(argmts$species,'_DifExp_tab.csv'))

write.csv(tab1, paste0(argmts$species,'_DifExp_tab.csv'), col.names = T)


if (argmts$ToxFam){
  print('repeat analisis for Toxin Families')
  for (i in 1:length(filesTF)){
    cmd = paste0("DE_",i,"<- read_csv('",filesTF[i],"',col_names = T)")
    eval(parse(text = cmd))
  }
  
  files1<- gsub('.csv','',filesTF)
  files1<- gsub('_tab','',files1)
  
  print('processing files:')
  print(files1)

  tab1<- data.frame(x)
  for (i in 1:length(files1)) {
    tab<-get(eval(paste0('DE_',i)))
    colnames(tab)[1]<- 'X1'
    print(colnames(tab))
    cmd = paste0(files1[i],'<-compare_tab(tab,x)')
    eval(parse(text = cmd))
    cmd1 = paste0("tab1<- data.frame(tab1,",files1[i],")")
    eval(parse(text = cmd1))
  }

  print('analisis finalized')
  
  tab1<-tab1[,-1]
  colnames(tab1)<-files1
  
  ind<- grep('EDGER',colnames(tab1))
  EdgeRn<-colnames(tab1)[ind]
  
  ind<- grep('DESeq',colnames(tab1))
  DESeqn<-colnames(tab1)[ind]
  
  
  ind2<- grep('SVL',colnames(tab1))
  b<-c()
  for (i in 1:nrow(tab1)){
    if (sum(tab1[i,ind2]) == 2){
      b<- c(b,1)
    }
    if (sum(tab1[i,ind2]) <= 1){
      b<-c(b,0)
    }
  }
  
  tab1<- mutate(tab1,SVL = b)
  
  if (argmts$Pop){
    ind1<- grep('Pop',colnames(tab1))
    a<-c()
    for (i in 1:nrow(tab1)){
      if (sum(tab1[i,ind1]) == 2){
        a<- c(a,1)
      }
      if (sum(tab1[i,ind1]) <= 1){
        a<-c(a,0)
      }
    }
    tab1<- mutate(tab1,Pop = a)
  }
  
  print('write .csv file as:')
  print(paste0(argmts$species,'_DifExp_TF_tab.csv'))
  
  write.csv(tab1, paste0(argmts$species,'_DifExp_TF_tab.csv'), col.names = T)
  
}



