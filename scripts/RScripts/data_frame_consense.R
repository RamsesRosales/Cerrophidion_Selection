#!/bin/env Rscript
#this R code will take a rsem data set, extension .genes.results

library(argparser)
library(tidyverse)

p<- arg_parser('Add the colum identifier with the individual identifier for a single file or for a list of files',
               'add_id')
p<- add_argument(p,short = '-n', arg = '--name', default = 'No', 
                 help = 'the name or identifier of the individual that was sequenced')
p<- add_argument(p, short = '-l', arg = '--list', default = 'No',
                 help = 'the with the name of the individuals')

argu<- parse_args(p)


if (argu$name != 'No'){
  print(paste('reading table',argu$name))
  
  rsemdf<-read.table(paste0(argu$name,'.genes.results'),header = T)
  print(paste('adding colum with identifier',argu$name))
  rsemdf<- rsemdf %>% mutate(ID = paste(argu$name))
  T_NT<-c()
  
  for (i in 1:length(rsemdf$TPM)){
    if( length(grep('TOXIN',rsemdf$gene_id[i])) != 0) {
      T_NT<-c(T_NT,'Toxin')
    }
    else{
      T_NT<-c(T_NT,'Nontoxin')
    }
  }
  TPM_df2 <- rsemdf %>% mutate( class = T_NT)
  d <-c("3FTx","AChe","AChE","BPP","CNP","CRISP","CTL","Ficolin","FusedToxin","GluCyc",
        "Goannatyrotoxin","HYAL","KUN","LAAO","MP","MYO","NGF","NUC","PDE","PLA2","PLB",
        "SVMMP","SVMPI","SVMPII","SVMPIII","SVSP","TruncHYAL","VDP4","VEGF",
        "VESP","Vespryn","VF","Waprin")
  
  ToxFam1<-rep('NT',length(TPM_df2$TPM))
  
  for (i in 1:length(TPM_df2$TPM)){
    if (TPM_df2$class[i] == 'Toxin'){
      ToxFam1[i]<- 'uncharacterised'
      for (j in d){
        if( length(grep(paste(j),TPM_df2$gene_id[i])) != 0) {
          ToxFam1[i]<- paste(j)
        }
        else{
        }
      }
    }
    else{
    }
  }
  
  TPM_df2 <- TPM_df2 %>%
    mutate(toxin_family = ToxFam1)
  write.csv(x = TPM_df2,file = paste0(argu$name,'_genes.results.csv'), row.names = F)
}else{
  if(argu$list != 'No'){
    identifiers<- scan(paste(argu$list), character(), quote = '')
    print(paste0('adding colum with identifiers and mergin data frames of ',length(identifiers),' elements in ',argu$list))
    db_list<- list()
    for (i in 1:length(identifiers)){
      new_element<-read.table(paste0(identifiers[i],'.genes.results'), header = T)
      new_element<- new_element %>% mutate(ID = paste(identifiers[i]))
      db_list[[paste(identifiers[i])]]<-new_element
    }
    rsem_df<- bind_rows(db_list)
    T_NT<-c()
    
    for (i in 1:length(rsem_df$TPM)){
      if( length(grep('TOXIN',rsem_df$gene_id[i])) != 0) {
        T_NT<-c(T_NT,'Toxin')
      }
      else{
        T_NT<-c(T_NT,'Nontoxin')
      }
    }
    TPM_df2 <- rsem_df %>% mutate( class = T_NT)
    d <-c("3FTx","AChe","AChE","BPP","CNP","CRISP","CTL","Ficolin","FusedToxin","GluCyc",
          "Goannatyrotoxin","HYAL","KUN","LAAO","MP","MYO","NGF","NUC","PDE","PLA2","PLB",
          "SVMMP","SVMPI","SVMPII","SVMPIII","SVSP","TruncHYAL","VDP4","VEGF",
          "VESP","Vespryn","VF","Waprin")
    
    ToxFam1<-rep('NT',length(TPM_df2$TPM))
    
    for (i in 1:length(TPM_df2$TPM)){
      if (TPM_df2$class[i] == 'Toxin'){
        ToxFam1[i]<- 'uncharacterised'
        for (j in d){
          if( length(grep(paste(j),TPM_df2$gene_id[i])) != 0) {
            ToxFam1[i]<- paste(j)
          }
          else{
          }
        }
      }
      else{
      }
    }
    
    TPM_df2 <- TPM_df2 %>%
      mutate(toxin_family = ToxFam1)
    write.csv(TPM_df2,paste0(argu$list,'_rsem_Cerrophidion.csv'),row.names = F)
  }
  else{
    print('no input data, use -n to add a name of the identifier of a single file <name>.genes.results
              or use -l to reference a file with a list of names')
  }
}

